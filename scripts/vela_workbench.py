"""Local Tk workbench for frozen sd-sim regression inputs."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
from queue import Empty, Queue
import re
import tempfile
from threading import Event, Thread

from vela_workbench_cases import prepare_cases, read_cases, run_case
from vela_workbench_text import CASE_NAMES, ZH


ACCEPTED_POINT = re.compile(rb"(?m)^stage=(\d+) bias=([+\-0-9.eE]+) V accepted\r?\n")


def last_accepted_point(log):
    try:
        with log.open("rb") as stream:
            size = stream.seek(0, os.SEEK_END)
            stream.seek(max(0, size - 8192))
            accepted = ACCEPTED_POINT.findall(stream.read())
    except OSError:
        return None
    if accepted:
        stage, bias = accepted[-1]
        return int(stage) + 1, float(bias)
    return None


def set_window_icon(root):
    from PIL import Image, ImageDraw, ImageTk

    # Draw the overlapping glyphs at a larger size for a crisp 16 px icon.
    background = "#10272B"
    image = Image.new("RGBA", (96, 96))
    draw = ImageDraw.Draw(image)
    draw.rounded_rectangle((3, 3, 92, 92), radius=20, fill=background)

    def stroke(points, color, width):
        draw.line(points, fill=color, width=width, joint="curve")
        radius = width // 2
        for x, y in (points[0], points[-1]):
            draw.ellipse((x - radius, y - radius, x + radius, y + radius), fill=color)

    letter_s = [
        (36, 20), (30, 17), (22, 17), (16, 21), (14, 27), (16, 32),
        (22, 36), (32, 39), (37, 44), (37, 50), (33, 55), (25, 58), (15, 55),
    ]
    letter_d = [
        (34, 25), (45, 25), (53, 27), (60, 33), (64, 41), (64, 53),
        (60, 61), (53, 67), (45, 69), (34, 69), (34, 25),
    ]
    stroke(letter_s, "#5FC4C1", 7)
    stroke([(x + 45, y + 20) for x, y in letter_s], "#70E6C9", 7)
    stroke(letter_d, background, 11)
    stroke(letter_d, "#E3FFFA", 7)
    icons = [
        ImageTk.PhotoImage(image.resize((size, size), Image.Resampling.LANCZOS), master=root)
        for size in (16, 32, 48)
    ]
    root.iconphoto(True, *icons)


def launch(launcher, cases, runs):
    import tkinter as tk
    from tkinter import ttk

    from matplotlib import rcParams
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
    from matplotlib.figure import Figure

    from vela_workbench_results import FIELDS, read_curve, read_field, read_stage_curves

    rcParams["font.sans-serif"] = ["Microsoft YaHei", *rcParams["font.sans-serif"]]
    root = tk.Tk()
    set_window_icon(root)
    root.geometry("1060x700")
    root.minsize(760, 480)
    chinese = False

    def tr(text, *args):
        label = ZH.get(text, text) if chinese else text
        return label.format(*args) if args else label

    def case_name(key):
        case = cases[key]
        return case.name_zh if chinese and case.name_zh else tr(CASE_NAMES.get(key, case.name))

    toolbar = ttk.Frame(root, padding=(12, 10))
    toolbar.pack(fill="x")
    current = ttk.Button(toolbar)
    batch = ttk.Button(toolbar)
    stop = ttk.Button(toolbar, state="disabled")
    folder = ttk.Button(toolbar, state="disabled")
    language = ttk.Button(toolbar)
    for button in (current, batch, stop, folder, language):
        button.pack(side="left", padx=(0, 8))

    pane = ttk.PanedWindow(root, orient="horizontal")
    pane.pack(fill="both", expand=True, padx=12, pady=(0, 12))
    sidebar = ttk.Frame(pane, padding=(0, 0, 8, 0), width=300)
    pane.add(sidebar)
    case_heading = ttk.Label(sidebar)
    case_heading.pack(anchor="w", pady=(0, 6))
    tree = ttk.Treeview(sidebar, columns=("status",), selectmode="extended", height=8)
    tree.column("#0", width=155)
    tree.column("status", width=155)
    case_scroll = ttk.Scrollbar(sidebar, command=tree.yview)
    case_scroll.pack(side="right", fill="y")
    tree.configure(yscrollcommand=case_scroll.set)
    tree.pack(side="left", fill="both", expand=True)
    for key in cases:
        tree.insert("", "end", iid=key)

    tabs = ttk.Notebook(pane)
    pane.add(tabs, weight=1)

    def plot_tab():
        tab = ttk.Frame(tabs, padding=7)
        tabs.add(tab)
        figure = Figure(figsize=(5.7, 4.5), layout="constrained")
        canvas = FigureCanvasTkAgg(figure, master=tab)
        navigation = NavigationToolbar2Tk(canvas, tab, pack_toolbar=False)
        navigation.pack(side="bottom", fill="x")
        return tab, figure, canvas, navigation

    curve_tab, curve_figure, curve_canvas, curve_tools = plot_tab()
    mesh_tab, mesh_figure, mesh_canvas, mesh_tools = plot_tab()
    cloud_tab, cloud_figure, cloud_canvas, cloud_tools = plot_tab()
    cloud_controls = ttk.Frame(cloud_tab)
    cloud_controls.pack(fill="x")
    field_heading = ttk.Label(cloud_controls)
    field_heading.pack(side="left", padx=(0, 8))
    field_picker = ttk.Combobox(cloud_controls, state="readonly", width=26)
    field_picker.pack(side="left")
    for canvas in (curve_canvas, mesh_canvas, cloud_canvas):
        canvas.get_tk_widget().pack(fill="both", expand=True)

    log_tab = ttk.Frame(tabs, padding=7)
    tabs.add(log_tab)
    log_picker = ttk.Combobox(log_tab, state="readonly")
    log_picker.pack(fill="x", pady=(0, 6))
    script_tab = ttk.Frame(tabs, padding=7)
    tabs.add(script_tab)
    script_heading = ttk.Label(script_tab)
    script_heading.pack(anchor="w", pady=(0, 6))

    def text_view(tab):
        frame = ttk.Frame(tab)
        frame.pack(fill="both", expand=True)
        frame.columnconfigure(0, weight=1)
        frame.rowconfigure(0, weight=1)
        view = tk.Text(frame, wrap="none", state="disabled")
        vertical = ttk.Scrollbar(frame, command=view.yview)
        horizontal = ttk.Scrollbar(frame, orient="horizontal", command=view.xview)
        view.configure(yscrollcommand=vertical.set, xscrollcommand=horizontal.set)
        view.grid(row=0, column=0, sticky="nsew")
        vertical.grid(row=0, column=1, sticky="ns")
        horizontal.grid(row=1, column=0, sticky="ew")
        return view

    log_view = text_view(log_tab)
    script_view = text_view(script_tab)

    progress_panel = ttk.Frame(root, padding=(12, 0, 12, 10))
    progress_panel.pack(fill="x")
    progress_bar = ttk.Progressbar(progress_panel, length=150)
    progress_bar.pack(side="left", padx=(0, 10))
    headline = ttk.Label(progress_panel)
    headline.pack(side="left")

    directories = {}
    events = Queue()
    stop_pending = Event()
    active = False
    statuses = dict.fromkeys(cases, ("Ready",))
    total = processed = canceled = 0
    running_key = None
    close_notice = False
    remembered_log = None
    field = None
    field_terminal = None
    field_keys = []
    view_key = None

    def focused():
        return tree.focus() or (tree.selection() or ("",))[0]

    def replace_text(view, content, *, follow=False):
        view.configure(state="normal")
        view.delete("1.0", "end")
        view.insert("1.0", content)
        if follow:
            view.see("end")
        view.configure(state="disabled")

    def log_update():
        nonlocal remembered_log
        directory = directories.get(focused())
        filename = log_picker.get()
        path = directory / filename if directory and filename else None
        if path and path.is_file():
            size = path.stat().st_size
            signature = (path, size, chinese)
            if signature == remembered_log:
                return
            with path.open("rb") as stream:
                offset = max(0, size - 64 * 1024)
                stream.seek(offset)
                content = stream.read().decode("utf-8", errors="replace")
            if offset:
                notice = tr("Showing the end of the log. Open the run folder for the full log.")
                content = notice + "\n\n" + content
        else:
            signature, content = None, ""
            if remembered_log is None:
                return
        remembered_log = signature
        replace_text(log_view, content, follow=True)

    def axes_for(figure, navigation):
        figure.clear()
        navigation.update()
        return figure.add_subplot()

    def empty_plot(axes, message):
        axes.set_axis_off()
        axes.text(.5, .5, message, ha="center", va="center", transform=axes.transAxes)

    def draw_field(_event=None):
        axes = axes_for(cloud_figure, cloud_tools)
        selection = field_picker.current()
        if field is None or selection < 0:
            empty_plot(axes, tr("No final field available"))
        else:
            name = field_keys[selection]
            x, y = field.points[:, :2].T
            image = axes.tripcolor(
                x, y, field.triangles, field.quantities[name], shading="gouraud",
            )
            unit = field.metadata.get(f"Unit:{name}", "")
            cloud_figure.colorbar(image, ax=axes, label=unit)
            axes.set(xlabel="x (µm)", ylabel="y (µm)", aspect="equal")
            stopped = field.metadata.get("GoalReached", field.metadata.get("goal_reached")) == 0
            bias = field.metadata.get(f"{field_terminal}_OuterVoltage_V")
            title = tr("Last accepted point" if stopped else "Final state")
            if bias is not None:
                title += " · " + tr("{} outer voltage {:.6g} V", field_terminal, bias)
            axes.set_title(title)
        cloud_canvas.draw_idle()

    def refresh(_event=None):
        nonlocal field, field_keys, field_terminal, view_key
        key = focused()
        directory = directories.get(key)
        folder.configure(state="normal" if directory else "disabled")
        case = cases[key]
        work = directory / case.work if directory else None
        plan = {}
        if work:
            try:
                plan = json.loads((work / "stage_plan.json").read_text(encoding="utf-8"))
            except (OSError, ValueError):
                pass
        stages = plan.get("stages", [])
        curve_path = directory / case.curve if directory and case.curve else None
        field_name = case.field or plan.get("outputs", {}).get("final_vtu")
        field_path = directory / field_name if directory and field_name else None
        paths = [curve_path, field_path]
        if work and not curve_path:
            paths.extend(work / stage["csv"] for stage in stages if "csv" in stage)
        identity = (
            key, directory, chinese,
            *(path.stat().st_mtime_ns if path and path.is_file() else None for path in paths),
        )
        log_names = []
        if directory:
            log_names = ["sd-sim.log"]
            log_names.extend(str(path.relative_to(directory)) for path in sorted(directory.glob("*_des.log")))
        selected_log = log_picker.get()
        log_picker.configure(values=log_names)
        log_picker.set(selected_log if selected_log in log_names else (log_names or [""])[0])
        log_update()
        if identity == view_key:
            return
        if view_key is None or view_key[0] != key or view_key[2] != chinese:
            script = tr("Original script unavailable. Prepare a new case catalog.")
            if case.original_script:
                try:
                    script = (case.source / case.original_script).read_text(encoding="utf-8")
                except (OSError, UnicodeError) as error:
                    script = tr("Could not read original script: {}", error)
            replace_text(script_view, script)
        view_key = identity
        curve = axes_for(curve_figure, curve_tools)
        mesh = axes_for(mesh_figure, mesh_tools)
        field = None
        field_terminal = None
        try:
            if curve_path:
                points = read_curve(curve_path) if curve_path.is_file() else []
                curves = [(None, "Anode", points)] if points else []
            else:
                curves = read_stage_curves(work, stages) if work else []
            if curves:
                for index, terminal, points in curves:
                    x, y = zip(*points)
                    label = tr("Stage {} · {}", index, terminal) if index is not None else terminal
                    curve.plot(x, y, marker="." if len(points) == 1 else None, label=label)
                curve.set(
                    xlabel=tr("Anode outer voltage (V)" if curve_path else "Swept terminal outer voltage (V)"),
                    ylabel=tr("Anode total current (A/µm)" if curve_path else "Swept terminal current (A/µm)"),
                )
                if not curve_path:
                    curve.legend()
                curve.grid(alpha=.25)
            else:
                empty_plot(curve, tr("No curve available"))
        except (OSError, ValueError, KeyError) as error:
            empty_plot(curve, tr("Could not read curve: {}", error))
        try:
            if field_path and field_path.is_file():
                field = read_field(field_path)
                stage_index = int(field.metadata.get("Stage", -1))
                if 0 <= stage_index < len(stages):
                    field_terminal = (stages[stage_index].get("sweep") or {}).get("contact")
                if field_terminal is None:
                    field_terminal = next(
                        (name.removesuffix("_OuterVoltage_V") for name in field.metadata
                         if name.endswith("_OuterVoltage_V") and name != "Cathode_OuterVoltage_V"),
                        None,
                    )
                x, y = field.points[:, :2].T
                mesh.triplot(x, y, field.triangles, linewidth=.45)
                mesh.set(xlabel="x (µm)", ylabel="y (µm)", aspect="equal")
            else:
                empty_plot(mesh, tr("No mesh available"))
        except (OSError, ValueError) as error:
            empty_plot(mesh, tr("Could not read mesh or field: {}", error))
        index = field_picker.current()
        selected_field = field_keys[index] if index >= 0 else None
        field_keys = list(field.quantities) if field else []
        field_picker.configure(values=[tr(FIELDS[name]) for name in field_keys])
        if field_keys:
            index = field_keys.index(selected_field) if selected_field in field_keys else 0
            field_picker.current(index)
        else:
            field_picker.set("")
        curve_canvas.draw_idle()
        mesh_canvas.draw_idle()
        draw_field()

    def work(keys):
        for key in keys:
            if stop_pending.is_set():
                events.put((key, None, ("Canceled in queue",)))
                continue
            directory = None
            try:
                directory = Path(tempfile.mkdtemp(prefix=f"{key}-", dir=runs))
                events.put((key, directory, ("Running",)))
                status = run_case(launcher, cases[key], directory)
            except (OSError, ValueError, KeyError) as error:
                status = ("Run failed: {}", error)
            events.put((key, directory, status))
        events.put((None, None, ()))

    def update_case(key, status):
        statuses[key] = status
        tree.item(key, text=case_name(key), values=(tr(*status),))

    def progress_text():
        if not total:
            return tr("Use Ctrl or Shift to select multiple cases")
        count = tr("{}/{} cases processed", processed, total)
        if not active:
            message = tr("Queue stopped" if canceled else "Batch finished")
        elif close_notice:
            message = tr("Close the window after the current run finishes")
        elif stop_pending.is_set():
            message = tr("Finishing the current case before stopping the queue")
        elif running_key:
            message = tr("Running {}", case_name(running_key))
            latest = last_accepted_point(directories[running_key] / "sd-sim.log")
            if latest:
                message += " · " + tr("Stage {} · last accepted bias {:.6g} V", *latest)
        else:
            message = tr("Starting solver")
        return " · ".join((count, message))

    def apply_language():
        root.title(tr("Semiconductor Device Simulator"))
        for button, label in (
            (current, "Run Current"), (batch, "Run Selected"),
            (stop, "Stop Queue"), (folder, "Open Run Folder"),
        ):
            button.configure(text=tr(label))
        language.configure(text="English" if chinese else "中文")
        case_heading.configure(text=tr("Built-in cases"))
        tree.heading("#0", text=tr("Case"))
        tree.heading("status", text=tr("Status"))
        for key, status in statuses.items():
            update_case(key, status)
        for tab, label in (
            (curve_tab, "Curve"), (mesh_tab, "Mesh"), (cloud_tab, "Final Field"),
            (log_tab, "Log"), (script_tab, "Sentaurus Script"),
        ):
            tabs.tab(tab, text=tr(label))
        field_heading.configure(text=tr("Field"))
        script_heading.configure(text=tr("Original Sentaurus input · read only"))
        headline.configure(text=progress_text())
        refresh()

    def switch_language():
        nonlocal chinese
        chinese = not chinese
        apply_language()

    def start(keys):
        nonlocal active, total, processed, canceled, close_notice
        if active or not keys:
            return
        active = True
        total = len(keys)
        processed = canceled = 0
        close_notice = False
        stop_pending.clear()
        for button in (current, batch):
            button.configure(state="disabled")
        stop.configure(state="normal")
        for key in keys:
            update_case(key, ("Queued",))
        progress_bar.configure(
            mode="indeterminate" if total == 1 else "determinate", maximum=total, value=0,
        )
        if total == 1:
            progress_bar.configure(maximum=100)
            progress_bar.start()
        headline.configure(text=progress_text())
        Thread(target=work, args=(keys,), daemon=True).start()

    def poll():
        nonlocal active, processed, canceled, running_key
        changed = False
        try:
            while True:
                key, directory, status = events.get_nowait()
                if key is None:
                    active = False
                    for button in (current, batch):
                        button.configure(state="normal")
                    stop.configure(state="disabled")
                    progress_bar.stop()
                    progress_bar.configure(mode="determinate", maximum=total, value=processed)
                else:
                    update_case(key, status)
                    if directory:
                        directories[key] = directory
                    if status[0] == "Running":
                        running_key = key
                    else:
                        running_key = None
                        processed += 1
                        canceled += status[0] == "Canceled in queue"
                        if total > 1:
                            progress_bar.configure(value=processed)
                    changed = True
        except Empty:
            pass
        if changed:
            refresh()
        else:
            log_update()
        headline.configure(text=progress_text())
        root.after(250, poll)

    def close():
        nonlocal close_notice
        if active:
            close_notice = True
            headline.configure(text=progress_text())
        else:
            root.destroy()

    current.configure(command=lambda: start((focused(),)))
    batch.configure(command=lambda: start(tree.selection()))
    stop.configure(command=stop_pending.set)
    folder.configure(command=lambda: os.startfile(directories[focused()]))
    language.configure(command=switch_language)
    root.protocol("WM_DELETE_WINDOW", close)
    tree.bind("<<TreeviewSelect>>", refresh)
    field_picker.bind("<<ComboboxSelected>>", draw_field)
    log_picker.bind("<<ComboboxSelected>>", lambda _: log_update())
    first = next(iter(cases))
    tree.focus(first)
    tree.selection_set(first)
    apply_language()
    root.after(250, poll)
    root.mainloop()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage-cases", type=Path, help="Prepare a new frozen case catalog and exit")
    parser.add_argument("--cases", type=Path, default=Path("cases"))
    parser.add_argument("--launcher", type=Path, default=Path("vela/sd-sim.exe"))
    data = Path(os.environ.get("LOCALAPPDATA", Path.home()))
    parser.add_argument("--runs", type=Path, default=data / "VelaTCAD" / "runs")
    args = parser.parse_args()
    if args.stage_cases:
        prepare_cases(args.stage_cases.resolve())
        return
    launcher = args.launcher.resolve()
    if not launcher.is_file():
        parser.error(f"sd-sim launcher does not exist: {launcher}")
    cases = read_cases(args.cases.resolve())
    if not cases:
        parser.error("No workbench cases are available")
    runs = args.runs.resolve()
    runs.mkdir(parents=True, exist_ok=True)
    launch(launcher, cases, runs)


if __name__ == "__main__":
    main()
