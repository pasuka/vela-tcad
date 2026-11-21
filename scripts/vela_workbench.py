"""Local Tk workbench for frozen sd-sim regression inputs."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
from queue import Empty, Queue
import subprocess
import tempfile
from threading import Event, Thread
from time import perf_counter

from vela_workbench_cases import prepare_cases, read_cases, read_run_status, run_case, verify_inputs
from vela_workbench_controls import action_icon, add_tooltip
from vela_workbench_settings import Settings, load_settings, save_language
from vela_workbench_text import CASE_NAMES, ZH


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


def launch(launcher, cases, runs, saved_run=None, *, config_path=None):
    import tkinter as tk
    from tkinter import filedialog, font, messagebox, ttk

    settings = Settings()
    settings_error = None
    if config_path:
        try:
            settings = load_settings(config_path)
        except (OSError, ValueError) as error:
            settings_error = error

    root = tk.Tk()
    set_window_icon(root)
    root.geometry("1150x720")
    root.minsize(650 + settings.font_size * 27, 500)
    chinese = settings.language == "zh-CN"
    style = ttk.Style(root)
    style.theme_use("clam")
    background, surface, ink = "#EDF1F5", "#FFFFFF", "#182A3D"
    muted, line, accent = "#4D6074", "#D9E2EB", "#245A9C"
    selection_fill, selection_ink = "#DCE9F8", "#193D69"
    root.configure(background=background)
    for name in ("TkDefaultFont", "TkTextFont", "TkMenuFont", "TkHeadingFont"):
        font.nametofont(name).configure(size=settings.font_size)
    mono = font.nametofont("TkFixedFont")
    mono.configure(size=settings.font_size)
    heading_font = font.nametofont("TkHeadingFont").copy()
    heading_font.configure(weight="bold")
    style.configure(".", background=background, foreground=ink, font="TkDefaultFont")
    style.configure("Surface.TFrame", background=surface)
    style.configure("Surface.TLabel", background=surface)
    style.configure("Heading.Surface.TLabel", font=heading_font)
    style.configure("Muted.TLabel", foreground=muted)
    style.configure("Tooltip.TLabel", background=selection_ink, foreground=surface, padding=(8, 5))
    style.configure("TButton", width=0, padding=(10, 6), background=surface, bordercolor=line)
    style.map("TButton", background=[("disabled", background), ("active", "#E6EEF7")],
              foreground=[("disabled", "#748397")])
    style.configure("Accent.TButton", background=accent, foreground=surface, bordercolor=accent)
    style.map("Accent.TButton",
              background=[("disabled", background), ("pressed", "#183F73"), ("active", "#1F4D87")],
              foreground=[("disabled", "#748397"), ("!disabled", surface)])
    style.configure("Treeview", rowheight=font.nametofont("TkDefaultFont").metrics("linespace") + 12,
                    background=surface, fieldbackground=surface, borderwidth=0)
    style.configure("Treeview.Heading", background="#E8EEF5", font=heading_font, padding=(6, 7))
    style.map("Treeview", background=[("selected", selection_fill)], foreground=[("selected", selection_ink)])
    style.configure("TNotebook", borderwidth=0, tabmargins=(0, 0, 0, 0))
    style.configure("TNotebook.Tab", padding=(9, 8), foreground=muted)
    style.map("TNotebook.Tab", background=[("selected", surface)], foreground=[("selected", ink)])
    style.configure("Horizontal.TProgressbar", background=accent, troughcolor=line, bordercolor=line)
    timing_path = os.environ.get("VELA_WORKBENCH_TIMINGS")

    def record_timing(event):
        if timing_path:
            with open(timing_path, "a", encoding="utf-8") as stream:
                stream.write(json.dumps({"event": event, "seconds": perf_counter()}) + "\n")

    def draw_plot(canvas):
        if timing_path:
            def after_draw(_event):
                canvas.mpl_disconnect(callback)
                root.after_idle(record_timing, "plot_ready")

            callback = canvas.mpl_connect("draw_event", after_draw)
        canvas.draw_idle()

    def tr(text, *args):
        label = ZH.get(text, text) if chinese else text
        return label.format(*args) if args else label

    def case_name(key):
        case = cases[key]
        return case.name_zh if chinese and case.name_zh else tr(CASE_NAMES.get(key, case.name))

    menubar = tk.Menu(root, tearoff=False)
    root.configure(menu=menubar)
    menus = {}
    for label in ("File", "Simulation", "Settings"):
        menu = tk.Menu(menubar, tearoff=False)
        menubar.add_cascade(menu=menu)
        menus[label] = menu
    actions = {}
    toolbar = ttk.Frame(root, padding=(16, 10))
    toolbar.pack(fill="x")

    def action_label(label):
        if label == "Language":
            return tr("Switch to English") if chinese else "Switch to Chinese"
        return tr(label)

    def add_action(group, label, icon, command):
        primary = icon == "run"
        size = max(18, settings.font_size + 6)
        normal = action_icon(root, icon, size, surface if primary else selection_ink)
        disabled = action_icon(root, icon, size, "#748397")
        button = ttk.Button(toolbar, image=(normal, "disabled", disabled), command=command,
                            style="Accent.TButton" if primary else "TButton")
        button.pack(side="left", padx=(0, 5))
        add_tooltip(button, lambda: action_label(label))
        menu = menus[group]
        menu.add_command(command=command)
        actions[label] = (menu, menu.index("end"), button, normal, disabled)

    pane = ttk.PanedWindow(root, orient="horizontal")
    pane.pack(fill="both", expand=True, padx=16, pady=(0, 8))
    sidebar = ttk.Frame(pane, padding=(10, 9), style="Surface.TFrame")
    pane.add(sidebar)
    case_heading = ttk.Label(sidebar, style="Heading.Surface.TLabel")
    case_heading.pack(anchor="w", pady=(1, 8))
    case_list = ttk.Frame(sidebar, style="Surface.TFrame")
    case_list.pack(fill="both", expand=True)
    tree = ttk.Treeview(case_list, columns=("status",), selectmode="extended", height=8)
    tree.column("#0", width=min(18 * settings.font_size, 230), minwidth=110)
    tree.column("status", width=min(13 * settings.font_size, 155), minwidth=65)
    case_scroll = ttk.Scrollbar(case_list, command=tree.yview)
    case_scroll.pack(side="right", fill="y")
    tree.configure(yscrollcommand=case_scroll.set)
    tree.pack(side="left", fill="both", expand=True)
    tree.tag_configure("alternate", background="#F4F7FB")
    for index, key in enumerate(cases):
        tree.insert("", "end", iid=key, tags=("alternate",) if index % 2 else ())

    results = ttk.Frame(pane, padding=(10, 0), style="Surface.TFrame")
    pane.add(results, weight=1)
    result_heading = ttk.Label(results, style="Heading.Surface.TLabel")
    result_heading.pack(fill="x", padx=(4, 0), pady=(10, 8))
    tabs = ttk.Notebook(results)
    tabs.pack(fill="both", expand=True)

    def plot_tab():
        tab = ttk.Frame(tabs, padding=(9, 10), style="Surface.TFrame")
        tabs.add(tab)
        return tab

    curve_tab, mesh_tab, cloud_tab = (plot_tab() for _ in range(3))
    plots = {}
    drawn = {}

    def plot_widgets(tab):
        if tab in plots:
            return plots[tab]
        data = Path(os.environ.get("LOCALAPPDATA", Path.home()))
        os.environ.setdefault("MPLCONFIGDIR", str(data / "VelaTCAD" / "matplotlib"))
        from matplotlib import rcParams
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
        from matplotlib.figure import Figure

        if not plots:
            rcParams["font.sans-serif"] = [
                "Microsoft YaHei", "Noto Sans CJK SC", "Noto Sans CJK JP", "AR PL UMing CN",
                *rcParams["font.sans-serif"],
            ]
            rcParams["font.size"] = settings.font_size
        figure = Figure(figsize=(5.7, 4.5), layout="constrained")
        canvas = FigureCanvasTkAgg(figure, master=tab)
        navigation = NavigationToolbar2Tk(canvas, tab, pack_toolbar=False)
        plots[tab] = figure, canvas, navigation
        return plots[tab]

    curve_controls = ttk.Frame(curve_tab, style="Surface.TFrame")
    curve_controls.pack(fill="x", pady=(0, 5))
    voltage_row = ttk.Frame(curve_controls, style="Surface.TFrame")
    voltage_row.pack(fill="x")
    voltage_heading = ttk.Label(voltage_row, style="Surface.TLabel")
    voltage_heading.pack(side="left", padx=(0, 8))
    voltage_picker = ttk.Combobox(voltage_row, state="readonly", width=16)
    voltage_picker.pack(side="left")
    voltage_notice = ttk.Label(curve_controls, style="Surface.TLabel", foreground=muted)
    voltage_notice.bind("<Configure>", lambda event: voltage_notice.configure(wraplength=event.width))

    mesh_controls = ttk.Frame(mesh_tab, style="Surface.TFrame")
    mesh_controls.pack(fill="x", pady=(0, 5))
    rotate_left = ttk.Button(mesh_controls, state="disabled")
    rotate_right = ttk.Button(mesh_controls, state="disabled")
    reset_rotation = ttk.Button(mesh_controls, state="disabled")
    for button in (rotate_left, rotate_right, reset_rotation):
        button.pack(side="left", padx=(0, 6))

    cloud_controls = ttk.Frame(cloud_tab, style="Surface.TFrame")
    cloud_controls.pack(fill="x")
    field_heading = ttk.Label(cloud_controls, style="Surface.TLabel")
    field_heading.pack(side="left", padx=(0, 8))
    field_picker = ttk.Combobox(cloud_controls, state="readonly", width=26)
    field_picker.pack(side="left")
    placeholders = {}
    for tab, message in (
        (curve_tab, "No curve available"), (mesh_tab, "No mesh available"),
        (cloud_tab, "No final field available"),
    ):
        label = ttk.Label(tab, anchor="center", style="Surface.TLabel")
        label.pack(expand=True, fill="both")
        placeholders[tab] = (label, message)

    log_tab = ttk.Frame(tabs, padding=(9, 10), style="Surface.TFrame")
    tabs.add(log_tab)
    log_heading = ttk.Label(log_tab, style="Surface.TLabel")
    log_heading.pack(anchor="w", pady=(0, 6))
    script_tab = ttk.Frame(tabs, padding=(9, 10), style="Surface.TFrame")
    tabs.add(script_tab)
    script_heading = ttk.Label(script_tab, style="Surface.TLabel")
    script_heading.pack(anchor="w", pady=(0, 6))

    def text_view(tab):
        frame = ttk.Frame(tab)
        frame.pack(fill="both", expand=True)
        frame.columnconfigure(0, weight=1)
        frame.rowconfigure(0, weight=1)
        view = tk.Text(frame, width=1, height=1, wrap="none", state="disabled", font=mono, background="#F6F8FB",
                       foreground=ink, selectbackground=selection_fill, relief="flat", padx=9, pady=8)
        vertical = ttk.Scrollbar(frame, command=view.yview)
        horizontal = ttk.Scrollbar(frame, orient="horizontal", command=view.xview)
        view.configure(yscrollcommand=vertical.set, xscrollcommand=horizontal.set)
        view.grid(row=0, column=0, sticky="nsew")
        vertical.grid(row=0, column=1, sticky="ns")
        horizontal.grid(row=1, column=0, sticky="ew")
        return view

    log_view = text_view(log_tab)
    script_view = text_view(script_tab)

    progress_panel = ttk.Frame(root, padding=(16, 4, 16, 12))
    progress_panel.pack(side="bottom", fill="x", before=pane)
    progress_panel.columnconfigure(0, weight=1)
    headline = ttk.Label(progress_panel, style="Muted.TLabel")
    headline.grid(row=0, column=0, sticky="ew", padx=(0, 16))
    headline.bind("<Configure>", lambda event: headline.configure(wraplength=event.width))
    progress_bar = ttk.Progressbar(progress_panel, length=165)
    progress_bar.grid(row=0, column=1, sticky="e")
    config_notice = ttk.Label(progress_panel, style="Muted.TLabel")
    if settings_error:
        config_notice.grid(row=1, column=0, columnspan=2, sticky="ew", pady=(4, 0))
        config_notice.bind("<Configure>", lambda event: config_notice.configure(wraplength=event.width))

    directories = {}
    events = Queue()
    stop_pending = Event()
    active = False
    statuses = dict.fromkeys(cases, ("Ready",))
    total = processed = canceled = 0
    running_key = None
    latest_progress = None
    close_notice = False
    remembered_log = None
    field = None
    field_terminal = None
    field_keys = []
    field_revision = None
    field_error = None
    view_key = None
    curve_axes = {}
    mesh_angles = {}

    def focused():
        return tree.focus() or (tree.selection() or ("",))[0]

    def sync_actions():
        states = (
            ("Run Current", not active),
            ("Run Selected", not active and bool(tree.selection())),
            ("Stop Queue", active and not stop_pending.is_set()),
            ("Open Saved Run", not active),
            ("Open Run Folder", focused() in directories),
            ("Open Settings File", config_path is not None),
        )
        for label, enabled in states:
            menu, index, button, *_ = actions[label]
            state = "normal" if enabled else "disabled"
            button.configure(state=state)
            menu.entryconfigure(index, state=state)

    def replace_text(view, content, *, follow=False):
        view.configure(state="normal")
        view.delete("1.0", "end")
        view.insert("1.0", content)
        if follow:
            view.see("end")
        view.configure(state="disabled")

    def log_update():
        nonlocal remembered_log
        if tabs.select() != str(log_tab):
            return
        key = focused()
        directory = directories.get(key)
        path = directory / cases[key].log if directory else None
        signature, content = None, ""
        if path is not None:
            try:
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
            except OSError:
                signature, content = None, ""
        if signature == remembered_log:
            return
        remembered_log = signature
        replace_text(log_view, content, follow=True)

    def axes_for(figure, navigation):
        figure.clear()
        navigation.update()
        return figure.add_subplot()

    def plot_message(tab, message=None):
        label, default = placeholders[tab]
        label.configure(text=message or tr(default))
        if tab in plots:
            _, canvas, navigation = plots[tab]
            canvas.get_tk_widget().pack_forget()
            navigation.pack_forget()
        label.pack(expand=True, fill="both")

    def show_plot(tab):
        placeholders[tab][0].pack_forget()
        figure, canvas, navigation = plot_widgets(tab)
        navigation.pack(side="bottom", fill="x")
        canvas.get_tk_widget().pack(fill="both", expand=True)
        return figure, canvas, axes_for(figure, navigation)

    def display_geometry(axes):
        from vela_workbench_results import display_coordinates

        turns = mesh_angles.get(directories.get(focused()), 0)
        x_label, y_label = ("x′ (µm)", "y′ (µm)") if turns else ("x (µm)", "y (µm)")
        axes.set(xlabel=x_label, ylabel=y_label, aspect="equal")
        return display_coordinates(field.points, turns)

    def draw_field(_event=None):
        selection = field_picker.current()
        if field is None or selection < 0:
            plot_message(cloud_tab)
            return
        cloud_figure, cloud_canvas, axes = show_plot(cloud_tab)
        name = field_keys[selection]
        x, y = display_geometry(axes)
        image = axes.tripcolor(
            x, y, field.triangles, field.quantities[name], shading="gouraud",
        )
        unit = field.metadata.get(f"Unit:{name}", "")
        cloud_figure.colorbar(image, ax=axes, label=unit)
        stopped = field.metadata.get("GoalReached", field.metadata.get("goal_reached")) == 0
        bias = field.metadata.get(f"{field_terminal}_OuterVoltage_V")
        title = tr("Last accepted point" if stopped else "Final state")
        if bias is not None:
            title += " · " + tr("{} outer voltage {:.6g} V", field_terminal, bias)
        axes.set_title(title)
        draw_plot(cloud_canvas)

    def refresh(_event=None):
        nonlocal view_key
        key = focused()
        directory = directories.get(key)
        sync_actions()
        voltage_picker.current(curve_axes.get(key, 0))
        case = cases[key]
        result_heading.configure(text=case_name(key))
        log_heading.configure(text=case.log if directory else "")
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
        curve_files = [curve_path] if curve_path else [
            work / stage["csv"] for stage in stages if work and "csv" in stage
        ]
        field_stamp = file_stamp(field_path)
        identity = (
            key, directory, chinese, field_stamp, *(file_stamp(path) for path in curve_files),
        )
        log_update()
        if view_key is not None and view_key[0] != key:
            drawn.clear()
        if view_key is None or view_key[0] != key or view_key[2] != chinese:
            script = tr("Original script unavailable. Prepare a new case catalog.")
            if case.original_script:
                try:
                    script = (case.source / case.original_script).read_text(encoding="utf-8")
                except (OSError, UnicodeError) as error:
                    script = tr("Could not read original script: {}", error)
            replace_text(script_view, script)
        view_key = identity
        selected = tabs.select()
        tab = tabs.nametowidget(selected) if selected else None
        if tab not in placeholders or drawn.get(tab) == identity:
            return
        drawn[tab] = identity
        if tab == curve_tab:
            draw_curve(curve_path, work, stages, curve_files)
        else:
            draw_mesh_or_field(tab, field_path, field_stamp, stages)

    def file_stamp(path):
        try:
            return path.stat().st_mtime_ns if path else None
        except OSError:
            return None

    def draw_curve(curve_path, work, stages, curve_files):
        voltage_notice.pack_forget()
        if not any(path.is_file() for path in curve_files):
            plot_message(curve_tab)
            return
        from vela_workbench_results import read_curve, read_stage_curves

        try:
            voltage = ("inner", "outer")[curve_axes.get(focused(), 0)]
            if curve_path:
                points = read_curve(curve_path, voltage) if curve_path.is_file() else []
                curves = [(None, "Anode", points)] if points else []
            else:
                curves = read_stage_curves(work, stages, voltage) if work else []
            missing_inner = voltage == "inner" and any(x is None for _, _, points in curves for x, _ in points)
            visible = [
                (index, terminal, points) for index, terminal, points in curves
                if any(x is not None for x, _ in points)
            ]
            if not visible:
                message = tr("Inner voltage unavailable. Select outer voltage.") if missing_inner else None
                plot_message(curve_tab, message)
                return
            if missing_inner:
                voltage_notice.configure(text=tr("Inner voltage unavailable for some points."))
                voltage_notice.pack(fill="x", pady=(4, 0))
            _, curve_canvas, curve = show_plot(curve_tab)
            for index, terminal, points in visible:
                x, y = zip(*points)
                label = tr("Stage {} · {}", index, terminal) if index is not None else terminal
                curve.plot([value if value is not None else float("nan") for value in x], y,
                           marker="." if sum(value is not None for value in x) == 1 else None, label=label)
            terminals = {terminal.casefold() for _, terminal, _ in visible}
            terminal = {"anode": "Anode", "drain": "Drain"}.get(next(iter(terminals)), "Swept terminal")
            if len(terminals) > 1:
                terminal = "Swept terminal"
            curve.set(
                xlabel=tr("{} inner voltage (V)" if voltage == "inner" else "{} outer voltage (V)", tr(terminal)),
                ylabel=tr("Anode total current (A/µm)" if curve_path else "Swept terminal current (A/µm)"),
            )
            if not curve_path:
                curve.legend()
            curve.grid(alpha=.25)
            draw_plot(curve_canvas)
        except (OSError, ValueError, KeyError) as error:
            plot_message(curve_tab, tr("Could not read curve: {}", error))

    def draw_mesh_or_field(tab, field_path, field_stamp, stages):
        nonlocal field, field_keys, field_terminal, field_revision, field_error
        revision = (field_path, field_stamp) if field_stamp is not None else None
        if revision != field_revision:
            field_revision = revision
            field = field_terminal = field_error = None
            if revision:
                from vela_workbench_results import read_field

                try:
                    field = read_field(field_path)
                except (OSError, ValueError) as error:
                    field_error = error
            if field:
                stage_index = int(field.metadata.get("Stage", -1))
                if 0 <= stage_index < len(stages):
                    field_terminal = (stages[stage_index].get("sweep") or {}).get("contact")
                if field_terminal is None:
                    field_terminal = next(
                        (name.removesuffix("_OuterVoltage_V") for name in field.metadata
                         if name.endswith("_OuterVoltage_V") and name != "Cathode_OuterVoltage_V"),
                        None,
                    )
        if tab == mesh_tab:
            for button in (rotate_left, rotate_right, reset_rotation):
                button.configure(state="normal" if field else "disabled")
            if field:
                _, mesh_canvas, mesh = show_plot(mesh_tab)
                x, y = display_geometry(mesh)
                mesh.triplot(x, y, field.triangles, linewidth=.45)
                draw_plot(mesh_canvas)
            else:
                message = tr("Could not read mesh or field: {}", field_error) if field_error else None
                plot_message(mesh_tab, message)
            return
        index = field_picker.current()
        selected_field = field_keys[index] if index >= 0 else None
        field_keys = list(field.quantities) if field else []
        if field_keys:
            from vela_workbench_results import FIELDS

            field_picker.configure(values=[tr(FIELDS[name]) for name in field_keys])
            index = field_keys.index(selected_field) if selected_field in field_keys else 0
            field_picker.current(index)
        else:
            field_picker.configure(values=())
            field_picker.set("")
        if field_error:
            plot_message(cloud_tab, tr("Could not read mesh or field: {}", field_error))
        else:
            draw_field()

    def select_voltage(_event):
        curve_axes[focused()] = voltage_picker.current()
        drawn.pop(curve_tab, None)
        refresh()

    def rotate_mesh(step=None):
        directory = directories.get(focused())
        if directory is None or field is None:
            return
        mesh_angles[directory] = 0 if step is None else (mesh_angles.get(directory, 0) + step) % 4
        drawn.pop(mesh_tab, None)
        drawn.pop(cloud_tab, None)
        refresh()

    def open_saved():
        key = focused()
        selected = filedialog.askdirectory(parent=root, initialdir=runs, title=tr("Open Saved Run"))
        if selected:
            load_saved(key, Path(selected))

    def load_saved(key, directory):
        record_timing("open_saved_run")
        case = cases[key]
        try:
            status = read_run_status(case, directory)
            verify_inputs(case, directory)
        except (OSError, ValueError, KeyError) as error:
            messagebox.showerror(tr("Open Saved Run"), tr("Could not open run for {}: {}", case_name(key), error),
                                 parent=root)
            return
        directories[key] = directory
        update_case(key, status)
        refresh()

    def work(keys):
        for key in keys:
            if stop_pending.is_set():
                events.put((key, None, ("Canceled in queue",)))
                continue
            directory = None
            try:
                directory = Path(tempfile.mkdtemp(prefix=f"{key}-", dir=runs))
                events.put((key, directory, ("Running",)))
                def on_progress(stage, bias):
                    events.put((key, directory, ("progress", stage, bias)))

                status = run_case(launcher, cases[key], directory, log_level=settings.log_level,
                                  on_progress=on_progress)
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
            if latest_progress is not None:
                message += " · " + tr("Stage {} · last accepted bias {:.6g} V", *latest_progress)
        else:
            message = tr("Starting solver")
        return " · ".join((count, message))

    def apply_language():
        root.title(tr("Semiconductor Device Simulator"))
        if settings_error:
            config_notice.configure(text=tr("Using defaults. Check the settings file: {}", settings_error))
        for index, label in enumerate(menus):
            menubar.entryconfigure(index, label=tr(label))
        for label, (menu, index, button, *_) in actions.items():
            text = action_label(label)
            menu.entryconfigure(index, label=text)
            button.configure(text=text)
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
        voltage_heading.configure(text=tr("Voltage axis"))
        voltage_picker.configure(values=(tr("Inner voltage"), tr("Outer voltage")))
        rotate_left.configure(text=tr("Rotate Left 90°"))
        rotate_right.configure(text=tr("Rotate Right 90°"))
        reset_rotation.configure(text=tr("Reset Rotation"))
        field_heading.configure(text=tr("Field"))
        script_heading.configure(text=tr("Original Sentaurus input · read only"))
        headline.configure(text=progress_text())
        refresh()

    def switch_language():
        nonlocal chinese
        chinese = not chinese
        apply_language()
        if config_path:
            try:
                save_language(config_path, "zh-CN" if chinese else "en")
            except (OSError, ValueError) as error:
                messagebox.showwarning(tr("Settings"), tr("Could not save the language: {}", error), parent=root)

    def open_settings():
        try:
            if os.name == "nt":
                subprocess.Popen(["notepad.exe", str(config_path)])
            else:
                subprocess.Popen(["xdg-open", str(config_path)])
        except OSError as error:
            messagebox.showerror(tr("Settings"), tr("Could not open the settings file: {}", error), parent=root)

    def open_run_folder():
        directory = directories[focused()]
        try:
            if os.name == "nt":
                os.startfile(directory)
            else:
                subprocess.Popen(["xdg-open", str(directory)])
        except OSError as error:
            messagebox.showerror(tr("Open Run Folder"), tr("Could not open the run folder: {}", error), parent=root)

    def start(keys):
        nonlocal active, total, processed, canceled, close_notice
        if active or not keys:
            return
        active = True
        total = len(keys)
        processed = canceled = 0
        close_notice = False
        stop_pending.clear()
        sync_actions()
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

    def stop_queue():
        stop_pending.set()
        sync_actions()
        headline.configure(text=progress_text())

    def poll():
        nonlocal active, processed, canceled, running_key, latest_progress
        changed = False
        try:
            while True:
                key, directory, status = events.get_nowait()
                if key is not None and status[0] == "progress":
                    if key == running_key and directory == directories[key]:
                        latest_progress = status[1:]
                    continue
                if key is None:
                    active = False
                    sync_actions()
                    progress_bar.stop()
                    progress_bar.configure(mode="determinate", maximum=total, value=processed)
                else:
                    update_case(key, status)
                    if directory:
                        directories[key] = directory
                    if status[0] == "Running":
                        running_key = key
                        latest_progress = None
                    else:
                        running_key = None
                        latest_progress = None
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

    add_action("Simulation", "Run Current", "run", lambda: start((focused(),)))
    add_action("Simulation", "Run Selected", "batch", lambda: start(tree.selection()))
    add_action("Simulation", "Stop Queue", "stop", stop_queue)
    ttk.Separator(toolbar, orient="vertical").pack(side="left", fill="y", padx=(1, 7), pady=3)
    add_action("File", "Open Saved Run", "saved", open_saved)
    add_action("File", "Open Run Folder", "folder", open_run_folder)
    ttk.Separator(toolbar, orient="vertical").pack(side="left", fill="y", padx=(1, 7), pady=3)
    add_action("Settings", "Open Settings File", "settings", open_settings)
    add_action("Settings", "Language", "language", switch_language)
    rotate_left.configure(command=lambda: rotate_mesh(1))
    rotate_right.configure(command=lambda: rotate_mesh(-1))
    reset_rotation.configure(command=rotate_mesh)
    root.protocol("WM_DELETE_WINDOW", close)
    tree.bind("<<TreeviewSelect>>", refresh)
    tabs.bind("<<NotebookTabChanged>>", refresh)
    voltage_picker.bind("<<ComboboxSelected>>", select_voltage)
    field_picker.bind("<<ComboboxSelected>>", draw_field)
    first = saved_run[0] if saved_run else next(iter(cases))
    tree.focus(first)
    tree.selection_set(first)
    apply_language()
    if timing_path:
        root.after_idle(record_timing, "window_ready")
    if saved_run:
        root.after_idle(lambda: root.after(0, load_saved, *saved_run))
    root.after(250, poll)
    root.mainloop()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage-cases", type=Path, help="Prepare a new frozen case catalog and exit")
    parser.add_argument("--cases", type=Path, default=Path("cases"))
    launcher_name = "sd-sim.exe" if os.name == "nt" else "sd-sim"
    parser.add_argument("--launcher", type=Path, default=Path("vela") / launcher_name)
    parser.add_argument("--open-run", nargs=2, metavar=("CASE", "DIRECTORY"),
                        help="Open the saved results for a built-in case at startup")
    data = Path(os.environ.get("LOCALAPPDATA", Path.home()))
    parser.add_argument("--runs", type=Path, default=data / "VelaTCAD" / "runs")
    parser.add_argument("--config", type=Path, default=data / "VelaTCAD" / "sds-settings.json",
                        help="Per-user JSON preferences file")
    args = parser.parse_args()
    if args.stage_cases:
        if "__compiled__" in globals():
            parser.error("Prepare case catalogs from the source checkout")
        prepare_cases(args.stage_cases.resolve())
        return
    launcher = args.launcher.resolve()
    if not launcher.is_file():
        parser.error(f"sd-sim launcher does not exist: {launcher}")
    cases = read_cases(args.cases.resolve())
    if not cases:
        parser.error("No workbench cases are available")
    saved_run = None
    if args.open_run:
        key, directory = args.open_run
        if key not in cases:
            parser.error(f"Unknown case ID: {key}")
        saved_run = key, Path(directory).resolve()
    runs = args.runs.resolve()
    runs.mkdir(parents=True, exist_ok=True)
    launch(launcher, cases, runs, saved_run, config_path=args.config.resolve())


if __name__ == "__main__":
    main()
