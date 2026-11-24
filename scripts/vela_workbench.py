"""Local Tk workbench for frozen sd-sim regression inputs."""
from __future__ import annotations

import argparse
import json
from math import isfinite
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
from vela_workbench_text import CASE_NAMES, GUIDE_EN, GUIDE_ZH, ZH


def solver_version(launcher):
    from zipfile import BadZipFile, ZipFile

    for directory in (launcher.parent / "bin", launcher.parent):
        metadata = directory / "frontend/vela_sdevice/_assets/version.json"
        archive = directory / "vela_sdevice.pyz"
        try:
            if metadata.is_file():
                version = json.loads(metadata.read_text(encoding="utf-8"))
            elif archive.is_file():
                with ZipFile(archive) as source:
                    version = json.loads(source.read("vela_sdevice/_assets/version.json"))
            else:
                continue
        except (OSError, KeyError, ValueError, BadZipFile):
            continue
        if isinstance(version, dict):
            return version
    return {}


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
    for label in ("File", "Simulation", "Settings", "Help"):
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
        menu = menus[group]
        menu.add_command(command=command)
        index = menu.index("end")
        if icon is None:
            actions[label] = (menu, index, None)
            return
        primary = icon == "run"
        size = max(18, settings.font_size + 6)
        normal = action_icon(root, icon, size, surface if primary else selection_ink)
        disabled = action_icon(root, icon, size, "#748397")
        button = ttk.Button(toolbar, image=(normal, "disabled", disabled), command=command,
                            style="Accent.TButton" if primary else "TButton")
        button.pack(side="left", padx=(0, 5))
        add_tooltip(button, lambda: action_label(label))
        actions[label] = (menu, index, button, normal, disabled)

    pane = ttk.PanedWindow(root, orient="horizontal")
    pane.pack(fill="both", expand=True, padx=16, pady=(0, 8))
    sidebar = ttk.Frame(pane, padding=(10, 9), style="Surface.TFrame")
    pane.add(sidebar)
    case_heading = ttk.Label(sidebar, style="Heading.Surface.TLabel")
    case_heading.pack(anchor="w", pady=(1, 8))
    device_options = (
        ("All", ""), ("PN2D", "sds-device-pn2d"), ("MOS / Gate", "sds-device-mos"),
        ("Schottky", "sds-device-schottky"), ("BJT", "sds-device-bjt"),
    )
    filter_row = ttk.Frame(sidebar, style="Surface.TFrame")
    filter_row.pack(fill="x", pady=(0, 5))
    device_heading = ttk.Label(filter_row, style="Surface.TLabel")
    device_heading.pack(side="left", padx=(0, 6))
    device_picker = ttk.Combobox(filter_row, state="readonly", width=13)
    device_picker.pack(side="right", fill="x", expand=True)
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
        figure = Figure(figsize=(5.7, 4.5), layout="compressed" if tab == cloud_tab else "constrained")
        canvas = FigureCanvasTkAgg(figure, master=tab)
        navigation = NavigationToolbar2Tk(canvas, tab, pack_toolbar=False)
        plot_zoom_buttons[tab] = add_plot_buttons(navigation, (
            ("+", "Zoom In", lambda: zoom_plot(tab, 0.8)),
            ("−", "Zoom Out", lambda: zoom_plot(tab, 1.25)),
        ))
        if tab == curve_tab:
            add_curve_tools(navigation)
        else:
            add_geometry_tools(navigation, tab)
        plots[tab] = figure, canvas, navigation
        return plots[tab]

    curve_notice = ttk.Label(curve_tab, style="Surface.TLabel", foreground=muted)
    curve_notice.bind("<Configure>", lambda event: curve_notice.configure(wraplength=event.width))
    curve_options = (
        ("Voltage", ("Inner voltage", "Outer voltage"), 12),
        ("X", ("Linear", "Log10"), 6),
        ("Y", ("Linear", "Log10"), 6),
    )
    curve_pickers = {}
    curve_style_options = (
        ("Line style", (("Solid", "-"), ("Dashed", "--"), ("Dash-dot", "-."), ("Dotted", ":")), 10),
        ("Marker", (("Automatic", "auto"), ("No marker", None), ("Circle", "o"),
                    ("Rectangle", "rectangle"), ("Triangle", "^"), ("Diamond", "D")), 11),
    )
    curve_style_pickers = {}
    plot_zoom_buttons = {}

    show_field_mesh = tk.BooleanVar(root, value=False)
    field_heading = field_picker = mesh_toggle = None

    def add_plot_buttons(navigation, actions):
        readout = next(widget for widget in navigation.pack_slaves() if widget.pack_info()["side"] == "right")
        buttons = []
        for icon, label, command in actions:
            button = tk.Button(navigation, text=icon, width=2, relief="flat",
                               background=navigation.cget("background"), foreground=ink,
                               activebackground=selection_fill, command=command)
            button.pack(side="left", before=readout, padx=(3, 0))
            add_tooltip(button, lambda label=label: tr(label))
            buttons.append(button)
        return buttons

    def plot_toolbar_row(navigation):
        controls = tk.Frame(navigation, background=navigation.cget("background"))
        # Reserve a row inside the toolbar before packing the native buttons and readout.
        controls.pack(side="bottom", before=navigation.pack_slaves()[0], fill="x", padx=6, pady=(0, 4))
        return controls

    def size_toolbar_row(navigation, controls):
        row_height = max(widget.winfo_reqheight() for widget in controls.pack_slaves())
        navigation.configure(height=int(navigation.cget("height")) + row_height + 4)

    def plot_picker(controls, label, width, command):
        heading = tk.Label(controls, text=tr(label), background=controls.cget("background"), foreground=ink)
        heading.pack(side="left", padx=(0, 4))
        picker = ttk.Combobox(controls, state="readonly", width=width)
        picker.pack(side="left", padx=(0, 8))
        picker.bind("<<ComboboxSelected>>", command)
        return heading, picker

    def add_curve_tools(navigation):
        controls = plot_toolbar_row(navigation)
        for label, _, width in curve_options:
            curve_pickers[label] = plot_picker(controls, label, width, select_curve_axes)
        size_toolbar_row(navigation, controls)
        appearance = plot_toolbar_row(navigation)
        for label, _, width in curve_style_options:
            curve_style_pickers[label] = plot_picker(appearance, label, width, select_curve_style)
        size_toolbar_row(navigation, appearance)

    def add_geometry_tools(navigation, tab):
        nonlocal field_heading, field_picker, mesh_toggle
        add_plot_buttons(navigation, (
            ("↶", "Rotate Left 90°", lambda: rotate_mesh(1)),
            ("↷", "Rotate Right 90°", lambda: rotate_mesh(-1)),
            ("0°", "Reset Rotation", rotate_mesh),
        ))
        if tab == mesh_tab:
            return
        controls = plot_toolbar_row(navigation)
        field_heading, field_picker = plot_picker(controls, "Field", 26, draw_field)
        style.configure("PlotToolbar.TCheckbutton", background=controls.cget("background"))
        mesh_toggle = ttk.Checkbutton(controls, text=tr("Show Mesh"), variable=show_field_mesh,
                                      style="PlotToolbar.TCheckbutton", command=toggle_field_mesh)
        mesh_toggle.pack(side="left")
        size_toolbar_row(navigation, controls)

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

    def text_view(tab, *, wrap="none", text_font=mono):
        frame = ttk.Frame(tab)
        frame.pack(fill="both", expand=True)
        frame.columnconfigure(0, weight=1)
        frame.rowconfigure(0, weight=1)
        view = tk.Text(frame, width=1, height=1, wrap=wrap, state="disabled", font=text_font, background="#F6F8FB",
                       foreground=ink, selectbackground=selection_fill, relief="flat", padx=9, pady=8)
        vertical = ttk.Scrollbar(frame, command=view.yview)
        view.configure(yscrollcommand=vertical.set)
        view.grid(row=0, column=0, sticky="nsew")
        vertical.grid(row=0, column=1, sticky="ns")
        if wrap == "none":
            horizontal = ttk.Scrollbar(frame, orient="horizontal", command=view.xview)
            view.configure(xscrollcommand=horizontal.set)
            horizontal.grid(row=1, column=0, sticky="ew")
        return view

    log_view = text_view(log_tab)
    script_view = text_view(script_tab)
    for tag, color in (
        ("keyword", accent), ("comment", "#566F56"), ("string", "#8A391A"), ("number", "#6D3F91"),
        ("directive", "#006A72"), ("placeholder", "#A2274D"),
    ):
        script_view.tag_configure(tag, foreground=color)
    script_view.tag_raise("sel")

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
    field_mesh_lines = []
    field_revision = None
    field_error = None
    view_key = None
    curve_axes = {}
    curve_styles = {}
    mesh_angles = {}

    def focused():
        return tree.focus() or (tree.selection() or ("",))[0]

    def filter_cases(_event=None):
        label = device_options[device_picker.current()][1]
        visible = [key for key, case in cases.items() if not label or label in case.labels]
        previous = focused()
        selected = set(tree.selection())
        tree.set_children("", *visible)
        for index, key in enumerate(visible):
            tree.item(key, tags=("alternate",) if index % 2 else ())
        current = previous if previous in visible else (visible[0] if visible else "")
        tree.focus(current)
        tree.selection_set([key for key in visible if key in selected] or ([current] if current else []))
        if current:
            tree.see(current)
        refresh()

    def sync_actions():
        states = (
            ("Run Current", not active and bool(focused())),
            ("Run Selected", not active and bool(tree.selection())),
            ("Stop Queue", active and not stop_pending.is_set()),
            ("Open Saved Run", not active and bool(focused())),
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

    def show_script(case):
        content = tr("Original script unavailable. Prepare a new case catalog.")
        if not case.original_script:
            replace_text(script_view, content)
            return
        try:
            content = (case.source / case.original_script).read_text(encoding="utf-8")
        except (OSError, UnicodeError) as error:
            replace_text(script_view, tr("Could not read original script: {}", error))
            return
        from vela_workbench_highlight import script_spans

        replace_text(script_view, content)
        for tag, row, start, end in script_spans(content):
            script_view.tag_add(tag, f"{row}.{start}", f"{row}.{end}")

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

    def finish_plot(tab):
        _, canvas, navigation = plots[tab]
        navigation.push_current()
        draw_plot(canvas)

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
        nonlocal field_mesh_lines
        field_mesh_lines = []
        selection = field_picker.current()
        if field is None or selection < 0:
            plot_message(cloud_tab)
            return
        from matplotlib.tri import Triangulation

        cloud_figure, _, axes = show_plot(cloud_tab)
        name = field_keys[selection]
        x, y = display_geometry(axes)
        mesh = Triangulation(x, y, triangles=field.triangles)
        image = axes.tripcolor(mesh, field.quantities[name], shading="gouraud")
        field_mesh_lines = axes.triplot(
            mesh, color="#E2E8F0", linewidth=.40, alpha=.55, visible=show_field_mesh.get(),
        )
        unit = field.metadata.get(f"Unit:{name}", "")
        cloud_figure.colorbar(image, ax=axes, label=unit, fraction=.055, pad=.03, panchor=False)
        stopped = field.metadata.get("GoalReached", field.metadata.get("goal_reached")) == 0
        bias = field.metadata.get(f"{field_terminal}_OuterVoltage_V")
        title = tr("Last accepted point" if stopped else "Final state")
        if bias is not None:
            title += "\n" + tr("{} outer voltage {:.6g} V", tr(field_terminal.capitalize()), bias)
        axes.set_title(title, fontsize="medium")
        finish_plot(cloud_tab)

    def toggle_field_mesh():
        for line in field_mesh_lines:
            line.set_visible(show_field_mesh.get())
        if field_mesh_lines:
            plots[cloud_tab][1].draw_idle()

    def refresh(_event=None):
        nonlocal view_key
        key = focused()
        directory = directories.get(key)
        sync_actions()
        selected = tabs.select()
        tab = tabs.nametowidget(selected) if selected else None
        if not key:
            view_key = None
            drawn.clear()
            result_heading.configure(text=tr("No matching cases"))
            log_heading.configure(text="")
            log_update()
            replace_text(script_view, "")
            if tab in placeholders:
                if tab == curve_tab:
                    curve_notice.pack_forget()
                plot_message(tab)
            return
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
            show_script(case)
        view_key = identity
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
        curve_notice.pack_forget()
        if not any(path.is_file() for path in curve_files):
            plot_message(curve_tab)
            return
        import numpy as np
        from vela_workbench_results import read_curve, read_stage_curves

        _, canvas, curve = show_plot(curve_tab)
        selections = curve_axes.get(focused(), (0, 0, 0))
        for (label, options, _), index in zip(curve_options, selections):
            picker = curve_pickers[label][1]
            picker.configure(values=tuple(tr(option) for option in options))
            picker.current(index)
        for (label, options, _), index in zip(curve_style_options, curve_styles.get(focused(), (0, 0))):
            picker = curve_style_pickers[label][1]
            picker.configure(values=tuple(tr(name) for name, _ in options))
            picker.current(index)
        for button in plot_zoom_buttons[curve_tab]:
            button.configure(state="disabled")

        def empty_curve(message):
            curve.set_axis_off()
            curve.text(.5, .5, message, transform=curve.transAxes, ha="center", va="center", wrap=True)
            draw_plot(canvas)

        try:
            voltage = ("inner", "outer")[selections[0]]
            if curve_path:
                points = read_curve(curve_path, voltage) if curve_path.is_file() else []
                curves = [(None, "Anode", points)] if points else []
            else:
                curves = read_stage_curves(work, stages, voltage) if work else []
        except (OSError, ValueError, KeyError) as error:
            empty_curve(tr("Could not read curve: {}", error))
            return

        log_axes = tuple(bool(selection) for selection in selections[1:])
        missing_inner = voltage == "inner" and any(x is None for _, _, points in curves for x, _ in points)
        omitted_by_log = False
        terminals = set()
        for index, terminal, points in curves:
            xy = np.asarray(points, dtype=float)
            finite = np.isfinite(xy).all(axis=1)
            visible = finite & (xy[:, log_axes] > 0).all(axis=1)
            omitted_by_log |= bool(np.any(finite & ~visible))
            if visible.any():
                xy[~visible] = np.nan
                label = tr("Stage {} · {}", index, terminal) if index is not None else terminal
                curve.plot(*xy.T, label=label)
                terminals.add(terminal.casefold())
        if not terminals:
            if omitted_by_log:
                empty_curve(tr("No points can be displayed with the selected log scales."))
            elif missing_inner:
                empty_curve(tr("Inner voltage unavailable. Select outer voltage."))
            else:
                empty_curve(tr("No curve available"))
            return

        notices = []
        if missing_inner:
            notices.append(tr("Inner voltage unavailable for some points."))
        if omitted_by_log:
            notices.append(tr("Log scales omit non-positive values."))
        if notices:
            curve_notice.configure(text=" ".join(notices))
            curve_notice.pack(fill="x", before=canvas.get_tk_widget(), pady=(0, 4))

        terminal = {"anode": "Anode", "drain": "Drain"}.get(next(iter(terminals)), "Swept terminal")
        if len(terminals) > 1:
            terminal = "Swept terminal"
        curve.set(
            xlabel=tr("{} inner voltage (V)" if voltage == "inner" else "{} outer voltage (V)", tr(terminal)),
            ylabel=tr("Anode total current (A/µm)" if curve_path else "Swept terminal current (A/µm)"),
            xscale="log" if log_axes[0] else "linear", yscale="log" if log_axes[1] else "linear",
        )
        apply_curve_style(curve)
        if not curve_path:
            curve.legend()
        curve.grid(alpha=.25)
        for button in plot_zoom_buttons[curve_tab]:
            button.configure(state="normal")
        finish_plot(curve_tab)

    def apply_curve_style(axes):
        line_index, marker_index = curve_styles.get(focused(), (0, 0))
        linestyle = curve_style_options[0][1][line_index][1]
        marker = curve_style_options[1][1][marker_index][1]
        if marker == "rectangle":
            from matplotlib.path import Path as MarkerPath

            marker = MarkerPath([(-.5, -.3), (.5, -.3), (.5, .3), (-.5, .3), (-.5, -.3)], closed=True)
        for line in axes.lines:
            line.set_linestyle(linestyle)
            if marker == "auto":
                line.set_marker("." if sum(isfinite(x) for x in line.get_xdata()) == 1 else "None")
            else:
                line.set_marker(marker if marker is not None else "None")

    def select_curve_style(_event):
        curve_styles[focused()] = tuple(picker.current() for _, picker in curve_style_pickers.values())
        figure, canvas, _ = plots[curve_tab]
        axes = figure.axes[0]
        apply_curve_style(axes)
        if axes.get_legend() is not None:
            axes.legend()
        draw_plot(canvas)

    def zoom_plot(tab, factor):
        import numpy as np

        figure, _, _ = plots[tab]
        axes = figure.axes[0]
        scale = axes.transScale
        limits = scale.transform(list(zip(axes.get_xlim(), axes.get_ylim())))
        center = limits[0] / 2 + limits[1] / 2
        with np.errstate(over="ignore", under="ignore", invalid="ignore"):
            updated = scale.inverted().transform(center + factor * (limits - center)).T
        if not np.isfinite(updated).all() or np.any(updated[:, 0] == updated[:, 1]):
            return
        for log, lower in zip((axes.get_xscale(), axes.get_yscale()), updated.min(axis=1)):
            if log == "log" and lower <= 0:
                return
        axes.set(xlim=updated[0], ylim=updated[1])
        finish_plot(tab)

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
            if field:
                _, _, mesh = show_plot(mesh_tab)
                x, y = display_geometry(mesh)
                mesh.triplot(x, y, field.triangles, linewidth=.45)
                finish_plot(mesh_tab)
            else:
                message = tr("Could not read mesh or field: {}", field_error) if field_error else None
                plot_message(mesh_tab, message)
            return
        if not field or not field.quantities:
            field_keys = []
            if field_picker is not None:
                field_picker.set("")
            message = tr("Could not read mesh or field: {}", field_error) if field_error else None
            plot_message(cloud_tab, message)
            return

        from vela_workbench_results import FIELDS

        plot_widgets(cloud_tab)
        index = field_picker.current()
        selected_field = field_keys[index] if index >= 0 else None
        field_keys = list(field.quantities)
        field_picker.configure(values=[tr(FIELDS[name]) for name in field_keys])
        index = field_keys.index(selected_field) if selected_field in field_keys else 0
        field_picker.current(index)
        draw_field()

    def select_curve_axes(_event):
        curve_axes[focused()] = tuple(picker.current() for _, picker in curve_pickers.values())
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
            if button is not None:
                button.configure(text=text)
        case_heading.configure(text=tr("Built-in cases"))
        selection = max(0, device_picker.current())
        device_heading.configure(text=tr("Device type"))
        device_picker.configure(values=[tr(name) for name, _ in device_options])
        device_picker.current(selection)
        tree.heading("#0", text=tr("Case"))
        tree.heading("status", text=tr("Status"))
        for key, status in statuses.items():
            update_case(key, status)
        for tab, label in (
            (curve_tab, "Curve"), (mesh_tab, "Mesh"), (cloud_tab, "Final Field"),
            (log_tab, "Log"), (script_tab, "Sentaurus Script"),
        ):
            tabs.tab(tab, text=tr(label))
        for pickers in (curve_pickers, curve_style_pickers):
            for label, (heading, _) in pickers.items():
                heading.configure(text=tr(label))
        if field_picker is not None:
            mesh_toggle.configure(text=tr("Show Mesh"))
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

    def open_guide():
        dialog = tk.Toplevel(root)
        dialog.title(tr("Documentation"))
        dialog.transient(root)
        dialog.geometry("620x450")
        dialog.minsize(400, 260)
        panel = ttk.Frame(dialog, padding=14)
        panel.pack(fill="both", expand=True)
        ttk.Button(panel, text=tr("Close"), command=dialog.destroy).pack(side="bottom", anchor="e", pady=(8, 0))
        view = text_view(panel, wrap="word", text_font="TkTextFont")
        replace_text(view, GUIDE_ZH if chinese else GUIDE_EN)
        dialog.bind("<Escape>", lambda _event: dialog.destroy())
        dialog.wait_visibility()
        dialog.grab_set()
        view.focus_set()

    def show_about():
        version = solver_version(launcher)
        text = "SDS · " + tr("Semiconductor Device Simulator")
        text += "\n\n" + tr("Desktop workbench for semiconductor device simulations with sd-sim.")
        if version.get("version") and version.get("commit"):
            text += "\n\n" + tr("Solver version {} ({})", version["version"], version["commit"][:12])
        messagebox.showinfo(tr("About SDS"), text, parent=root)

    add_action("Simulation", "Run Current", "run", lambda: start((focused(),)))
    add_action("Simulation", "Run Selected", "batch", lambda: start(tree.selection()))
    add_action("Simulation", "Stop Queue", "stop", stop_queue)
    ttk.Separator(toolbar, orient="vertical").pack(side="left", fill="y", padx=(1, 7), pady=3)
    add_action("File", "Open Saved Run", "saved", open_saved)
    add_action("File", "Open Run Folder", "folder", open_run_folder)
    menus["File"].add_separator()
    add_action("File", "Exit", None, close)
    ttk.Separator(toolbar, orient="vertical").pack(side="left", fill="y", padx=(1, 7), pady=3)
    add_action("Settings", "Open Settings File", "settings", open_settings)
    add_action("Settings", "Language", "language", switch_language)
    add_action("Help", "Documentation", None, open_guide)
    add_action("Help", "About SDS", None, show_about)
    root.protocol("WM_DELETE_WINDOW", close)
    device_picker.bind("<<ComboboxSelected>>", filter_cases)
    tree.bind("<<TreeviewSelect>>", refresh)
    tabs.bind("<<NotebookTabChanged>>", refresh)
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
