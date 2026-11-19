"""Icons and hover hints for the SDS command bar."""


def action_icon(master, name, size, color):
    from PIL import Image, ImageDraw, ImageTk

    image = Image.new("RGBA", (48, 48))
    draw = ImageDraw.Draw(image)

    def stroke(points):
        draw.line(points, fill=color, width=4, joint="curve")

    if name == "saved":
        stroke([(12, 42), (12, 6), (29, 6), (37, 14), (37, 42), (12, 42)])
        stroke([(28, 6), (28, 15), (36, 15)])
        stroke([(19, 24), (30, 24)])
        stroke([(19, 32), (30, 32)])
    elif name == "folder":
        stroke([(7, 36), (7, 11), (21, 11), (25, 16), (40, 16), (40, 21)])
        stroke([(6, 38), (14, 22), (44, 22), (37, 38), (6, 38)])
    elif name == "run":
        draw.polygon([(15, 7), (39, 24), (15, 41)], fill=color)
    elif name == "batch":
        draw.polygon([(7, 10), (26, 24), (7, 38)], fill=color)
        draw.polygon([(24, 10), (43, 24), (24, 38)], fill=color)
    elif name == "stop":
        draw.rounded_rectangle((11, 11, 37, 37), radius=3, fill=color)
    elif name == "settings":
        draw.ellipse((10, 10, 38, 38), outline=color, width=5)
        draw.ellipse((20, 20, 28, 28), outline=color, width=3)
        for points in (
            [(24, 4), (24, 11)], [(24, 37), (24, 44)], [(4, 24), (11, 24)], [(37, 24), (44, 24)],
            [(10, 10), (15, 15)], [(33, 33), (38, 38)], [(10, 38), (15, 33)], [(33, 15), (38, 10)],
        ):
            stroke(points)
    elif name == "language":
        stroke([(4, 33), (12, 11), (20, 33)])
        stroke([(7, 26), (17, 26)])
        stroke([(23, 15), (43, 15)])
        stroke([(32, 8), (34, 14)])
        stroke([(39, 16), (35, 26), (23, 37)])
        stroke([(27, 20), (32, 29), (43, 37)])
    return ImageTk.PhotoImage(image.resize((size, size), Image.Resampling.LANCZOS), master=master)


def add_tooltip(widget, get_text):
    import tkinter as tk
    from tkinter import ttk

    timer = None
    popup = None

    def show():
        nonlocal timer, popup
        timer = None
        popup = tk.Toplevel(widget)
        popup.withdraw()
        popup.overrideredirect(True)
        popup.transient(widget.winfo_toplevel())
        ttk.Label(popup, text=get_text(), style="Tooltip.TLabel").pack()
        popup.update_idletasks()
        x = widget.winfo_rootx()
        y = widget.winfo_rooty() + widget.winfo_height() + 5
        popup.geometry(f"+{max(0, x)}+{max(0, y)}")
        popup.deiconify()

    def hide(_event=None):
        nonlocal timer, popup
        if timer is not None:
            widget.after_cancel(timer)
            timer = None
        if popup is not None:
            popup.destroy()
            popup = None

    def schedule(event):
        nonlocal timer
        hide(event)
        timer = widget.after(500, show)

    for sequence in ("<Enter>", "<FocusIn>"):
        widget.bind(sequence, schedule, add="+")
    for sequence in ("<Leave>", "<FocusOut>", "<ButtonPress>", "<Destroy>"):
        widget.bind(sequence, hide, add="+")
