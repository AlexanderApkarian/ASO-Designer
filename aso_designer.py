from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

from aso_logic import (
    CHEMISTRY_PRESETS,
    DEFAULT_CHEMISTRY,
    LINKAGE_OPTIONS,
    NUCLEOBASE_MODIFICATION_OPTIONS,
    PENALTY_BASE_MODES,
    PENALTY_POSITION_MODES,
    RIBOSE_MODIFICATION_OPTIONS,
    AsoInputError,
    AsoInputs,
    PenaltyAsoInputs,
    build_gapmer_pattern,
    chemistry_optimization_walk,
    chemistry_display_label,
    combine_base_modification,
    convert_sequences_to_idt,
    generate_design,
    generate_penalty_design,
    resolve_chemistry,
    mutation_header_indexes,
    normalise_mutation_type,
    split_base_modification,
    summarise_linkages,
)


EXAMPLE_SEQUENCE = (
    "AUGCUACGUAUGCUACGUAUGGCAUCGUAUGCUACGUAUGCUACGUA"
)
APP_NAME = "ASO Designer - by Alexander Apkarian"
BUBBLE_EXPORT_SCALE = 4
BUBBLE_EXPORT_DPI = 600


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Design ASO microwalks and export IDT order codes.")
    parser.add_argument("--export", help="Write an .xlsx result workbook instead of opening the GUI.")
    parser.add_argument("--example", action="store_true", help="Use a small built-in example input.")
    parser.add_argument("--target-gene", default="")
    parser.add_argument("--snp-identifier", default="")
    parser.add_argument("--chemistry-number", default="")
    parser.add_argument("--aso-chemistry", default=DEFAULT_CHEMISTRY, choices=list(CHEMISTRY_PRESETS))
    parser.add_argument("--gap-length", type=int, default=12)
    parser.add_argument("--wing-length", type=int, default=3)
    parser.add_argument("--wing-chemistry", default="LNA", choices=list(RIBOSE_MODIFICATION_OPTIONS) + ["None"])
    parser.add_argument("--backbone-modification", default="PS", choices=["PS", "PO", "MIXED", "None"])
    parser.add_argument("--mutation-type", default="Insertion", choices=["Insertion", "Deletion", "Substitution"])
    parser.add_argument("--mutation-length", type=int, default=1)
    parser.add_argument("--mutation-start", type=int, default=21)
    parser.add_argument("--microwalk-step-size", type=int, default=1)
    parser.add_argument("--rna-sequence", default="")
    parser.add_argument("--rna-file", help="Read the RNA sequence from a text file.")
    return parser.parse_args(argv)


def inputs_from_args(args: argparse.Namespace) -> AsoInputs:
    sequence = args.rna_sequence
    if args.rna_file:
        sequence = Path(args.rna_file).read_text(encoding="utf-8")
    if args.example and not sequence:
        sequence = EXAMPLE_SEQUENCE
    target_gene = args.target_gene or ("GENE" if args.example else "")
    snp_identifier = args.snp_identifier or ("Example" if args.example else "")
    chemistry_number = args.chemistry_number or ("C1" if args.example else "")
    return AsoInputs(
        target_gene=target_gene,
        snp_identifier=snp_identifier,
        chemistry_number=chemistry_number,
        aso_chemistry=args.aso_chemistry,
        gap_length=args.gap_length,
        wing_length=args.wing_length,
        wing_chemistry=args.wing_chemistry,
        backbone_modification=args.backbone_modification,
        mutation_type=args.mutation_type,
        mutation_length=args.mutation_length,
        mutation_start=args.mutation_start,
        microwalk_step_size=args.microwalk_step_size,
        rna_sequence=sequence,
    )


class AsoDesignerApp:
    def __init__(self) -> None:
        import tkinter as tk
        from tkinter import filedialog, messagebox, ttk

        self.tk = tk
        self.ttk = ttk
        self.filedialog = filedialog
        self.messagebox = messagebox

        self.root = tk.Tk()
        self.root.title(APP_NAME)
        self.root.geometry("1220x820")
        self.root.minsize(980, 640)

        self.vars: dict[str, tk.StringVar] = {}
        self.last_result = None
        self.show_bubble_base_letters = tk.BooleanVar(value=True)
        self.show_chemopt_bubble_base_letters = tk.BooleanVar(value=True)
        self.show_penalty_bubble_base_letters = tk.BooleanVar(value=True)
        self.custom_base_modifications: tuple[str, ...] = ()
        self.custom_linkages: tuple[str, ...] = ()
        self.converter_custom_base_modifications: tuple[str, ...] = ()
        self.converter_custom_linkages: tuple[str, ...] = ()
        self.penalty_custom_base_modifications: tuple[str, ...] = ()
        self.penalty_custom_linkages: tuple[str, ...] = ()
        self.last_idt_conversion_rows = ()
        self.last_idt_conversion_chemistry = None
        self.last_penalty_result = None
        self.chemopt_base_modifications: tuple[str, ...] = ()
        self.chemopt_linkages: tuple[str, ...] = ()
        self.chemopt_motif_positions: tuple[int, ...] = ()
        self.last_chemopt_rows = ()
        self.last_chemopt_chemistry_label = ""
        self.custom_editor_window = None
        self.converter_custom_editor_window = None
        self.penalty_custom_editor_window = None
        self.chemopt_editor_window = None
        self._home_click_times: dict[str, float] = {}
        self._variant_bubble_dirty = False
        self._variant_render_token = 0

        self._build_ui()
        self._apply_chemistry_preset()
        self._apply_converter_chemistry_preset()
        self._apply_penalty_chemistry_preset()
        self._apply_chemopt_chemistry_preset()

    def run(self) -> None:
        self.root.mainloop()

    def _var(self, name: str, value: str = ""):
        var = self.tk.StringVar(value=value)
        self.vars[name] = var
        return var

    def _build_ui(self) -> None:
        tk = self.tk
        ttk = self.ttk

        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)

        self.home_frame = ttk.Frame(self.root)
        self.variant_frame = ttk.Frame(self.root)
        self.chemistry_optimization_frame = ttk.Frame(self.root)
        self.idt_converter_frame = ttk.Frame(self.root)
        self.penalty_design_frame = ttk.Frame(self.root)
        self.home_frame.grid(row=0, column=0, sticky="nsew")
        self.variant_frame.grid(row=0, column=0, sticky="nsew")
        self.chemistry_optimization_frame.grid(row=0, column=0, sticky="nsew")
        self.idt_converter_frame.grid(row=0, column=0, sticky="nsew")
        self.penalty_design_frame.grid(row=0, column=0, sticky="nsew")

        self._build_home_page(self.home_frame)
        self._build_variant_tool(self.variant_frame)
        self._build_chemistry_optimization_tool(self.chemistry_optimization_frame)
        self._build_idt_converter_tool(self.idt_converter_frame)
        self._build_penalty_design_tool(self.penalty_design_frame)
        self._show_home()

    def _build_home_page(self, parent) -> None:
        ttk = self.ttk

        parent.columnconfigure(0, weight=1)
        parent.rowconfigure(0, weight=1)
        parent.rowconfigure(2, weight=1)

        content = ttk.Frame(parent, padding=36)
        content.grid(row=1, column=0, sticky="ew")
        content.columnconfigure(0, weight=1)

        ttk.Label(
            content,
            text="ASO Designer",
            font=("TkDefaultFont", 28, "bold"),
            foreground="#1f4e79",
            anchor="center",
        ).grid(row=0, column=0, sticky="ew", pady=(0, 8))
        ttk.Label(
            content,
            text="By Alexander Apkarian",
            font=("TkDefaultFont", 13),
            foreground="#6b7280",
            anchor="center",
        ).grid(row=1, column=0, sticky="ew", pady=(0, 28))

        tiles = ttk.Frame(content)
        tiles.grid(row=2, column=0, sticky="ew")
        for column in range(4):
            tiles.columnconfigure(column, weight=1, uniform="home_tools")

        self._home_tool_tile(
            tiles,
            0,
            "Variant Microwalk Tool",
            self._show_variant_tool,
            "variant",
        )
        self._home_tool_tile(
            tiles,
            1,
            "Chemistry Optimisation Tool",
            self._show_chemistry_optimization_tool,
            "chemistry",
        )
        self._home_tool_tile(
            tiles,
            2,
            "IDT notation converter",
            self._show_idt_converter_tool,
            "converter",
        )
        self._home_tool_tile(
            tiles,
            3,
            "Penalty ASO Design Tool",
            self._show_penalty_design_tool,
            "penalty",
        )

    def _home_tool_tile(self, parent, column: int, title: str, command, icon_kind: str) -> None:
        tk = self.tk
        tile = tk.Frame(
            parent,
            background="#ffffff",
            highlightbackground="#d1d5db",
            highlightcolor="#1f4e79",
            highlightthickness=1,
            padx=14,
            pady=14,
            cursor="hand2",
        )
        tile.grid(row=0, column=column, sticky="nsew", padx=8, pady=4)
        tile.columnconfigure(0, weight=1)

        canvas = tk.Canvas(tile, width=240, height=128, background="#ffffff", highlightthickness=0, cursor="hand2")
        canvas.grid(row=0, column=0, pady=(0, 12))
        self._draw_home_tool_icon(canvas, icon_kind)

        label = tk.Label(
            tile,
            text=title,
            background="#ffffff",
            foreground="#111827",
            font=("TkDefaultFont", 14, "bold"),
            justify="center",
            wraplength=230,
            cursor="hand2",
        )
        label.grid(row=1, column=0, sticky="ew")

        def open_tool(_event=None):
            now = time.monotonic()
            last_click = self._home_click_times.get(title, 0.0)
            if now - last_click < 0.35:
                return "break"
            self._home_click_times[title] = now
            self.root.after_idle(command)
            return "break"

        def highlight(_event=None):
            tile.configure(highlightbackground="#1f4e79")

        def unhighlight(_event=None):
            tile.configure(highlightbackground="#d1d5db")

        for widget in (tile, canvas, label):
            widget.bind("<Button-1>", open_tool)
            widget.bind("<ButtonRelease-1>", open_tool)
            widget.bind("<Enter>", highlight)
            widget.bind("<Leave>", unhighlight)

    def _draw_home_tool_icon(self, canvas, icon_kind: str) -> None:
        canvas.delete("all")
        if icon_kind == "variant":
            self._draw_variant_home_icon(canvas)
        elif icon_kind == "chemistry":
            self._draw_chemistry_home_icon(canvas)
        elif icon_kind == "penalty":
            self._draw_penalty_home_icon(canvas)
        else:
            self._draw_converter_home_icon(canvas)

    def _draw_variant_home_icon(self, canvas) -> None:
        pitch = 20
        start_x = 66
        r = 8
        row_colors = ("#cfe2f3", "#d9ead3", "#f4cccc")
        for row, y in enumerate((40, 66, 92)):
            offset = row * 14
            canvas.create_text(40, y, text=f"A{row + 1}", font=("TkDefaultFont", 8, "bold"), fill="#4b5563")
            for idx in range(6):
                x = start_x + offset + idx * pitch
                canvas.create_oval(x - r, y - r, x + r, y + r, fill=row_colors[row], outline="#374151")

    def _draw_chemistry_home_icon(self, canvas) -> None:
        pitch = 20
        start_x = 58
        r = 8
        for row, y in enumerate((38, 66, 94)):
            canvas.create_text(34, y, text=f"C{row + 1}", font=("TkDefaultFont", 8, "bold"), fill="#4b5563")
            motif_start = 2 + row
            for idx in range(8):
                x = start_x + idx * pitch
                fill = "#cfe2f3" if idx < 2 or idx > 5 else "#f8fafc"
                if motif_start <= idx <= motif_start + 1:
                    fill = "#f4cccc"
                canvas.create_oval(x - r, y - r, x + r, y + r, fill=fill, outline="#374151")
                if idx < 7:
                    diamond_x = x + pitch / 2
                    canvas.create_polygon(
                        diamond_x,
                        y - 5,
                        diamond_x + 5,
                        y,
                        diamond_x,
                        y + 5,
                        diamond_x - 5,
                        y,
                        fill="#111827",
                        outline="#111827",
                    )

    def _draw_converter_home_icon(self, canvas) -> None:
        pitch = 20
        start_x = 24
        y = 46
        for idx, base in enumerate("ACGU"):
            x = start_x + idx * pitch
            canvas.create_oval(x - 8, y - 8, x + 8, y + 8, fill="#f8fafc", outline="#374151")
            canvas.create_text(x, y, text=base, font=("TkDefaultFont", 8, "bold"), fill="#111827")

        canvas.create_line(112, y, 150, y, fill="#1f4e79", width=3, arrow="last")
        canvas.create_rectangle(164, 25, 222, 67, fill="#d9ead3", outline="#6aa84f", width=2)
        canvas.create_text(193, 39, text="IDT", font=("TkDefaultFont", 13, "bold"), fill="#274e13")
        canvas.create_text(193, 56, text="/MOE/*", font=("TkDefaultFont", 8, "bold"), fill="#274e13")
        canvas.create_line(42, 88, 202, 88, fill="#9ca3af", width=2)
        for idx, token in enumerate(("/5", "A", "*", "C", "/3")):
            x = 62 + idx * 30
            canvas.create_text(x, 104, text=token, font=("TkDefaultFont", 8, "bold"), fill="#4b5563")

    def _draw_penalty_home_icon(self, canvas) -> None:
        pitch = 20
        start_x = 54
        r = 8
        for row, y in enumerate((40, 68, 96)):
            offset = row * 10
            canvas.create_text(34, y, text=f"P{row + 1}", font=("TkDefaultFont", 8, "bold"), fill="#4b5563")
            penalty_idx = 3 + row
            for idx in range(7):
                x = start_x + offset + idx * pitch
                fill = "#f8fafc"
                outline = "#374151"
                width = 1
                if idx == penalty_idx:
                    fill = "#f6c343"
                    outline = "#c00000"
                    width = 2
                canvas.create_oval(x - r, y - r, x + r, y + r, fill=fill, outline=outline, width=width)
                canvas.create_text(x, y, text="G" if idx == penalty_idx else "A", font=("TkDefaultFont", 7, "bold"), fill="#111827")

    def _show_home(self) -> None:
        self.home_frame.tkraise()

    def _show_variant_tool(self) -> None:
        self.variant_frame.tkraise()

    def _show_chemistry_optimization_tool(self) -> None:
        self.chemistry_optimization_frame.tkraise()

    def _show_idt_converter_tool(self) -> None:
        self.idt_converter_frame.tkraise()

    def _show_penalty_design_tool(self) -> None:
        self.penalty_design_frame.tkraise()

    def _show_tool_coming_soon(self, tool_name: str) -> None:
        self.messagebox.showinfo(APP_NAME, f"{tool_name} is not built yet.")

    def _build_variant_tool(self, parent) -> None:
        tk = self.tk
        ttk = self.ttk

        parent.columnconfigure(0, weight=0)
        parent.columnconfigure(1, weight=1)
        parent.rowconfigure(0, weight=1)
        parent.rowconfigure(1, weight=0)

        left = ttk.Frame(parent, padding=12)
        left.grid(row=0, column=0, sticky="ns")
        left.columnconfigure(1, weight=1)

        right = ttk.Frame(parent, padding=(0, 12, 12, 12))
        right.grid(row=0, column=1, sticky="nsew")
        right.rowconfigure(0, weight=1)
        right.columnconfigure(0, weight=1)

        credit = ttk.Label(parent, text="By Alexander Apkarian", font=("TkDefaultFont", 11), foreground="#6b7280")
        credit.grid(row=1, column=0, columnspan=2, sticky="e", padx=12, pady=(0, 6))

        row = 0
        nav_row = ttk.Frame(left)
        nav_row.grid(row=row, column=0, columnspan=2, sticky="ew", pady=(0, 8))
        nav_row.columnconfigure(1, weight=1)
        ttk.Button(nav_row, text="Home", command=self._show_home).grid(row=0, column=0, sticky="w", padx=(0, 8))
        ttk.Label(nav_row, text="Variant Microwalk Tool", font=("TkDefaultFont", 14, "bold")).grid(
            row=0,
            column=1,
            sticky="e",
        )
        row += 1

        title = ttk.Label(left, text="Inputs", font=("TkDefaultFont", 14, "bold"))
        title.grid(row=row, column=0, columnspan=2, sticky="w", pady=(0, 8))
        row += 1

        self._entry(left, row, "Target Gene Name", "target_gene")
        row += 1
        self._entry(left, row, "SNP Identifier", "snp_identifier")
        row += 1
        self._entry(left, row, "Chemistry Identifier", "chemistry_number")
        row += 1

        self._combo(left, row, "ASO Chemistry", "aso_chemistry", list(CHEMISTRY_PRESETS), DEFAULT_CHEMISTRY)
        self.vars["aso_chemistry"].trace_add("write", lambda *_: self._on_chemistry_changed())
        row += 1
        self.custom_button = ttk.Button(left, text="Edit Custom Pattern", command=self._open_custom_chemistry_editor)
        self.custom_button.grid(row=row, column=1, sticky="ew", pady=3)
        row += 1
        self.gap_entry = self._entry(left, row, "Central Gap Length", "gap_length", "12")
        row += 1
        self.wing_entry = self._entry(left, row, "Modified Wing Length", "wing_length", "3")
        row += 1
        self.wing_combo = self._combo(left, row, "Modified Wing Chemistry", "wing_chemistry", list(RIBOSE_MODIFICATION_OPTIONS) + ["None"], "LNA")
        row += 1
        self.backbone_combo = self._combo(left, row, "Backbone Modification", "backbone_modification", ["PS", "PO", "MIXED", "None"], "PS")
        row += 1

        ttk.Separator(left).grid(row=row, column=0, columnspan=2, sticky="ew", pady=10)
        row += 1
        self._combo(left, row, "Mutation Type", "mutation_type", ["Insertion", "Deletion", "Substitution"], "Insertion")
        row += 1
        self._entry(left, row, "Indel/Substitution Length", "mutation_length", "1")
        row += 1
        self._entry(left, row, "Indel/Substitution Start Position (first base = 0)", "mutation_start", "21")
        row += 1
        self._entry(left, row, "Microwalk Step Size", "microwalk_step_size", "1")
        row += 1

        ttk.Label(left, text="RNA 5' to 3' (hyphens optional)").grid(row=row, column=0, columnspan=2, sticky="w")
        row += 1
        self.sequence_text = tk.Text(left, width=44, height=12, wrap="word", font=("Menlo", 11))
        self.sequence_text.grid(row=row, column=0, columnspan=2, sticky="nsew")
        left.rowconfigure(row, weight=1)
        row += 1

        button_row = ttk.Frame(left)
        button_row.grid(row=row, column=0, columnspan=2, sticky="ew", pady=(10, 0))
        button_row.columnconfigure(0, weight=1)
        button_row.columnconfigure(1, weight=1)
        ttk.Button(button_row, text="Calculate", command=self.calculate).grid(row=0, column=0, sticky="ew", padx=(0, 6))
        ttk.Button(button_row, text="Export Excel", command=self.export_excel).grid(row=0, column=1, sticky="ew")

        notebook = ttk.Notebook(right)
        self.variant_notebook = notebook
        notebook.grid(row=0, column=0, sticky="nsew")

        results_frame = ttk.Frame(notebook)
        results_frame.rowconfigure(0, weight=1)
        results_frame.columnconfigure(0, weight=1)
        notebook.add(results_frame, text="ASO Table")

        self.table_text = tk.Text(results_frame, wrap="none", font=("Menlo", 11), undo=False)
        self.table_text.configure(tabs=("160", "520", "760", "1020"))
        self.table_text.tag_configure("table_header", background="#1f4e79", foreground="#ffffff")
        self.table_text.tag_configure("display_mutation", foreground="#c00000")
        table_y = ttk.Scrollbar(results_frame, orient="vertical", command=self.table_text.yview)
        table_x = ttk.Scrollbar(results_frame, orient="horizontal", command=self.table_text.xview)
        self.table_text.configure(yscrollcommand=table_y.set, xscrollcommand=table_x.set)
        self.table_text.grid(row=0, column=0, sticky="nsew")
        table_y.grid(row=0, column=1, sticky="ns")
        table_x.grid(row=1, column=0, sticky="ew")

        align_frame = ttk.Frame(notebook)
        align_frame.rowconfigure(0, weight=1)
        align_frame.columnconfigure(0, weight=1)
        notebook.add(align_frame, text="Alignment")
        self.alignment_canvas = tk.Canvas(align_frame, background="#ffffff", highlightthickness=0)
        self.alignment_canvas.grid(row=0, column=0, sticky="nsew")
        align_y = ttk.Scrollbar(align_frame, orient="vertical", command=self.alignment_canvas.yview)
        align_x = ttk.Scrollbar(align_frame, orient="horizontal", command=self.alignment_canvas.xview)
        self.alignment_canvas.configure(yscrollcommand=align_y.set, xscrollcommand=align_x.set)
        align_y.grid(row=0, column=1, sticky="ns")
        align_x.grid(row=1, column=0, sticky="ew")

        bubble_frame = ttk.Frame(notebook)
        self.variant_bubble_frame = bubble_frame
        bubble_frame.rowconfigure(1, weight=1)
        bubble_frame.columnconfigure(0, weight=1)
        notebook.add(bubble_frame, text="Visualisation")
        notebook.bind("<<NotebookTabChanged>>", lambda _event: self._render_variant_bubble_if_visible())
        bubble_toolbar = ttk.Frame(bubble_frame, padding=(0, 0, 0, 6))
        bubble_toolbar.grid(row=0, column=0, columnspan=2, sticky="ew")
        bubble_toolbar.columnconfigure(0, weight=1)
        ttk.Checkbutton(
            bubble_toolbar,
            text="Show base letters",
            variable=self.show_bubble_base_letters,
            command=self._refresh_bubble_figure,
        ).grid(row=0, column=1, sticky="e", padx=(0, 12))
        ttk.Button(bubble_toolbar, text="Save HD Image", command=self.save_bubble_figure_image).grid(row=0, column=2, sticky="e")
        self.bubble_canvas = tk.Canvas(bubble_frame, background="#ffffff", highlightthickness=0)
        self.bubble_canvas.grid(row=1, column=0, sticky="nsew")
        bubble_y = ttk.Scrollbar(bubble_frame, orient="vertical", command=self.bubble_canvas.yview)
        bubble_x = ttk.Scrollbar(bubble_frame, orient="horizontal", command=self.bubble_canvas.xview)
        self.bubble_canvas.configure(yscrollcommand=bubble_y.set, xscrollcommand=bubble_x.set)
        bubble_y.grid(row=1, column=1, sticky="ns")
        bubble_x.grid(row=2, column=0, sticky="ew")

    def _build_chemistry_optimization_tool(self, parent) -> None:
        tk = self.tk
        ttk = self.ttk

        parent.columnconfigure(0, weight=0)
        parent.columnconfigure(1, weight=1)
        parent.rowconfigure(0, weight=1)
        parent.rowconfigure(1, weight=0)

        left = ttk.Frame(parent, padding=12)
        left.grid(row=0, column=0, sticky="ns")
        left.columnconfigure(1, weight=1)

        right = ttk.Frame(parent, padding=(0, 12, 12, 12))
        right.grid(row=0, column=1, sticky="nsew")
        right.rowconfigure(1, weight=1)
        right.columnconfigure(0, weight=1)

        credit = ttk.Label(parent, text="By Alexander Apkarian", font=("TkDefaultFont", 11), foreground="#6b7280")
        credit.grid(row=1, column=0, columnspan=2, sticky="e", padx=12, pady=(0, 6))

        row = 0
        nav_row = ttk.Frame(left)
        nav_row.grid(row=row, column=0, columnspan=2, sticky="ew", pady=(0, 8))
        nav_row.columnconfigure(1, weight=1)
        ttk.Button(nav_row, text="Home", command=self._show_home).grid(row=0, column=0, sticky="w", padx=(0, 8))
        ttk.Label(nav_row, text="Chemistry Optimisation Tool", font=("TkDefaultFont", 14, "bold")).grid(
            row=0,
            column=1,
            sticky="e",
        )
        row += 1

        self._var("chemopt_aso_chemistry", "5-10-5 MOE/DNA")
        self._var("chemopt_gap_length", "10")
        self._var("chemopt_wing_length", "5")
        self._var("chemopt_wing_chemistry", "MOE")
        self._var("chemopt_backbone_modification", "PS")
        self._var("chemopt_step_size", "1")
        self.chemopt_custom_button = ttk.Button(
            left,
            text="Select Chemistry and Walk Motif",
            command=self._open_chemopt_chemistry_editor,
        )
        self.chemopt_custom_button.grid(row=row, column=0, columnspan=2, sticky="ew", pady=(4, 3))
        row += 1

        self.chemopt_motif_status = ttk.Label(
            left,
            text="Walk motif: not selected",
            foreground="#4b5563",
            wraplength=330,
        )
        self.chemopt_motif_status.grid(row=row, column=0, columnspan=2, sticky="w", pady=(6, 0))
        row += 1

        self._entry(left, row, "Microwalk Step Size", "chemopt_step_size", "1")
        row += 1

        ttk.Separator(left).grid(row=row, column=0, columnspan=2, sticky="ew", pady=10)
        row += 1
        ttk.Label(left, text="ASO sequence 5' to 3'").grid(row=row, column=0, columnspan=2, sticky="w")
        row += 1
        ttk.Label(
            left,
            text="Enter the fixed ASO sequence. The selected chemistry motif will be walked along the central core only.",
            foreground="#4b5563",
            wraplength=330,
        ).grid(row=row, column=0, columnspan=2, sticky="w", pady=(0, 4))
        row += 1
        self.chemopt_sequence_text = tk.Text(left, width=44, height=6, wrap="none", font=("Menlo", 11))
        self.chemopt_sequence_text.grid(row=row, column=0, columnspan=2, sticky="nsew")
        left.rowconfigure(row, weight=1)
        row += 1

        ttk.Button(left, text="Generate Chemistry Walk", command=self.generate_chemistry_optimization).grid(
            row=row,
            column=0,
            columnspan=2,
            sticky="ew",
            pady=(10, 0),
        )

        ttk.Label(right, text="Chemistry optimisation output", font=("TkDefaultFont", 14, "bold")).grid(
            row=0,
            column=0,
            sticky="w",
            pady=(0, 8),
        )
        chemopt_notebook = ttk.Notebook(right)
        chemopt_notebook.grid(row=1, column=0, sticky="nsew")

        output_frame = ttk.Frame(chemopt_notebook)
        output_frame.rowconfigure(0, weight=1)
        output_frame.columnconfigure(0, weight=1)
        chemopt_notebook.add(output_frame, text="Output Table")
        self.chemopt_output_text = tk.Text(output_frame, wrap="none", font=("Menlo", 11), undo=False)
        self.chemopt_output_text.configure(tabs=("90", "250", "410", "650", "960"))
        output_y = ttk.Scrollbar(output_frame, orient="vertical", command=self.chemopt_output_text.yview)
        output_x = ttk.Scrollbar(output_frame, orient="horizontal", command=self.chemopt_output_text.xview)
        self.chemopt_output_text.configure(yscrollcommand=output_y.set, xscrollcommand=output_x.set)
        self.chemopt_output_text.grid(row=0, column=0, sticky="nsew")
        output_y.grid(row=0, column=1, sticky="ns")
        output_x.grid(row=1, column=0, sticky="ew")

        chemopt_bubble_frame = ttk.Frame(chemopt_notebook)
        chemopt_bubble_frame.rowconfigure(1, weight=1)
        chemopt_bubble_frame.columnconfigure(0, weight=1)
        chemopt_notebook.add(chemopt_bubble_frame, text="Visualisation")
        chemopt_bubble_toolbar = ttk.Frame(chemopt_bubble_frame, padding=(0, 0, 0, 6))
        chemopt_bubble_toolbar.grid(row=0, column=0, columnspan=2, sticky="ew")
        chemopt_bubble_toolbar.columnconfigure(0, weight=1)
        ttk.Checkbutton(
            chemopt_bubble_toolbar,
            text="Show base letters",
            variable=self.show_chemopt_bubble_base_letters,
            command=self._refresh_chemopt_bubble_figure,
        ).grid(row=0, column=1, sticky="e", padx=(0, 12))
        ttk.Button(
            chemopt_bubble_toolbar,
            text="Save HD Image",
            command=self.save_chemopt_bubble_figure_image,
        ).grid(row=0, column=2, sticky="e")
        self.chemopt_bubble_canvas = tk.Canvas(chemopt_bubble_frame, background="#ffffff", highlightthickness=0)
        self.chemopt_bubble_canvas.grid(row=1, column=0, sticky="nsew")
        chemopt_bubble_y = ttk.Scrollbar(
            chemopt_bubble_frame,
            orient="vertical",
            command=self.chemopt_bubble_canvas.yview,
        )
        chemopt_bubble_x = ttk.Scrollbar(
            chemopt_bubble_frame,
            orient="horizontal",
            command=self.chemopt_bubble_canvas.xview,
        )
        self.chemopt_bubble_canvas.configure(yscrollcommand=chemopt_bubble_y.set, xscrollcommand=chemopt_bubble_x.set)
        chemopt_bubble_y.grid(row=1, column=1, sticky="ns")
        chemopt_bubble_x.grid(row=2, column=0, sticky="ew")

    def _build_idt_converter_tool(self, parent) -> None:
        tk = self.tk
        ttk = self.ttk

        parent.columnconfigure(0, weight=0)
        parent.columnconfigure(1, weight=1)
        parent.rowconfigure(0, weight=1)
        parent.rowconfigure(1, weight=0)

        left = ttk.Frame(parent, padding=12)
        left.grid(row=0, column=0, sticky="ns")
        left.columnconfigure(1, weight=1)

        right = ttk.Frame(parent, padding=(0, 12, 12, 12))
        right.grid(row=0, column=1, sticky="nsew")
        right.rowconfigure(1, weight=1)
        right.columnconfigure(0, weight=1)

        credit = ttk.Label(parent, text="By Alexander Apkarian", font=("TkDefaultFont", 11), foreground="#6b7280")
        credit.grid(row=1, column=0, columnspan=2, sticky="e", padx=12, pady=(0, 6))

        row = 0
        nav_row = ttk.Frame(left)
        nav_row.grid(row=row, column=0, columnspan=2, sticky="ew", pady=(0, 8))
        nav_row.columnconfigure(1, weight=1)
        ttk.Button(nav_row, text="Home", command=self._show_home).grid(row=0, column=0, sticky="w", padx=(0, 8))
        ttk.Label(nav_row, text="IDT notation converter", font=("TkDefaultFont", 14, "bold")).grid(
            row=0,
            column=1,
            sticky="e",
        )
        row += 1

        self._combo(
            left,
            row,
            "ASO Chemistry",
            "converter_aso_chemistry",
            list(CHEMISTRY_PRESETS),
            DEFAULT_CHEMISTRY,
        )
        self.vars["converter_aso_chemistry"].trace_add("write", lambda *_: self._on_converter_chemistry_changed())
        row += 1
        self.converter_custom_button = ttk.Button(
            left,
            text="Edit Custom Pattern",
            command=self._open_converter_custom_chemistry_editor,
        )
        self.converter_custom_button.grid(row=row, column=1, sticky="ew", pady=3)
        row += 1
        self.converter_gap_entry = self._entry(left, row, "Central Gap Length", "converter_gap_length", "12")
        row += 1
        self.converter_wing_entry = self._entry(left, row, "Modified Wing Length", "converter_wing_length", "3")
        row += 1
        self.converter_wing_combo = self._combo(
            left,
            row,
            "Modified Wing Chemistry",
            "converter_wing_chemistry",
            list(RIBOSE_MODIFICATION_OPTIONS) + ["None"],
            "LNA",
        )
        row += 1
        self.converter_backbone_combo = self._combo(
            left,
            row,
            "Backbone Modification",
            "converter_backbone_modification",
            ["PS", "PO", "MIXED", "None"],
            "PS",
        )
        row += 1

        ttk.Separator(left).grid(row=row, column=0, columnspan=2, sticky="ew", pady=10)
        row += 1
        ttk.Label(left, text="ASO sequences 5' to 3'").grid(row=row, column=0, columnspan=2, sticky="w")
        row += 1
        ttk.Label(
            left,
            text="Enter one clean sequence per line. Hyphens, spaces, U, and T are tolerated.",
            foreground="#4b5563",
            wraplength=330,
        ).grid(row=row, column=0, columnspan=2, sticky="w", pady=(0, 4))
        row += 1
        self.idt_sequence_text = tk.Text(left, width=44, height=14, wrap="none", font=("Menlo", 11))
        self.idt_sequence_text.grid(row=row, column=0, columnspan=2, sticky="nsew")
        left.rowconfigure(row, weight=1)
        row += 1

        converter_button_row = ttk.Frame(left)
        converter_button_row.grid(row=row, column=0, columnspan=2, sticky="ew", pady=(10, 0))
        converter_button_row.columnconfigure(0, weight=1)
        converter_button_row.columnconfigure(1, weight=1)
        ttk.Button(converter_button_row, text="Convert to IDT Notation", command=self.convert_idt_notation).grid(
            row=0,
            column=0,
            sticky="ew",
            padx=(0, 6),
        )
        ttk.Button(converter_button_row, text="Export Excel", command=self.export_idt_converter_excel).grid(
            row=0,
            column=1,
            sticky="ew",
        )

        ttk.Label(right, text="IDT notation output", font=("TkDefaultFont", 14, "bold")).grid(
            row=0,
            column=0,
            sticky="w",
            pady=(0, 8),
        )
        output_frame = ttk.Frame(right)
        output_frame.grid(row=1, column=0, sticky="nsew")
        output_frame.rowconfigure(0, weight=1)
        output_frame.columnconfigure(0, weight=1)
        self.idt_output_text = tk.Text(output_frame, wrap="none", font=("Menlo", 11), undo=False)
        self.idt_output_text.configure(tabs=("70", "260", "520"))
        output_y = ttk.Scrollbar(output_frame, orient="vertical", command=self.idt_output_text.yview)
        output_x = ttk.Scrollbar(output_frame, orient="horizontal", command=self.idt_output_text.xview)
        self.idt_output_text.configure(yscrollcommand=output_y.set, xscrollcommand=output_x.set)
        self.idt_output_text.grid(row=0, column=0, sticky="nsew")
        output_y.grid(row=0, column=1, sticky="ns")
        output_x.grid(row=1, column=0, sticky="ew")

    def _build_penalty_design_tool(self, parent) -> None:
        tk = self.tk
        ttk = self.ttk

        parent.columnconfigure(0, weight=0)
        parent.columnconfigure(1, weight=1)
        parent.rowconfigure(0, weight=1)
        parent.rowconfigure(1, weight=0)

        left = ttk.Frame(parent, padding=12)
        left.grid(row=0, column=0, sticky="ns")
        left.columnconfigure(1, weight=1)

        right = ttk.Frame(parent, padding=(0, 12, 12, 12))
        right.grid(row=0, column=1, sticky="nsew")
        right.rowconfigure(1, weight=1)
        right.columnconfigure(0, weight=1)

        credit = ttk.Label(parent, text="By Alexander Apkarian", font=("TkDefaultFont", 11), foreground="#6b7280")
        credit.grid(row=1, column=0, columnspan=2, sticky="e", padx=12, pady=(0, 6))

        row = 0
        nav_row = ttk.Frame(left)
        nav_row.grid(row=row, column=0, columnspan=2, sticky="ew", pady=(0, 8))
        nav_row.columnconfigure(1, weight=1)
        ttk.Button(nav_row, text="Home", command=self._show_home).grid(row=0, column=0, sticky="w", padx=(0, 8))
        ttk.Label(nav_row, text="Penalty ASO Design Tool", font=("TkDefaultFont", 14, "bold")).grid(
            row=0,
            column=1,
            sticky="e",
        )
        row += 1

        self._entry(left, row, "Target Gene Name", "penalty_target_gene")
        row += 1
        self._entry(left, row, "Target Identifier", "penalty_target_identifier")
        row += 1
        self._entry(left, row, "Chemistry Identifier", "penalty_chemistry_number")
        row += 1

        self._combo(
            left,
            row,
            "ASO Chemistry",
            "penalty_aso_chemistry",
            list(CHEMISTRY_PRESETS),
            DEFAULT_CHEMISTRY,
        )
        self.vars["penalty_aso_chemistry"].trace_add("write", lambda *_: self._on_penalty_chemistry_changed())
        row += 1
        self.penalty_custom_button = ttk.Button(
            left,
            text="Edit Custom Pattern",
            command=self._open_penalty_custom_chemistry_editor,
        )
        self.penalty_custom_button.grid(row=row, column=1, sticky="ew", pady=3)
        row += 1
        self.penalty_gap_entry = self._entry(left, row, "Central Gap Length", "penalty_gap_length", "12")
        row += 1
        self.penalty_wing_entry = self._entry(left, row, "Modified Wing Length", "penalty_wing_length", "3")
        row += 1
        self.penalty_wing_combo = self._combo(
            left,
            row,
            "Modified Wing Chemistry",
            "penalty_wing_chemistry",
            list(RIBOSE_MODIFICATION_OPTIONS) + ["None"],
            "LNA",
        )
        row += 1
        self.penalty_backbone_combo = self._combo(
            left,
            row,
            "Backbone Modification",
            "penalty_backbone_modification",
            ["PS", "PO", "MIXED", "None"],
            "PS",
        )
        row += 1

        ttk.Separator(left).grid(row=row, column=0, columnspan=2, sticky="ew", pady=10)
        row += 1
        self._entry(left, row, "First Parent ASO Start Position (RNA 3' to 5', first base = 0)", "penalty_parent_start", "0")
        row += 1
        self._entry(left, row, "Number of Parent ASOs", "penalty_parent_count", "3")
        row += 1
        self._entry(left, row, "Microwalk Step Size", "penalty_step_size", "1")
        row += 1
        self._combo(
            left,
            row,
            "Penalty Positions",
            "penalty_position_mode",
            list(PENALTY_POSITION_MODES),
            "Central core positions",
        )
        row += 1
        self._entry(left, row, "Selected ASO Positions", "penalty_selected_positions")
        row += 1
        self._combo(
            left,
            row,
            "Penalty Base Options",
            "penalty_base_mode",
            list(PENALTY_BASE_MODES),
            "Like-for-like only",
        )
        row += 1

        ttk.Separator(left).grid(row=row, column=0, columnspan=2, sticky="ew", pady=10)
        row += 1
        ttk.Label(left, text="RNA 5' to 3' (hyphens optional)").grid(row=row, column=0, columnspan=2, sticky="w")
        row += 1
        self.penalty_sequence_text = tk.Text(left, width=44, height=10, wrap="word", font=("Menlo", 11))
        self.penalty_sequence_text.grid(row=row, column=0, columnspan=2, sticky="nsew")
        left.rowconfigure(row, weight=1)
        row += 1

        penalty_button_row = ttk.Frame(left)
        penalty_button_row.grid(row=row, column=0, columnspan=2, sticky="ew", pady=(10, 0))
        penalty_button_row.columnconfigure(0, weight=1)
        penalty_button_row.columnconfigure(1, weight=1)
        ttk.Button(penalty_button_row, text="Generate Penalty ASOs", command=self.generate_penalty_design).grid(
            row=0,
            column=0,
            sticky="ew",
            padx=(0, 6),
        )
        ttk.Button(penalty_button_row, text="Export Excel", command=self.export_penalty_design_excel).grid(
            row=0,
            column=1,
            sticky="ew",
        )

        ttk.Label(right, text="Penalty ASO output", font=("TkDefaultFont", 14, "bold")).grid(
            row=0,
            column=0,
            sticky="w",
            pady=(0, 8),
        )
        penalty_notebook = ttk.Notebook(right)
        penalty_notebook.grid(row=1, column=0, sticky="nsew")

        output_frame = ttk.Frame(penalty_notebook)
        output_frame.rowconfigure(0, weight=1)
        output_frame.columnconfigure(0, weight=1)
        penalty_notebook.add(output_frame, text="Output Table")
        self.penalty_output_text = tk.Text(output_frame, wrap="none", font=("Menlo", 11), undo=False)
        self.penalty_output_text.tag_configure("priority_recommended", foreground="#006100")
        self.penalty_output_text.tag_configure("priority_alternative", foreground="#9a6700")
        self.penalty_output_text.tag_configure("priority_lower", foreground="#c00000")
        output_y = ttk.Scrollbar(output_frame, orient="vertical", command=self.penalty_output_text.yview)
        output_x = ttk.Scrollbar(output_frame, orient="horizontal", command=self.penalty_output_text.xview)
        self.penalty_output_text.configure(yscrollcommand=output_y.set, xscrollcommand=output_x.set)
        self.penalty_output_text.grid(row=0, column=0, sticky="nsew")
        output_y.grid(row=0, column=1, sticky="ns")
        output_x.grid(row=1, column=0, sticky="ew")

        penalty_bubble_frame = ttk.Frame(penalty_notebook)
        penalty_bubble_frame.rowconfigure(1, weight=1)
        penalty_bubble_frame.columnconfigure(0, weight=1)
        penalty_notebook.add(penalty_bubble_frame, text="Visualisation")
        penalty_bubble_toolbar = ttk.Frame(penalty_bubble_frame, padding=(0, 0, 0, 6))
        penalty_bubble_toolbar.grid(row=0, column=0, columnspan=2, sticky="ew")
        penalty_bubble_toolbar.columnconfigure(0, weight=1)
        ttk.Checkbutton(
            penalty_bubble_toolbar,
            text="Show base letters",
            variable=self.show_penalty_bubble_base_letters,
            command=self._refresh_penalty_bubble_figure,
        ).grid(row=0, column=1, sticky="e", padx=(0, 12))
        ttk.Button(
            penalty_bubble_toolbar,
            text="Save HD Image",
            command=self.save_penalty_bubble_figure_image,
        ).grid(row=0, column=2, sticky="e")
        self.penalty_bubble_canvas = tk.Canvas(penalty_bubble_frame, background="#ffffff", highlightthickness=0)
        self.penalty_bubble_canvas.grid(row=1, column=0, sticky="nsew")
        penalty_bubble_y = ttk.Scrollbar(
            penalty_bubble_frame,
            orient="vertical",
            command=self.penalty_bubble_canvas.yview,
        )
        penalty_bubble_x = ttk.Scrollbar(
            penalty_bubble_frame,
            orient="horizontal",
            command=self.penalty_bubble_canvas.xview,
        )
        self.penalty_bubble_canvas.configure(yscrollcommand=penalty_bubble_y.set, xscrollcommand=penalty_bubble_x.set)
        penalty_bubble_y.grid(row=1, column=1, sticky="ns")
        penalty_bubble_x.grid(row=2, column=0, sticky="ew")

    def _entry(self, parent, row: int, label: str, name: str, value: str = ""):
        ttk = self.ttk
        ttk.Label(parent, text=label).grid(row=row, column=0, sticky="w", pady=3)
        entry = ttk.Entry(parent, textvariable=self._var(name, value), width=24)
        entry.grid(row=row, column=1, sticky="ew", pady=3)
        return entry

    def _refresh_bubble_figure(self) -> None:
        if self.last_result is not None:
            self._variant_bubble_dirty = True
            self._render_variant_bubble_if_visible()

    def _variant_visualisation_is_selected(self) -> bool:
        notebook = getattr(self, "variant_notebook", None)
        bubble_frame = getattr(self, "variant_bubble_frame", None)
        if notebook is None or bubble_frame is None:
            return False
        try:
            return notebook.select() == str(bubble_frame)
        except Exception:
            return False

    def _render_variant_bubble_if_visible(self) -> None:
        if not self._variant_bubble_dirty or self.last_result is None:
            return
        if not self._variant_visualisation_is_selected():
            return
        self._variant_bubble_dirty = False
        self._render_bubble_figure(self.last_result)

    def _combo(self, parent, row: int, label: str, name: str, values: list[str], value: str):
        ttk = self.ttk
        ttk.Label(parent, text=label).grid(row=row, column=0, sticky="w", pady=3)
        combo = ttk.Combobox(parent, textvariable=self._var(name, value), values=values, state="readonly", width=22)
        combo.grid(row=row, column=1, sticky="ew", pady=3)
        return combo

    def _initial_custom_pattern(self) -> tuple[list[str], list[str]]:
        if self.custom_base_modifications:
            return list(self.custom_base_modifications), list(self.custom_linkages)
        try:
            base_mods, linkages = build_gapmer_pattern(
                int(self.vars["gap_length"].get()),
                int(self.vars["wing_length"].get()),
                self.vars["wing_chemistry"].get(),
                self.vars["backbone_modification"].get(),
            )
        except Exception:
            base_mods, linkages = build_gapmer_pattern(12, 3, "LNA", "PS")
        return list(base_mods), list(linkages)

    def _set_custom_pattern(self, base_mods: list[str], linkages: list[str]) -> None:
        self.custom_base_modifications = tuple(base_mods)
        self.custom_linkages = tuple(linkages)
        self.vars["gap_length"].set(str(len(base_mods)))
        self.vars["wing_length"].set("0")
        self.vars["wing_chemistry"].set("Custom")
        linkage_set = set(linkages)
        if linkage_set == {"PS"}:
            backbone = "PS"
        elif linkage_set <= {"PO"}:
            backbone = "PO"
        else:
            backbone = "mixed PS/PO"
        self.vars["backbone_modification"].set(backbone)

    def _open_custom_chemistry_editor(
        self,
        *,
        chemistry_var_name: str = "aso_chemistry",
        initial_pattern_func=None,
        apply_func=None,
        window_attr: str = "custom_editor_window",
        title: str = "Custom Chemistry Pattern",
        allow_motif_selection: bool = False,
        initial_motif_positions: tuple[int, ...] = (),
        require_custom: bool = True,
        template_callback=None,
    ) -> None:
        if require_custom and self.vars[chemistry_var_name].get() != "Custom":
            return
        editor_window = getattr(self, window_attr, None)
        if editor_window is not None and editor_window.winfo_exists():
            editor_window.lift()
            editor_window.focus_force()
            return

        tk = self.tk
        ttk = self.ttk
        if initial_pattern_func is None:
            initial_pattern_func = self._initial_custom_pattern
        if apply_func is None:
            apply_func = self._set_custom_pattern
        base_mods, linkages = initial_pattern_func()

        win = tk.Toplevel(self.root)
        setattr(self, window_attr, win)
        win.title(title)
        win.geometry("980x430")
        win.minsize(760, 360)
        win.transient(self.root)
        win.columnconfigure(0, weight=1)
        win.rowconfigure(2, weight=1)

        def close_window() -> None:
            setattr(self, window_attr, None)
            win.destroy()

        win.protocol("WM_DELETE_WINDOW", close_window)

        top = ttk.Frame(win, padding=(12, 12, 12, 6))
        top.grid(row=0, column=0, sticky="ew")
        top.columnconfigure(11, weight=1)

        length_var = tk.StringVar(value=str(len(base_mods)))
        ribose_tool = tk.StringVar(value="DNA")
        nucleobase_tool = tk.StringVar(value="None")
        linkage_tool = tk.StringVar(value="PS")
        motif_mode = tk.BooleanVar(value=False)
        motif_positions = set(int(position) for position in initial_motif_positions)
        template_metadata = None

        ttk.Label(top, text="ASO length").grid(row=0, column=0, sticky="w", padx=(0, 6))
        length_spin = ttk.Spinbox(top, from_=1, to=80, textvariable=length_var, width=6)
        length_spin.grid(row=0, column=1, sticky="w", padx=(0, 18))
        ttk.Label(top, text="Ribose").grid(row=0, column=2, sticky="w", padx=(0, 6))
        ttk.Combobox(top, textvariable=ribose_tool, values=list(RIBOSE_MODIFICATION_OPTIONS), state="readonly", width=10).grid(
            row=0, column=3, sticky="w", padx=(0, 18)
        )
        ttk.Label(top, text="Base").grid(row=0, column=4, sticky="w", padx=(0, 6))
        ttk.Combobox(
            top,
            textvariable=nucleobase_tool,
            values=list(NUCLEOBASE_MODIFICATION_OPTIONS),
            state="readonly",
            width=8,
        ).grid(row=0, column=5, sticky="w", padx=(0, 18))
        ttk.Label(top, text="Linkage").grid(row=0, column=6, sticky="w", padx=(0, 6))
        ttk.Combobox(top, textvariable=linkage_tool, values=list(LINKAGE_OPTIONS), state="readonly", width=8).grid(
            row=0, column=7, sticky="w"
        )
        if allow_motif_selection:
            ttk.Checkbutton(top, text="Select walk motif", variable=motif_mode).grid(
                row=0,
                column=8,
                sticky="w",
                padx=(18, 0),
            )
            instruction = (
                "Choose ribose, base, or linkage options, then click circles or diamonds to edit the starting chemistry. "
                "Turn on Select walk motif, then click the core bases that should slide along the central gap."
            )
        else:
            instruction = (
                "Choose ribose, base, or linkage options, then click individual circles or diamonds below to apply them. "
                "Use a template as a starting point if helpful."
            )
        ttk.Label(
            top,
            text=instruction,
            foreground="#4b5563",
            wraplength=780,
        ).grid(row=1, column=0, columnspan=10, sticky="w", pady=(8, 0))

        template_row = ttk.Frame(win, padding=(12, 0, 12, 6))
        template_row.grid(row=1, column=0, sticky="ew")

        canvas_frame = ttk.Frame(win, padding=(12, 0, 12, 8))
        canvas_frame.grid(row=2, column=0, sticky="nsew")
        canvas_frame.rowconfigure(0, weight=1)
        canvas_frame.columnconfigure(0, weight=1)
        canvas = tk.Canvas(canvas_frame, height=220, background="#ffffff", highlightthickness=1, highlightbackground="#d1d5db")
        canvas.grid(row=0, column=0, sticky="nsew")
        canvas_x = ttk.Scrollbar(canvas_frame, orient="horizontal", command=canvas.xview)
        canvas_x.grid(row=1, column=0, sticky="ew")
        canvas.configure(xscrollcommand=canvas_x.set)

        bottom = ttk.Frame(win, padding=(12, 0, 12, 12))
        bottom.grid(row=3, column=0, sticky="ew")
        bottom.columnconfigure(0, weight=1)

        ribose_colors = {
            "DNA": "#f8fafc",
            "LNA": "#f4cccc",
            "MOE": "#cfe2f3",
            "2'OMe": "#d9ead3",
            "2'F": "#d9d2e9",
        }
        ribose_text = {
            "DNA": "DNA",
            "LNA": "LNA",
            "MOE": "MOE",
            "2'OMe": "OMe",
            "2'F": "F",
        }

        def popup_ribose_mod(modification: str) -> str:
            ribose_mod, _base_mod = split_base_modification(modification)
            return ribose_mod or "DNA"

        def popup_has_5mec(modification: str) -> bool:
            _ribose_mod, base_mod = split_base_modification(modification)
            return base_mod == "5MeC"

        def normalise_lengths() -> None:
            nonlocal linkages
            needed = max(0, len(base_mods) - 1)
            if len(linkages) > needed:
                linkages = linkages[:needed]
            while len(linkages) < needed:
                linkages.append(linkages[-1] if linkages else "PS")

        def mark_template_custom() -> None:
            if template_metadata is not None:
                template_metadata["label"] = "Custom"

        def resize_from_length(*_args) -> None:
            nonlocal template_metadata
            try:
                new_length = max(1, min(80, int(length_var.get())))
            except Exception:
                new_length = len(base_mods)
            old_length = len(base_mods)
            if new_length > old_length:
                base_mods.extend(["DNA"] * (new_length - old_length))
            elif new_length < old_length:
                del base_mods[new_length:]
                motif_positions.intersection_update(range(new_length))
            normalise_lengths()
            if template_metadata is not None:
                template_metadata = {
                    "label": "Custom",
                    "gap_length": len(base_mods),
                    "wing_length": 0,
                    "wing_chemistry": "Custom",
                    "backbone_modification": summarise_linkages(tuple(linkages)),
                }
            length_var.set(str(len(base_mods)))
            draw()

        def set_gapmer(label: str, gap: int, wing: int, wing_chemistry: str, backbone: str, **kwargs) -> None:
            nonlocal base_mods, linkages, template_metadata
            base_tuple, linkage_tuple = build_gapmer_pattern(gap, wing, wing_chemistry, backbone, **kwargs)
            base_mods = list(base_tuple)
            linkages = list(linkage_tuple)
            motif_positions.clear()
            template_metadata = {
                "label": label,
                "gap_length": gap,
                "wing_length": wing,
                "wing_chemistry": wing_chemistry,
                "backbone_modification": summarise_linkages(tuple(linkages)),
            }
            length_var.set(str(len(base_mods)))
            draw()

        def set_all(label: str, base_mod: str, linkage: str) -> None:
            nonlocal base_mods, linkages, template_metadata
            try:
                length = max(1, min(80, int(length_var.get())))
            except Exception:
                length = len(base_mods)
            base_mods = [base_mod] * length
            linkages = [linkage] * max(0, length - 1)
            motif_positions.clear()
            template_metadata = {
                "label": label,
                "gap_length": length,
                "wing_length": 0,
                "wing_chemistry": "Custom",
                "backbone_modification": summarise_linkages(tuple(linkages)),
            }
            length_var.set(str(length))
            draw()

        template_buttons = [
            ("5-10-5 MOE", lambda: set_gapmer("5-10-5 MOE/DNA", 10, 5, "MOE", "PS")),
            (
                "KT777/valeriasen",
                lambda: set_gapmer(
                    "KT777/valeriasen",
                    10,
                    5,
                    "MOE",
                    "MIXED",
                    methyl_c=True,
                    linkage_pattern="KT777",
                ),
            ),
            ("4-10-4 LNA", lambda: set_gapmer("4-10-4 LNA/DNA", 10, 4, "LNA", "PS")),
            ("3-12-3 LNA", lambda: set_gapmer("3-12-3 LNA/DNA", 12, 3, "LNA", "PS")),
            ("3-10-3 LNA", lambda: set_gapmer("3-10-3 LNA/DNA", 10, 3, "LNA", "PS")),
            ("3-9-3 LNA", lambda: set_gapmer("3-9-3 LNA/DNA", 9, 3, "LNA", "PS")),
            ("All DNA", lambda: set_all("Custom", "DNA", "PS")),
        ]
        for idx, (label, command) in enumerate(template_buttons):
            ttk.Button(template_row, text=label, command=command).grid(
                row=idx // 4,
                column=idx % 4,
                sticky="ew",
                padx=(0, 6),
                pady=3,
            )

        def draw() -> None:
            canvas.delete("all")
            normalise_lengths()
            step = 54
            radius = 18
            start_x = 48
            y = 102
            total_width = start_x * 2 + max(0, len(base_mods) - 1) * step

            for idx, linkage in enumerate(linkages):
                x = start_x + idx * step + step / 2
                fill = "#111827" if linkage == "PS" else "#ffffff"
                text_fill = "#ffffff" if linkage == "PS" else "#374151"
                tag = f"linkage_{idx}"
                canvas.create_polygon(
                    x,
                    y - 12,
                    x + 12,
                    y,
                    x,
                    y + 12,
                    x - 12,
                    y,
                    fill=fill,
                    outline="#6b7280",
                    width=1.5,
                    tags=(tag,),
                )
                canvas.create_text(x, y, text=linkage, font=("TkDefaultFont", 8, "bold"), fill=text_fill, tags=(tag,))
                canvas.tag_bind(tag, "<Button-1>", lambda _event, i=idx: set_linkage(i))

            for idx, mod in enumerate(base_mods):
                x = start_x + idx * step
                tag = f"base_{idx}"
                selected = idx in motif_positions
                ribose_mod = popup_ribose_mod(mod)
                has_5mec = popup_has_5mec(mod)
                canvas.create_oval(
                    x - radius,
                    y - radius,
                    x + radius,
                    y + radius,
                    fill=ribose_colors.get(ribose_mod, "#f8fafc"),
                    outline="#1f4e79" if selected else "#374151",
                    width=3 if selected else 1.5,
                    tags=(tag,),
                )
                canvas.create_text(
                    x,
                    y,
                    text=ribose_text.get(ribose_mod, ribose_mod),
                    font=("TkDefaultFont", 9, "bold"),
                    fill="#111827",
                    tags=(tag,),
                )
                if has_5mec:
                    canvas.create_oval(
                        x + 8,
                        y - 17,
                        x + 18,
                        y - 7,
                        fill="#f6c343",
                        outline="#111827",
                        tags=(tag,),
                    )
                if selected:
                    canvas.create_text(x, y - 31, text="walk", font=("TkDefaultFont", 8, "bold"), fill="#1f4e79")
                canvas.create_text(x, y + 34, text=str(idx + 1), font=("TkDefaultFont", 9), fill="#6b7280")
                canvas.tag_bind(tag, "<Button-1>", lambda _event, i=idx: set_base(i))

            canvas.configure(scrollregion=(0, 0, total_width + 40, 210))

        def set_base(index: int) -> None:
            if allow_motif_selection and motif_mode.get():
                if index in motif_positions:
                    motif_positions.remove(index)
                else:
                    motif_positions.add(index)
            else:
                base_mods[index] = combine_base_modification(ribose_tool.get(), nucleobase_tool.get())
                mark_template_custom()
            draw()

        def set_linkage(index: int) -> None:
            linkages[index] = linkage_tool.get()
            mark_template_custom()
            draw()

        def apply_pattern() -> None:
            if template_callback is not None and template_metadata is not None:
                template_callback(**template_metadata)
            if allow_motif_selection:
                apply_func(base_mods, linkages, tuple(sorted(motif_positions)))
            else:
                apply_func(base_mods, linkages)
            close_window()

        length_spin.configure(command=resize_from_length)
        length_spin.bind("<Return>", resize_from_length)
        length_spin.bind("<FocusOut>", resize_from_length)

        ttk.Button(bottom, text="Apply", command=apply_pattern).grid(row=0, column=1, sticky="e", padx=(0, 6))
        ttk.Button(bottom, text="Cancel", command=close_window).grid(row=0, column=2, sticky="e")

        draw()

    def _on_chemistry_changed(self) -> None:
        self._apply_chemistry_preset()
        if self.vars["aso_chemistry"].get() == "Custom":
            self.root.after(80, self._open_custom_chemistry_editor)

    def _apply_chemistry_preset(self) -> None:
        label = self.vars["aso_chemistry"].get()
        preset = CHEMISTRY_PRESETS.get(label)
        custom = preset is None
        if preset is not None:
            self.vars["gap_length"].set(str(preset["gap_length"]))
            self.vars["wing_length"].set(str(preset["wing_length"]))
            self.vars["wing_chemistry"].set(str(preset["wing_chemistry"]))
            self.vars["backbone_modification"].set(str(preset["backbone_modification"]))
            self.custom_base_modifications = ()
            self.custom_linkages = ()
        state = "disabled"
        combo_state = "disabled"
        self.gap_entry.configure(state=state)
        self.wing_entry.configure(state=state)
        self.wing_combo.configure(state=combo_state)
        self.backbone_combo.configure(state=combo_state)
        self.custom_button.configure(state="normal" if custom else "disabled")

    def _initial_chemopt_pattern(self) -> tuple[list[str], list[str]]:
        if self.chemopt_base_modifications:
            return list(self.chemopt_base_modifications), list(self.chemopt_linkages)
        try:
            chemistry = resolve_chemistry(self._collect_chemopt_chemistry_inputs())
            return list(chemistry.base_modifications), list(chemistry.linkages)
        except Exception:
            base_mods, linkages = build_gapmer_pattern(10, 5, "MOE", "PS")
            return list(base_mods), list(linkages)

    def _set_chemopt_pattern(
        self,
        base_mods: list[str],
        linkages: list[str],
        motif_positions: tuple[int, ...],
    ) -> None:
        self.chemopt_base_modifications = tuple(base_mods)
        self.chemopt_linkages = tuple(linkages)
        self.chemopt_motif_positions = tuple(sorted(motif_positions))
        try:
            wing_len = int(self.vars["chemopt_wing_length"].get())
            gap_len = int(self.vars["chemopt_gap_length"].get())
        except Exception:
            wing_len = 0
            gap_len = 0
        if len(base_mods) != wing_len * 2 + gap_len:
            self.vars["chemopt_gap_length"].set(str(len(base_mods)))
            self.vars["chemopt_wing_length"].set("0")
            self.vars["chemopt_wing_chemistry"].set("Custom")
        if not self._chemopt_pattern_matches_label(self.vars["chemopt_aso_chemistry"].get(), base_mods, linkages):
            self.vars["chemopt_aso_chemistry"].set("Custom")
        linkage_set = set(linkages)
        if linkage_set == {"PS"}:
            backbone = "PS"
        elif linkage_set <= {"PO"}:
            backbone = "PO"
        else:
            backbone = "mixed PS/PO"
        self.vars["chemopt_backbone_modification"].set(backbone)
        self._update_chemopt_motif_status()

    def _chemopt_pattern_matches_label(self, label: str, base_mods: list[str], linkages: list[str]) -> bool:
        preset = CHEMISTRY_PRESETS.get(label)
        if preset is None:
            return False
        try:
            preset_base_mods, preset_linkages = build_gapmer_pattern(
                int(preset["gap_length"]),
                int(preset["wing_length"]),
                str(preset["wing_chemistry"]),
                str(preset["backbone_modification"]),
                methyl_c=bool(preset.get("methyl_c", False)),
                linkage_pattern=str(preset.get("linkage_pattern", "")),
            )
        except Exception:
            return False
        return tuple(base_mods) == preset_base_mods and tuple(linkages) == preset_linkages

    def _set_chemopt_template_metadata(
        self,
        label: str,
        gap_length: int,
        wing_length: int,
        wing_chemistry: str,
        backbone_modification: str,
    ) -> None:
        self.vars["chemopt_aso_chemistry"].set(label if label in CHEMISTRY_PRESETS else "Custom")
        self.vars["chemopt_gap_length"].set(str(gap_length))
        self.vars["chemopt_wing_length"].set(str(wing_length))
        self.vars["chemopt_wing_chemistry"].set(str(wing_chemistry))
        self.vars["chemopt_backbone_modification"].set(str(backbone_modification))

    def _open_chemopt_chemistry_editor(self) -> None:
        self._open_custom_chemistry_editor(
            chemistry_var_name="chemopt_aso_chemistry",
            initial_pattern_func=self._initial_chemopt_pattern,
            apply_func=self._set_chemopt_pattern,
            window_attr="chemopt_editor_window",
            title="Chemistry Optimisation Pattern",
            allow_motif_selection=True,
            initial_motif_positions=self.chemopt_motif_positions,
            require_custom=False,
            template_callback=self._set_chemopt_template_metadata,
        )

    def _on_chemopt_chemistry_changed(self) -> None:
        self._apply_chemopt_chemistry_preset()
        if self.vars["chemopt_aso_chemistry"].get() == "Custom":
            self.root.after(80, self._open_chemopt_chemistry_editor)

    def _apply_chemopt_chemistry_preset(self) -> None:
        label = self.vars["chemopt_aso_chemistry"].get()
        preset = CHEMISTRY_PRESETS.get(label)
        if preset is not None:
            self.vars["chemopt_gap_length"].set(str(preset["gap_length"]))
            self.vars["chemopt_wing_length"].set(str(preset["wing_length"]))
            self.vars["chemopt_wing_chemistry"].set(str(preset["wing_chemistry"]))
            self.vars["chemopt_backbone_modification"].set(str(preset["backbone_modification"]))
            self.chemopt_base_modifications = ()
            self.chemopt_linkages = ()
            self.chemopt_motif_positions = ()
        if hasattr(self, "chemopt_custom_button"):
            self.chemopt_custom_button.configure(state="normal")
        self._update_chemopt_motif_status()

    def _collect_chemopt_chemistry_inputs(self) -> AsoInputs:
        return AsoInputs(
            aso_chemistry=self.vars["chemopt_aso_chemistry"].get(),
            gap_length=int(self.vars["chemopt_gap_length"].get()),
            wing_length=int(self.vars["chemopt_wing_length"].get()),
            wing_chemistry=self.vars["chemopt_wing_chemistry"].get(),
            backbone_modification=self.vars["chemopt_backbone_modification"].get(),
            custom_base_modifications=(
                self.chemopt_base_modifications if self.vars["chemopt_aso_chemistry"].get() == "Custom" else ()
            ),
            custom_linkages=(
                self.chemopt_linkages if self.vars["chemopt_aso_chemistry"].get() == "Custom" else ()
            ),
        )

    def _collect_chemopt_pattern(self):
        wing_len = int(self.vars["chemopt_wing_length"].get())
        gap_len = int(self.vars["chemopt_gap_length"].get())
        if self.chemopt_base_modifications:
            base_mods = self.chemopt_base_modifications
            linkages = self.chemopt_linkages
            label = (
                f"{self.vars['chemopt_aso_chemistry'].get()}, edited chemistry walk pattern, "
                f"{summarise_linkages(linkages)} backbone modification"
            )
            return label, base_mods, linkages, wing_len, gap_len

        chemistry = resolve_chemistry(self._collect_chemopt_chemistry_inputs())
        return (
            chemistry_display_label(chemistry),
            chemistry.base_modifications,
            chemistry.linkages,
            chemistry.wing_length,
            chemistry.gap_length,
        )

    def _update_chemopt_motif_status(self) -> None:
        if not hasattr(self, "chemopt_motif_status"):
            return
        if not self.chemopt_motif_positions:
            self.chemopt_motif_status.configure(text="Walk motif: not selected")
            return
        aso_positions = tuple(position + 1 for position in self.chemopt_motif_positions)
        try:
            wing_len = int(self.vars["chemopt_wing_length"].get())
            gap_positions = tuple(position - wing_len + 1 for position in self.chemopt_motif_positions)
            gap_positions = tuple(position for position in gap_positions if position > 0)
        except Exception:
            gap_positions = ()
        gap_text = f"; core positions {self._format_position_tuple(gap_positions)}" if gap_positions else ""
        self.chemopt_motif_status.configure(
            text=f"Walk motif: ASO positions {self._format_position_tuple(aso_positions)}{gap_text}"
        )

    def generate_chemistry_optimization(self) -> None:
        try:
            label, base_mods, linkages, wing_len, gap_len = self._collect_chemopt_pattern()
            try:
                step_size = int(self.vars["chemopt_step_size"].get())
            except Exception:
                raise AsoInputError("Microwalk step size must be a whole number.")
            rows = chemistry_optimization_walk(
                self.chemopt_sequence_text.get("1.0", "end"),
                base_mods,
                linkages,
                wing_len,
                gap_len,
                self.chemopt_motif_positions,
                step_size,
            )
        except Exception as exc:
            self.messagebox.showerror(APP_NAME, str(exc))
            return
        self.last_chemopt_rows = rows
        self.last_chemopt_chemistry_label = label
        self._render_chemopt_output(label, rows)
        self._render_chemopt_bubble_figure(rows)

    def _render_chemopt_output(self, chemistry_label: str, rows) -> None:
        text = self.chemopt_output_text
        text.configure(state="normal")
        text.delete("1.0", "end")
        text.configure(tabs=self._chemopt_tab_stops(rows))
        text.insert("end", f"Base chemistry\t{chemistry_label}\n")
        text.insert("end", f"Variants generated\t{len(rows)}\n\n")
        text.insert(
            "end",
            "Variant\tMotif ASO positions\tMotif core positions\tSequence 5' to 3'\tChemistry pattern\tIDT notation\n",
        )
        for row in rows:
            text.insert(
                "end",
                (
                    f"{row.row_number}\t"
                    f"{self._format_position_tuple(row.motif_aso_positions)}\t"
                    f"{self._format_position_tuple(row.motif_gap_positions)}\t"
                    f"{row.clean_sequence}\t"
                    f"{self._format_base_mod_pattern(row.base_modifications)}\t"
                    f"{row.idt_code}\n"
                ),
            )
        text.configure(state="disabled")

    def _chemopt_tab_stops(self, rows) -> tuple[str, ...]:
        from tkinter import font as tkfont

        table_font = tkfont.Font(font=self.chemopt_output_text.cget("font"))
        padding = 28
        columns = [
            ("Variant", [str(row.row_number) for row in rows]),
            ("Motif ASO positions", [self._format_position_tuple(row.motif_aso_positions) for row in rows]),
            ("Motif core positions", [self._format_position_tuple(row.motif_gap_positions) for row in rows]),
            ("Sequence 5' to 3'", [row.clean_sequence for row in rows]),
            ("Chemistry pattern", [self._format_base_mod_pattern(row.base_modifications) for row in rows]),
        ]

        tab_positions: list[str] = []
        x = 0
        for header, values in columns:
            longest = max((header, *values), key=table_font.measure)
            x += table_font.measure(longest) + padding
            tab_positions.append(str(x))
        return tuple(tab_positions)

    def _format_position_tuple(self, positions: tuple[int, ...]) -> str:
        if not positions:
            return ""
        sorted_positions = sorted(positions)
        ranges: list[str] = []
        start = prev = sorted_positions[0]
        for position in sorted_positions[1:]:
            if position == prev + 1:
                prev = position
                continue
            ranges.append(str(start) if start == prev else f"{start}-{prev}")
            start = prev = position
        ranges.append(str(start) if start == prev else f"{start}-{prev}")
        return ", ".join(ranges)

    def _format_base_mod_pattern(self, base_modifications: tuple[str, ...]) -> str:
        if not base_modifications:
            return ""
        chunks: list[str] = []
        start = 1
        current = base_modifications[0]
        for idx, modification in enumerate(base_modifications[1:], start=2):
            if modification == current:
                continue
            label = str(start) if start == idx - 1 else f"{start}-{idx - 1}"
            chunks.append(f"{label} {current}")
            start = idx
            current = modification
        label = str(start) if start == len(base_modifications) else f"{start}-{len(base_modifications)}"
        chunks.append(f"{label} {current}")
        return " | ".join(chunks)

    def _refresh_chemopt_bubble_figure(self) -> None:
        if self.last_chemopt_rows:
            self._render_chemopt_bubble_figure(self.last_chemopt_rows)

    def _chemopt_bubble_base_letters_visible(self) -> bool:
        variable = getattr(self, "show_chemopt_bubble_base_letters", None)
        if variable is None:
            return True
        return bool(variable.get())

    def _chemopt_bubble_row_label(self, row) -> str:
        return f"Variant {row.row_number} (core {self._format_position_tuple(row.motif_gap_positions)})"

    def _chemopt_bubble_layout(self, rows) -> dict[str, int | float]:
        radius = 12
        pitch = radius * 2
        row_h = 30
        seq_len = len(rows[0].clean_sequence) if rows else 0
        labels = ["ASO sequence 5' - 3'", "Chemistry variants"] + [
            self._chemopt_bubble_row_label(row) for row in rows
        ]
        label_w = max(260, max([len(label) for label in labels] or [0]) * 7 + 24)
        left = 14
        top = 16
        grid_left = left + label_w
        position_y = top + 8
        sequence_y = top + 32
        direction_y = sequence_y + 36
        row_top = direction_y + 24
        grid_right = grid_left + seq_len * pitch
        legend_x = grid_right + 56
        legend_y = row_top + 84
        width = max(grid_right + 48, legend_x + 250, left + 840)
        height = max(row_top + len(rows) * row_h + 40, legend_y + 190)
        return {
            "left_label_x": left + label_w - 8,
            "grid_left": grid_left,
            "position_y": position_y,
            "sequence_y": sequence_y,
            "direction_y": direction_y,
            "row_top": row_top,
            "row_h": row_h,
            "radius": radius,
            "pitch": pitch,
            "diamond_r": 5,
            "legend_x": legend_x,
            "legend_y": legend_y,
            "width": width,
            "height": height,
        }

    def _chemopt_used_modifications(self, rows) -> tuple[list[str], list[str]]:
        base_set: set[str] = set()
        linkage_set: set[str] = set()
        methyl_used = False
        for row in rows:
            for idx, base in enumerate(row.clean_sequence):
                if idx < len(row.base_modifications):
                    base_set.add(self._bubble_ribose_mod(row.base_modifications[idx]))
                    methyl_used = methyl_used or self._bubble_has_5mec(base, row.base_modifications[idx])
            linkage_set.update(row.linkages)
        used_bases = [mod for mod in RIBOSE_MODIFICATION_OPTIONS if mod in base_set]
        if methyl_used:
            used_bases.append("5MeC")
        used_linkages = [linkage for linkage in LINKAGE_OPTIONS if linkage in linkage_set]
        return used_bases, used_linkages

    def _render_chemopt_bubble_figure(self, rows) -> None:
        canvas = self.chemopt_bubble_canvas
        canvas.delete("all")
        if not rows:
            canvas.configure(scrollregion=(0, 0, 0, 0))
            return

        layout = self._chemopt_bubble_layout(rows)
        styles = self._bubble_styles()
        used_base_mods, used_linkages = self._chemopt_used_modifications(rows)
        show_letters = self._chemopt_bubble_base_letters_visible()
        sequence = rows[0].clean_sequence

        self._draw_bubble_canvas_label(canvas, layout["left_label_x"], layout["sequence_y"], "ASO sequence 5' - 3'")
        for idx, base in enumerate(sequence):
            x = layout["grid_left"] + idx * layout["pitch"] + layout["radius"]
            canvas.create_text(
                x,
                layout["position_y"],
                text=str(idx + 1),
                anchor="s",
                font=("Menlo", 8),
                fill="#111827",
            )
            self._draw_bubble_canvas_circle(
                canvas,
                x,
                layout["sequence_y"],
                layout["radius"],
                base if show_letters else "",
                styles["rna_fill"],
                "#111827",
            )

        self._draw_bubble_canvas_label(
            canvas,
            layout["left_label_x"],
            layout["direction_y"],
            "Chemistry variants",
        )

        for row in rows:
            y = layout["row_top"] + (row.row_number - 1) * layout["row_h"]
            self._draw_bubble_canvas_label(
                canvas,
                layout["left_label_x"],
                y,
                self._chemopt_bubble_row_label(row),
                font_size=9,
            )
            centers: dict[int, tuple[float, float]] = {}
            for idx, base in enumerate(row.clean_sequence):
                if idx >= len(row.base_modifications):
                    continue
                x = layout["grid_left"] + idx * layout["pitch"] + layout["radius"]
                modification = row.base_modifications[idx]
                mod = self._bubble_ribose_mod(modification)
                self._draw_bubble_canvas_circle(
                    canvas,
                    x,
                    y,
                    layout["radius"],
                    base if show_letters else "",
                    styles["base_colors"].get(mod, styles["base_colors"]["DNA"]),
                    styles["base_text_colors"].get(mod, "#111827"),
                    self._bubble_has_5mec(base, modification),
                )
                centers[idx] = (x, y)
            for idx, (x, y) in centers.items():
                if idx >= len(row.linkages) or idx + 1 not in centers:
                    continue
                self._draw_linkage_canvas_diamond(
                    canvas,
                    x + layout["radius"],
                    y,
                    layout["diamond_r"],
                    styles["linkage_colors"].get(row.linkages[idx], "#ffffff"),
                    "#111827",
                )

        legend_bottom = self._draw_bubble_canvas_legend(
            canvas,
            layout["legend_x"],
            layout["legend_y"],
            used_base_mods,
            used_linkages,
            styles,
        )
        canvas.configure(scrollregion=(0, 0, layout["width"], max(layout["height"], legend_bottom + 24)))

    def save_chemopt_bubble_figure_image(self) -> None:
        if not self.last_chemopt_rows:
            self.messagebox.showerror(APP_NAME, "Generate a chemistry walk before saving a bubble figure.")
            return
        path = self.filedialog.asksaveasfilename(
            title="Save high-resolution chemistry optimisation bubble figure",
            defaultextension=".png",
            filetypes=[("PNG image", "*.png")],
            initialfile="chemistry_optimisation_bubble_figure_600dpi.png",
        )
        if not path:
            return
        try:
            self._save_chemopt_bubble_figure_png(self.last_chemopt_rows, Path(path))
        except ModuleNotFoundError as exc:
            if exc.name == "PIL":
                self.messagebox.showerror(
                    APP_NAME,
                    "PNG export requires Pillow. Install requirements again, then reopen the app.",
                )
            else:
                self.messagebox.showerror(APP_NAME, str(exc))
            return
        except Exception as exc:
            self.messagebox.showerror(APP_NAME, str(exc))
            return
        self.messagebox.showinfo(APP_NAME, f"Saved high-resolution PNG:\n{path}")

    def _save_chemopt_bubble_figure_png(self, rows, output_path: Path, show_letters: bool | None = None) -> Path:
        from PIL import Image, ImageDraw, ImageFont

        layout = self._chemopt_bubble_layout(rows)
        styles = self._bubble_styles()
        used_base_mods, used_linkages = self._chemopt_used_modifications(rows)
        if show_letters is None:
            show_letters = self._chemopt_bubble_base_letters_visible()
        scale = BUBBLE_EXPORT_SCALE
        width = int(layout["width"] * scale)
        height = int(layout["height"] * scale)
        image = Image.new("RGB", (width, height), "white")
        draw = ImageDraw.Draw(image)

        def s(value: float) -> int:
            return int(round(float(value) * scale))

        def font(size: int, bold: bool = False):
            scaled_size = size * scale
            candidates = [
                "/System/Library/Fonts/Supplemental/Arial Bold.ttf" if bold else "/System/Library/Fonts/Supplemental/Arial.ttf",
                "arialbd.ttf" if bold else "arial.ttf",
                "DejaVuSans-Bold.ttf" if bold else "DejaVuSans.ttf",
            ]
            for candidate in candidates:
                try:
                    return ImageFont.truetype(candidate, scaled_size)
                except Exception:
                    continue
            return ImageFont.load_default()

        label_font = font(13, True)
        row_label_font = font(12, True)
        base_font = font(12, True)
        pos_font = font(10, False)
        legend_title_font = font(13, True)
        legend_font = font(12, False)

        def text_right(x: float, y: float, value: str, fnt, fill: str = "#1f4e79") -> None:
            bbox = draw.textbbox((0, 0), value, font=fnt)
            draw.text((s(x) - (bbox[2] - bbox[0]), s(y) - (bbox[3] - bbox[1]) / 2), value, font=fnt, fill=fill)

        def text_center(x: float, y: float, value: str, fnt, fill: str = "#111827") -> None:
            bbox = draw.textbbox((0, 0), value, font=fnt)
            draw.text(
                (s(x) - (bbox[2] - bbox[0]) / 2, s(y) - (bbox[3] - bbox[1]) / 2 - scale),
                value,
                font=fnt,
                fill=fill,
            )

        def circle(
            x: float,
            y: float,
            radius: float,
            value: str,
            fill: str,
            text_fill: str,
            methyl_marker: bool = False,
        ) -> None:
            draw.ellipse(
                (s(x - radius), s(y - radius), s(x + radius), s(y + radius)),
                fill=fill,
                outline="#111827",
                width=s(2),
            )
            if value:
                text_center(x, y, value, base_font, text_fill)
            if methyl_marker:
                marker_radius = max(3, radius * 0.34)
                draw.ellipse(
                    (
                        s(x + radius * 0.32),
                        s(y - radius * 0.92),
                        s(x + radius * 0.32 + marker_radius),
                        s(y - radius * 0.92 + marker_radius),
                    ),
                    fill=styles["methyl_marker"],
                    outline="#111827",
                    width=max(1, s(1)),
                )

        def diamond(x: float, y: float, radius: float, fill: str) -> None:
            points = [(s(x), s(y - radius)), (s(x + radius), s(y)), (s(x), s(y + radius)), (s(x - radius), s(y))]
            draw.polygon(points, fill=fill)
            draw.line(points + [points[0]], fill="#111827", width=max(1, s(1.3)))

        sequence = rows[0].clean_sequence
        text_right(layout["left_label_x"], layout["sequence_y"], "ASO sequence 5' - 3'", label_font)
        for idx, base in enumerate(sequence):
            x = layout["grid_left"] + idx * layout["pitch"] + layout["radius"]
            text_center(x, layout["position_y"], str(idx + 1), pos_font)
            circle(x, layout["sequence_y"], layout["radius"], base if show_letters else "", styles["rna_fill"], "#111827")

        text_right(layout["left_label_x"], layout["direction_y"], "Chemistry variants", label_font)
        for row in rows:
            y = layout["row_top"] + (row.row_number - 1) * layout["row_h"]
            text_right(layout["left_label_x"], y, self._chemopt_bubble_row_label(row), row_label_font)
            centers: dict[int, tuple[float, float]] = {}
            for idx, base in enumerate(row.clean_sequence):
                if idx >= len(row.base_modifications):
                    continue
                x = layout["grid_left"] + idx * layout["pitch"] + layout["radius"]
                modification = row.base_modifications[idx]
                mod = self._bubble_ribose_mod(modification)
                circle(
                    x,
                    y,
                    layout["radius"],
                    base if show_letters else "",
                    styles["base_colors"].get(mod, styles["base_colors"]["DNA"]),
                    styles["base_text_colors"].get(mod, "#111827"),
                    self._bubble_has_5mec(base, modification),
                )
                centers[idx] = (x, y)
            for idx, (x, y) in centers.items():
                if idx >= len(row.linkages) or idx + 1 not in centers:
                    continue
                diamond(x + layout["radius"], y, layout["diamond_r"], styles["linkage_colors"].get(row.linkages[idx], "#ffffff"))

        self._draw_bubble_png_legend(
            draw,
            layout["legend_x"],
            layout["legend_y"],
            used_base_mods,
            used_linkages,
            styles,
            legend_title_font,
            legend_font,
            scale,
        )

        output_path.parent.mkdir(parents=True, exist_ok=True)
        image.save(output_path, dpi=(BUBBLE_EXPORT_DPI, BUBBLE_EXPORT_DPI))
        return output_path

    def _initial_converter_custom_pattern(self) -> tuple[list[str], list[str]]:
        if self.converter_custom_base_modifications:
            return list(self.converter_custom_base_modifications), list(self.converter_custom_linkages)
        try:
            base_mods, linkages = build_gapmer_pattern(
                int(self.vars["converter_gap_length"].get()),
                int(self.vars["converter_wing_length"].get()),
                self.vars["converter_wing_chemistry"].get(),
                self.vars["converter_backbone_modification"].get(),
            )
        except Exception:
            base_mods, linkages = build_gapmer_pattern(12, 3, "LNA", "PS")
        return list(base_mods), list(linkages)

    def _set_converter_custom_pattern(self, base_mods: list[str], linkages: list[str]) -> None:
        self.converter_custom_base_modifications = tuple(base_mods)
        self.converter_custom_linkages = tuple(linkages)
        self.vars["converter_gap_length"].set(str(len(base_mods)))
        self.vars["converter_wing_length"].set("0")
        self.vars["converter_wing_chemistry"].set("Custom")
        linkage_set = set(linkages)
        if linkage_set == {"PS"}:
            backbone = "PS"
        elif linkage_set <= {"PO"}:
            backbone = "PO"
        else:
            backbone = "mixed PS/PO"
        self.vars["converter_backbone_modification"].set(backbone)

    def _open_converter_custom_chemistry_editor(self) -> None:
        self._open_custom_chemistry_editor(
            chemistry_var_name="converter_aso_chemistry",
            initial_pattern_func=self._initial_converter_custom_pattern,
            apply_func=self._set_converter_custom_pattern,
            window_attr="converter_custom_editor_window",
            title="Custom IDT Converter Chemistry Pattern",
        )

    def _on_converter_chemistry_changed(self) -> None:
        self._apply_converter_chemistry_preset()
        if self.vars["converter_aso_chemistry"].get() == "Custom":
            self.root.after(80, self._open_converter_custom_chemistry_editor)

    def _apply_converter_chemistry_preset(self) -> None:
        label = self.vars["converter_aso_chemistry"].get()
        preset = CHEMISTRY_PRESETS.get(label)
        custom = preset is None
        if preset is not None:
            self.vars["converter_gap_length"].set(str(preset["gap_length"]))
            self.vars["converter_wing_length"].set(str(preset["wing_length"]))
            self.vars["converter_wing_chemistry"].set(str(preset["wing_chemistry"]))
            self.vars["converter_backbone_modification"].set(str(preset["backbone_modification"]))
            self.converter_custom_base_modifications = ()
            self.converter_custom_linkages = ()
        state = "disabled"
        combo_state = "disabled"
        self.converter_gap_entry.configure(state=state)
        self.converter_wing_entry.configure(state=state)
        self.converter_wing_combo.configure(state=combo_state)
        self.converter_backbone_combo.configure(state=combo_state)
        self.converter_custom_button.configure(state="normal" if custom else "disabled")

    def _collect_converter_chemistry_inputs(self) -> AsoInputs:
        return AsoInputs(
            aso_chemistry=self.vars["converter_aso_chemistry"].get(),
            gap_length=int(self.vars["converter_gap_length"].get()),
            wing_length=int(self.vars["converter_wing_length"].get()),
            wing_chemistry=self.vars["converter_wing_chemistry"].get(),
            backbone_modification=self.vars["converter_backbone_modification"].get(),
            custom_base_modifications=(
                self.converter_custom_base_modifications
                if self.vars["converter_aso_chemistry"].get() == "Custom"
                else ()
            ),
            custom_linkages=(
                self.converter_custom_linkages if self.vars["converter_aso_chemistry"].get() == "Custom" else ()
            ),
        )

    def convert_idt_notation(self) -> None:
        conversion = self._build_idt_conversion(show_errors=True)
        if conversion is None:
            return
        chemistry, rows = conversion
        self.last_idt_conversion_chemistry = chemistry
        self.last_idt_conversion_rows = rows
        self._render_idt_conversion_output(chemistry, rows)

    def _build_idt_conversion(self, *, show_errors: bool = True):
        try:
            chemistry = resolve_chemistry(self._collect_converter_chemistry_inputs())
        except Exception as exc:
            if show_errors:
                self.messagebox.showerror(APP_NAME, str(exc))
            return None

        rows = convert_sequences_to_idt(self.idt_sequence_text.get("1.0", "end"), chemistry)
        if not rows:
            if show_errors:
                self.messagebox.showerror(APP_NAME, "Enter at least one ASO sequence to convert.")
            return None
        return chemistry, rows

    def _render_idt_conversion_output(self, chemistry, rows) -> None:
        text = self.idt_output_text
        text.configure(state="normal")
        text.delete("1.0", "end")
        text.configure(tabs=self._idt_converter_tab_stops(rows))
        text.insert("end", f"Chemistry\t{chemistry_display_label(chemistry)}\n")
        text.insert("end", f"Expected length\t{chemistry.aso_length} bases\n\n")
        text.insert("end", "No.\tSequence 5' to 3'\tClean sequence\tIDT notation\n")
        for row in rows:
            text.insert(
                "end",
                f"{row.row_number}\t{row.input_sequence}\t{row.clean_sequence}\t{row.idt_code}\n",
            )
        text.configure(state="disabled")

    def _idt_converter_tab_stops(self, rows) -> tuple[str, ...]:
        from tkinter import font as tkfont

        table_font = tkfont.Font(font=self.idt_output_text.cget("font"))
        padding = 28
        columns = [
            ("No.", [str(row.row_number) for row in rows]),
            ("Sequence 5' to 3'", [row.input_sequence for row in rows]),
            ("Clean sequence", [row.clean_sequence for row in rows]),
        ]

        tab_positions: list[str] = []
        x = 0
        for header, values in columns:
            longest = max((header, *values), key=table_font.measure)
            x += table_font.measure(longest) + padding
            tab_positions.append(str(x))
        return tuple(tab_positions)

    def export_idt_converter_excel(self) -> None:
        conversion = self._build_idt_conversion(show_errors=True)
        if conversion is None:
            return
        chemistry, rows = conversion
        self.last_idt_conversion_chemistry = chemistry
        self.last_idt_conversion_rows = rows
        self._render_idt_conversion_output(chemistry, rows)

        path = self.filedialog.asksaveasfilename(
            title="Export IDT notation workbook",
            defaultextension=".xlsx",
            filetypes=[("Excel workbook", "*.xlsx")],
            initialfile="idt_notation_output.xlsx",
        )
        if not path:
            return
        try:
            from excel_export import export_idt_conversion_to_xlsx

            output = export_idt_conversion_to_xlsx(rows, chemistry, path)
        except Exception as exc:
            self.messagebox.showerror(APP_NAME, str(exc))
            return
        self.messagebox.showinfo(APP_NAME, f"Exported {output}")

    def _initial_penalty_custom_pattern(self) -> tuple[list[str], list[str]]:
        if self.penalty_custom_base_modifications:
            return list(self.penalty_custom_base_modifications), list(self.penalty_custom_linkages)
        try:
            base_mods, linkages = build_gapmer_pattern(
                int(self.vars["penalty_gap_length"].get()),
                int(self.vars["penalty_wing_length"].get()),
                self.vars["penalty_wing_chemistry"].get(),
                self.vars["penalty_backbone_modification"].get(),
            )
        except Exception:
            base_mods, linkages = build_gapmer_pattern(12, 3, "LNA", "PS")
        return list(base_mods), list(linkages)

    def _set_penalty_custom_pattern(self, base_mods: list[str], linkages: list[str]) -> None:
        self.penalty_custom_base_modifications = tuple(base_mods)
        self.penalty_custom_linkages = tuple(linkages)
        self.vars["penalty_gap_length"].set(str(len(base_mods)))
        self.vars["penalty_wing_length"].set("0")
        self.vars["penalty_wing_chemistry"].set("Custom")
        linkage_set = set(linkages)
        if linkage_set == {"PS"}:
            backbone = "PS"
        elif linkage_set <= {"PO"}:
            backbone = "PO"
        else:
            backbone = "mixed PS/PO"
        self.vars["penalty_backbone_modification"].set(backbone)

    def _open_penalty_custom_chemistry_editor(self) -> None:
        self._open_custom_chemistry_editor(
            chemistry_var_name="penalty_aso_chemistry",
            initial_pattern_func=self._initial_penalty_custom_pattern,
            apply_func=self._set_penalty_custom_pattern,
            window_attr="penalty_custom_editor_window",
            title="Custom Penalty ASO Chemistry Pattern",
        )

    def _on_penalty_chemistry_changed(self) -> None:
        self._apply_penalty_chemistry_preset()
        if self.vars["penalty_aso_chemistry"].get() == "Custom":
            self.root.after(80, self._open_penalty_custom_chemistry_editor)

    def _apply_penalty_chemistry_preset(self) -> None:
        if "penalty_aso_chemistry" not in self.vars:
            return
        label = self.vars["penalty_aso_chemistry"].get()
        preset = CHEMISTRY_PRESETS.get(label)
        custom = preset is None
        if preset is not None:
            self.vars["penalty_gap_length"].set(str(preset["gap_length"]))
            self.vars["penalty_wing_length"].set(str(preset["wing_length"]))
            self.vars["penalty_wing_chemistry"].set(str(preset["wing_chemistry"]))
            self.vars["penalty_backbone_modification"].set(str(preset["backbone_modification"]))
            self.penalty_custom_base_modifications = ()
            self.penalty_custom_linkages = ()
        state = "disabled"
        combo_state = "disabled"
        self.penalty_gap_entry.configure(state=state)
        self.penalty_wing_entry.configure(state=state)
        self.penalty_wing_combo.configure(state=combo_state)
        self.penalty_backbone_combo.configure(state=combo_state)
        self.penalty_custom_button.configure(state="normal" if custom else "disabled")

    def _collect_penalty_inputs(self) -> PenaltyAsoInputs:
        return PenaltyAsoInputs(
            target_gene=self.vars["penalty_target_gene"].get(),
            target_identifier=self.vars["penalty_target_identifier"].get(),
            chemistry_number=self.vars["penalty_chemistry_number"].get(),
            aso_chemistry=self.vars["penalty_aso_chemistry"].get(),
            gap_length=int(self.vars["penalty_gap_length"].get()),
            wing_length=int(self.vars["penalty_wing_length"].get()),
            wing_chemistry=self.vars["penalty_wing_chemistry"].get(),
            backbone_modification=self.vars["penalty_backbone_modification"].get(),
            custom_base_modifications=(
                self.penalty_custom_base_modifications
                if self.vars["penalty_aso_chemistry"].get() == "Custom"
                else ()
            ),
            custom_linkages=(
                self.penalty_custom_linkages if self.vars["penalty_aso_chemistry"].get() == "Custom" else ()
            ),
            parent_start=int(self.vars["penalty_parent_start"].get()),
            parent_count=int(self.vars["penalty_parent_count"].get()),
            microwalk_step_size=int(self.vars["penalty_step_size"].get()),
            penalty_position_mode=self.vars["penalty_position_mode"].get(),
            selected_penalty_positions=self.vars["penalty_selected_positions"].get(),
            penalty_base_mode=self.vars["penalty_base_mode"].get(),
            rna_sequence=self.penalty_sequence_text.get("1.0", "end").strip(),
        )

    def generate_penalty_design(self, show_errors: bool = True) -> None:
        try:
            result = generate_penalty_design(self._collect_penalty_inputs())
        except Exception as exc:
            if show_errors:
                self.messagebox.showerror(APP_NAME, str(exc))
            return
        self.last_penalty_result = result
        self._render_penalty_output(result)
        self._render_penalty_bubble_figure(result)

    def _render_penalty_output(self, result) -> None:
        text = self.penalty_output_text
        text.configure(state="normal")
        text.delete("1.0", "end")
        text.configure(tabs=self._penalty_tab_stops(result))
        text.insert("end", f"Chemistry\t{chemistry_display_label(result.chemistry)}\n")
        text.insert("end", f"Penalty ASOs generated\t{len(result.rows)}\n\n")
        text.insert(
            "end",
            (
                "No.\tParent ASO ID\tPenalty ASO ID\tParent start\tPenalty ASO pos\t"
                "RNA 3' to 5' pos\tRNA 5' to 3' pos\tRNA context\tCanonical\tPenalty\t"
                "Mismatch\tPriority\tScore\tSequence 5' to 3'\tIDT notation\tReason\n"
            ),
        )
        for row in result.rows:
            line_before_priority = (
                f"{row.row_number}\t"
                f"{row.parent_aso_id}\t"
                f"{row.penalty_aso_id}\t"
                f"{row.parent_start}\t"
                f"{row.penalty_aso_position}\t"
                f"{row.target_position_3to5}\t"
                f"{row.target_position_5to3}\t"
                f"{row.local_rna_context}\t"
                f"{row.canonical_aso_base}\t"
                f"{row.penalty_aso_base}\t"
                f"{row.mismatch_pair}\t"
            )
            text.insert("end", line_before_priority)
            tag = {
                "Recommended": "priority_recommended",
                "Alternative": "priority_alternative",
                "Lower priority": "priority_lower",
            }.get(row.priority, "")
            text.insert("end", row.priority, tag)
            text.insert(
                "end",
                f"\t{row.score}\t{row.clean_sequence}\t{row.idt_code}\t{row.reason}\n",
            )
        text.configure(state="disabled")

    def _penalty_tab_stops(self, result) -> tuple[str, ...]:
        from tkinter import font as tkfont

        table_font = tkfont.Font(font=self.penalty_output_text.cget("font"))
        padding = 28
        columns = [
            ("No.", [str(row.row_number) for row in result.rows]),
            ("Parent ASO ID", [row.parent_aso_id for row in result.rows]),
            ("Penalty ASO ID", [row.penalty_aso_id for row in result.rows]),
            ("Parent start", [str(row.parent_start) for row in result.rows]),
            ("Penalty ASO pos", [str(row.penalty_aso_position) for row in result.rows]),
            ("RNA 3' to 5' pos", [str(row.target_position_3to5) for row in result.rows]),
            ("RNA 5' to 3' pos", [str(row.target_position_5to3) for row in result.rows]),
            ("RNA context", [row.local_rna_context for row in result.rows]),
            ("Canonical", [row.canonical_aso_base for row in result.rows]),
            ("Penalty", [row.penalty_aso_base for row in result.rows]),
            ("Mismatch", [row.mismatch_pair for row in result.rows]),
            ("Priority", [row.priority for row in result.rows]),
            ("Score", [str(row.score) for row in result.rows]),
            ("Sequence 5' to 3'", [row.clean_sequence for row in result.rows]),
            ("IDT notation", [row.idt_code for row in result.rows]),
        ]

        tab_positions: list[str] = []
        x = 0
        for header, values in columns:
            longest = max((header, *values), key=table_font.measure)
            x += table_font.measure(longest) + padding
            tab_positions.append(str(x))
        return tuple(tab_positions)

    def export_penalty_design_excel(self) -> None:
        self.generate_penalty_design(show_errors=True)
        if self.last_penalty_result is None:
            return
        path = self.filedialog.asksaveasfilename(
            title="Export penalty ASO workbook",
            defaultextension=".xlsx",
            filetypes=[("Excel workbook", "*.xlsx")],
            initialfile="penalty_aso_design_output.xlsx",
        )
        if not path:
            return
        try:
            from excel_export import export_penalty_design_to_xlsx

            output = export_penalty_design_to_xlsx(self.last_penalty_result, path)
        except Exception as exc:
            self.messagebox.showerror(APP_NAME, str(exc))
            return
        self.messagebox.showinfo(APP_NAME, f"Exported {output}")

    def _collect_inputs(self) -> AsoInputs:
        return AsoInputs(
            target_gene=self.vars["target_gene"].get(),
            snp_identifier=self.vars["snp_identifier"].get(),
            chemistry_number=self.vars["chemistry_number"].get(),
            aso_chemistry=self.vars["aso_chemistry"].get(),
            gap_length=int(self.vars["gap_length"].get()),
            wing_length=int(self.vars["wing_length"].get()),
            wing_chemistry=self.vars["wing_chemistry"].get(),
            backbone_modification=self.vars["backbone_modification"].get(),
            custom_base_modifications=(
                self.custom_base_modifications if self.vars["aso_chemistry"].get() == "Custom" else ()
            ),
            custom_linkages=(
                self.custom_linkages if self.vars["aso_chemistry"].get() == "Custom" else ()
            ),
            mutation_type=self.vars["mutation_type"].get(),
            mutation_length=int(self.vars["mutation_length"].get()),
            mutation_start=int(self.vars["mutation_start"].get()),
            microwalk_step_size=self.vars["microwalk_step_size"].get(),
            rna_sequence=self.sequence_text.get("1.0", "end").strip(),
        )

    def calculate(self, show_errors: bool = True) -> None:
        try:
            result = generate_design(self._collect_inputs())
        except Exception as exc:
            if show_errors:
                self.messagebox.showerror(APP_NAME, str(exc))
            return
        self.last_result = result
        self._show_result(result)

    def _show_result(self, result) -> None:
        self._variant_render_token += 1
        render_token = self._variant_render_token
        self._variant_bubble_dirty = True

        steps = (
            lambda: self._render_aso_table(result),
            lambda: self._render_alignment(result),
            self._render_variant_bubble_if_visible,
        )

        def run_step(index: int = 0) -> None:
            if render_token != self._variant_render_token:
                return
            if index >= len(steps):
                return
            steps[index]()
            if index + 1 < len(steps):
                self.root.after(1, lambda: run_step(index + 1))

        self.root.after_idle(run_step)

    def _render_aso_table(self, result) -> None:
        text = self.table_text
        text.configure(state="normal")
        text.delete("1.0", "end")
        text.configure(tabs=self._table_tab_stops(result))

        headers = (
            "ASO ID\t"
            "Sequence 5' to 3' with IDT Codes\t"
            "Highlighted Sequence 5' to 3'\t"
            "Indel/Substitution Start Position\t"
            "Chemistry\n"
        )
        text.insert("end", headers, "table_header")

        for row in result.rows:
            text.insert("end", f"{row.aso_id}\t{row.idt_code}\t")
            self._insert_display_sequence(text, row.display_sequence, row.display_spans)
            text.insert("end", f"\t{row.starting_position}\t{row.chemistry}\n")

        text.configure(state="disabled")

    def _table_tab_stops(self, result) -> tuple[str, ...]:
        from tkinter import font as tkfont

        table_font = tkfont.Font(font=self.table_text.cget("font"))
        padding = 28
        columns = [
            ("ASO ID", [row.aso_id for row in result.rows]),
            ("Sequence 5' to 3' with IDT Codes", [row.idt_code for row in result.rows]),
            ("Highlighted Sequence 5' to 3'", [row.display_sequence for row in result.rows]),
            (
                "Indel/Substitution Start Position",
                [str(row.starting_position) for row in result.rows],
            ),
        ]

        tab_positions: list[str] = []
        x = 0
        for header, values in columns:
            longest = max((header, *values), key=table_font.measure)
            x += table_font.measure(longest) + padding
            tab_positions.append(str(x))
        return tuple(tab_positions)

    def _insert_display_sequence(self, text, sequence: str, spans) -> None:
        cursor = 0
        for span in spans:
            if span.start > cursor:
                text.insert("end", sequence[cursor : span.start])
            if span.kind == "mutation":
                text.insert("end", sequence[span.start : span.end], "display_mutation")
            else:
                text.insert("end", sequence[span.start : span.end])
            cursor = span.end
        if cursor < len(sequence):
            text.insert("end", sequence[cursor:])

    def _render_alignment(self, result) -> None:
        canvas = self.alignment_canvas
        canvas.delete("all")

        cell_w = max(34, 12 * len(str(max(result.header_positions, default=1))) + 14)
        cell_h = 24
        max_label_chars = max(
            [len("RNA 3' - 5'"), len("ASO sequences 5' - 3'"), len("pos")]
            + [len(row.aso_id) for row in result.rows]
        )
        label_w = max(180, max_label_chars * 8 + 24)
        left = 12
        top = 12
        highlighted_header_indexes = mutation_header_indexes(result)

        self._draw_alignment_label(canvas, left, top, label_w, cell_h, "pos")
        self._draw_alignment_label(canvas, left, top + cell_h, label_w, cell_h, "RNA 3' - 5'")
        for idx, position in enumerate(result.header_positions):
            x = left + label_w + idx * cell_w
            self._draw_alignment_cell(canvas, x, top, cell_w, cell_h, str(position), "#e7e6e6", "#000000")
            if idx in highlighted_header_indexes:
                fill = "#fce4d6"
                text_color = "#c00000"
            else:
                fill = "#1f4e79"
                text_color = "#ffffff"
            self._draw_alignment_cell(canvas, x, top + cell_h, cell_w, cell_h, result.header_bases[idx], fill, text_color)

        aso_direction_y = top + 2 * cell_h + 10
        canvas.create_text(
            left + label_w - 6,
            aso_direction_y + cell_h / 2,
            text="ASO sequences 5' - 3'",
            anchor="e",
            font=("Menlo", 11, "bold"),
            fill="#1f4e79",
        )

        row_top = aso_direction_y + cell_h
        for row in result.rows:
            y = row_top + (row.row_number - 1) * cell_h
            self._draw_alignment_label(canvas, left, y, label_w, cell_h, row.aso_id)
            for idx, grid_value in enumerate(row.grid_cells):
                x = left + label_w + idx * cell_w
                used = grid_value != "##"
                fill = "#d9ead3" if used else "#e7e6e6"
                value = grid_value if used else ""
                self._draw_alignment_cell(canvas, x, y, cell_w, cell_h, value, fill, "#000000")

        bottom = row_top + len(result.rows) * cell_h
        for boundary in self._alignment_indel_boundaries(result):
            x = left + label_w + boundary * cell_w
            canvas.create_line(x, top, x, bottom, fill="#c00000", width=3)

        canvas.configure(scrollregion=canvas.bbox("all"))

    def _draw_alignment_label(self, canvas, x: int, y: int, width: int, height: int, label: str) -> None:
        canvas.create_rectangle(x, y, x + width, y + height, fill="#ffffff", outline="#d9d9d9")
        canvas.create_text(x + width - 6, y + height / 2, text=label, anchor="e", font=("Menlo", 11, "bold"), fill="#1f4e79")

    def _draw_alignment_cell(
        self,
        canvas,
        x: int,
        y: int,
        width: int,
        height: int,
        value: str,
        fill: str,
        text_color: str,
    ) -> None:
        canvas.create_rectangle(x, y, x + width, y + height, fill=fill, outline="#bfbfbf")
        canvas.create_text(x + width / 2, y + height / 2, text=value, font=("Menlo", 11, "bold"), fill=text_color)

    def _alignment_indel_boundaries(self, result) -> list[int]:
        mutation_type = normalise_mutation_type(result.inputs.mutation_type)
        if mutation_type == "DELETION":
            boundary = result.mutation_start_reversed - result.crop_start
            if 0 <= boundary <= len(result.header_positions):
                return [boundary]
            return []

        if mutation_type in {"INSERTION", "SUBSTITUTION"}:
            highlighted = mutation_header_indexes(result)
            if highlighted:
                return [min(highlighted), max(highlighted) + 1]
        return []

    def _render_bubble_figure(self, result) -> None:
        canvas = self.bubble_canvas
        canvas.delete("all")

        layout = self._bubble_layout(result)
        styles = self._bubble_styles()
        used_base_mods, used_linkages = self._bubble_used_modifications(result)
        show_letters = self._bubble_base_letters_visible()

        self._draw_bubble_canvas_label(canvas, layout["left_label_x"], layout["rna_y"], "RNA 3' - 5'")
        for idx, position in enumerate(result.header_positions):
            x = layout["grid_left"] + idx * layout["pitch"] + layout["radius"]
            canvas.create_text(
                x,
                layout["position_y"],
                text=str(position),
                anchor="s",
                font=("Menlo", 8),
                fill="#111827",
            )
            self._draw_bubble_canvas_circle(
                canvas,
                x,
                layout["rna_y"],
                layout["radius"],
                result.header_bases[idx] if show_letters else "",
                styles["rna_fill"],
                "#111827",
            )

        self._draw_bubble_canvas_label(
            canvas,
            layout["left_label_x"],
            layout["direction_y"],
            "ASO sequences 5' - 3'",
        )

        for row in result.rows:
            y = layout["row_top"] + (row.row_number - 1) * layout["row_h"]
            self._draw_bubble_canvas_label(canvas, layout["left_label_x"], y, row.aso_id, font_size=9)
            drawn_centers: dict[int, tuple[float, float]] = {}

            for idx, grid_value in enumerate(row.grid_cells):
                if grid_value == "##":
                    continue
                local_idx = idx - row.starting_position
                if not 0 <= local_idx < len(result.chemistry.base_modifications):
                    continue

                x = layout["grid_left"] + idx * layout["pitch"] + layout["radius"]
                modification = result.chemistry.base_modifications[local_idx]
                base_mod = self._bubble_ribose_mod(modification)
                fill = styles["base_colors"].get(base_mod, styles["base_colors"]["DNA"])
                text_color = styles["base_text_colors"].get(base_mod, "#111827")
                self._draw_bubble_canvas_circle(
                    canvas,
                    x,
                    y,
                    layout["radius"],
                    grid_value if show_letters else "",
                    fill,
                    text_color,
                    self._bubble_has_5mec(grid_value, modification),
                )
                drawn_centers[local_idx] = (x, y)

            for local_idx, (x, y) in drawn_centers.items():
                if local_idx >= len(result.chemistry.linkages) or local_idx + 1 not in drawn_centers:
                    continue
                linkage = result.chemistry.linkages[local_idx]
                self._draw_linkage_canvas_diamond(
                    canvas,
                    x + layout["radius"],
                    y,
                    layout["diamond_r"],
                    styles["linkage_colors"].get(linkage, "#ffffff"),
                    "#111827",
                )

        legend_bottom = self._draw_bubble_canvas_legend(
            canvas,
            layout["legend_x"],
            layout["legend_y"],
            used_base_mods,
            used_linkages,
            styles,
        )
        canvas.configure(scrollregion=(0, 0, layout["width"], max(layout["height"], legend_bottom + 24)))

    def _bubble_styles(self) -> dict:
        base_colors = {
            "DNA": "#ffffff",
            "LNA": "#2f80ed",
            "MOE": "#d92d20",
            "2'OMe": "#2e7d32",
            "2'F": "#7e57c2",
        }
        base_text_colors = {
            "DNA": "#111827",
            "LNA": "#ffffff",
            "MOE": "#ffffff",
            "2'OMe": "#ffffff",
            "2'F": "#ffffff",
        }
        linkage_colors = {
            "PS": "#111827",
            "PO": "#ffffff",
        }
        return {
            "base_colors": base_colors,
            "base_text_colors": base_text_colors,
            "linkage_colors": linkage_colors,
            "rna_fill": "#d9d9d9",
            "methyl_marker": "#f6c343",
        }

    def _bubble_layout(self, result) -> dict[str, int | float]:
        radius = 12
        pitch = radius * 2
        row_h = 30
        label_w = max(
            235,
            max([len("RNA 3' - 5'"), len("ASO sequences 5' - 3'")] + [len(row.aso_id) for row in result.rows])
            * 7
            + 24,
        )
        left = 14
        top = 16
        grid_left = left + label_w
        position_y = top + 8
        rna_y = top + 32
        direction_y = rna_y + 36
        row_top = direction_y + 24
        grid_right = grid_left + len(result.header_positions) * pitch
        legend_x = grid_right + 56
        legend_y = row_top + 84
        width = max(grid_right + 48, legend_x + 250, left + 840)
        height = max(row_top + len(result.rows) * row_h + 40, legend_y + 190)
        return {
            "left": left,
            "top": top,
            "label_w": label_w,
            "left_label_x": left + label_w - 8,
            "grid_left": grid_left,
            "position_y": position_y,
            "rna_y": rna_y,
            "direction_y": direction_y,
            "row_top": row_top,
            "row_h": row_h,
            "radius": radius,
            "pitch": pitch,
            "diamond_r": 5,
            "legend_x": legend_x,
            "legend_y": legend_y,
            "width": width,
            "height": height,
        }

    def _effective_bubble_base_mod(self, base: str, modification: str) -> str:
        ribose_mod, base_mod = split_base_modification(modification)
        if base_mod == "5MeC" and base.upper() != "C":
            return ribose_mod or "DNA"
        if base_mod == "5MeC":
            return combine_base_modification(ribose_mod or "DNA", "5MeC")
        return ribose_mod or modification

    def _bubble_ribose_mod(self, modification: str) -> str:
        ribose_mod, _base_mod = split_base_modification(modification)
        return ribose_mod or "DNA"

    def _bubble_has_5mec(self, base: str, modification: str) -> bool:
        _ribose_mod, base_mod = split_base_modification(modification)
        return base_mod == "5MeC" and base.upper() == "C"

    def _bubble_base_letters_visible(self) -> bool:
        variable = getattr(self, "show_bubble_base_letters", None)
        if variable is None:
            return True
        return bool(variable.get())

    def _bubble_used_modifications(self, result) -> tuple[list[str], list[str]]:
        base_set: set[str] = set()
        linkage_set: set[str] = set()
        methyl_used = False
        for row in result.rows:
            for idx, grid_value in enumerate(row.grid_cells):
                if grid_value == "##":
                    continue
                local_idx = idx - row.starting_position
                if not 0 <= local_idx < len(result.chemistry.base_modifications):
                    continue
                modification = result.chemistry.base_modifications[local_idx]
                base_set.add(self._bubble_ribose_mod(modification))
                methyl_used = methyl_used or self._bubble_has_5mec(grid_value, modification)
                if (
                    local_idx < len(result.chemistry.linkages)
                    and idx + 1 < len(row.grid_cells)
                    and row.grid_cells[idx + 1] != "##"
                ):
                    linkage_set.add(result.chemistry.linkages[local_idx])
        used_bases = [mod for mod in RIBOSE_MODIFICATION_OPTIONS if mod in base_set]
        if methyl_used:
            used_bases.append("5MeC")
        used_linkages = [linkage for linkage in LINKAGE_OPTIONS if linkage in linkage_set]
        return used_bases, used_linkages

    def _refresh_penalty_bubble_figure(self) -> None:
        if self.last_penalty_result is not None:
            self._render_penalty_bubble_figure(self.last_penalty_result)

    def _penalty_bubble_base_letters_visible(self) -> bool:
        variable = getattr(self, "show_penalty_bubble_base_letters", None)
        if variable is None:
            return True
        return bool(variable.get())

    def _penalty_bubble_layout(self, result) -> dict[str, int | float]:
        radius = 12
        pitch = radius * 2
        row_h = 30
        labels = ["RNA 3' - 5'", "Penalty ASOs 5' - 3'"] + [row.penalty_aso_id for row in result.rows]
        label_w = max(300, max([len(label) for label in labels] or [0]) * 7 + 24)
        left = 14
        top = 16
        grid_left = left + label_w
        position_y = top + 8
        rna_y = top + 32
        direction_y = rna_y + 36
        row_top = direction_y + 24
        grid_right = grid_left + len(result.header_positions) * pitch
        legend_x = grid_right + 56
        legend_y = row_top + 84
        width = max(grid_right + 48, legend_x + 280, left + 900)
        height = max(row_top + len(result.rows) * row_h + 40, legend_y + 240)
        return {
            "left_label_x": left + label_w - 8,
            "grid_left": grid_left,
            "position_y": position_y,
            "rna_y": rna_y,
            "direction_y": direction_y,
            "row_top": row_top,
            "row_h": row_h,
            "radius": radius,
            "pitch": pitch,
            "diamond_r": 5,
            "legend_x": legend_x,
            "legend_y": legend_y,
            "width": width,
            "height": height,
        }

    def _penalty_used_modifications(self, result) -> tuple[list[str], list[str]]:
        base_set: set[str] = set()
        linkage_set: set[str] = set()
        methyl_used = False
        for row in result.rows:
            for idx, grid_value in enumerate(row.grid_cells):
                if grid_value == "##":
                    continue
                local_idx = result.header_positions[idx] - row.parent_start
                if not 0 <= local_idx < len(result.chemistry.base_modifications):
                    continue
                modification = result.chemistry.base_modifications[local_idx]
                base_set.add(self._bubble_ribose_mod(modification))
                methyl_used = methyl_used or self._bubble_has_5mec(grid_value, modification)
                if (
                    local_idx < len(result.chemistry.linkages)
                    and idx + 1 < len(row.grid_cells)
                    and row.grid_cells[idx + 1] != "##"
                ):
                    linkage_set.add(result.chemistry.linkages[local_idx])
        used_bases = [mod for mod in RIBOSE_MODIFICATION_OPTIONS if mod in base_set]
        if methyl_used:
            used_bases.append("5MeC")
        used_linkages = [linkage for linkage in LINKAGE_OPTIONS if linkage in linkage_set]
        return used_bases, used_linkages

    def _draw_penalty_marker_legend(self, canvas, x: float, y: float) -> float:
        canvas.create_oval(x, y - 9, x + 18, y + 9, fill="#f6c343", outline="#c00000", width=2.4)
        canvas.create_text(
            x + 28,
            y,
            text="Intentional penalty base",
            anchor="w",
            font=("Menlo", 10),
            fill="#111827",
        )
        return y + 32

    def _render_penalty_bubble_figure(self, result) -> None:
        canvas = self.penalty_bubble_canvas
        canvas.delete("all")
        if not result.rows:
            canvas.configure(scrollregion=(0, 0, 0, 0))
            return

        layout = self._penalty_bubble_layout(result)
        styles = self._bubble_styles()
        used_base_mods, used_linkages = self._penalty_used_modifications(result)
        show_letters = self._penalty_bubble_base_letters_visible()

        self._draw_bubble_canvas_label(canvas, layout["left_label_x"], layout["rna_y"], "RNA 3' - 5'")
        for idx, position in enumerate(result.header_positions):
            x = layout["grid_left"] + idx * layout["pitch"] + layout["radius"]
            canvas.create_text(
                x,
                layout["position_y"],
                text=str(position),
                anchor="s",
                font=("Menlo", 8),
                fill="#111827",
            )
            self._draw_bubble_canvas_circle(
                canvas,
                x,
                layout["rna_y"],
                layout["radius"],
                result.header_bases[idx] if show_letters else "",
                styles["rna_fill"],
                "#111827",
            )

        self._draw_bubble_canvas_label(
            canvas,
            layout["left_label_x"],
            layout["direction_y"],
            "Penalty ASOs 5' - 3'",
        )

        for row in result.rows:
            y = layout["row_top"] + (row.row_number - 1) * layout["row_h"]
            self._draw_bubble_canvas_label(canvas, layout["left_label_x"], y, row.penalty_aso_id, font_size=9)
            drawn_centers: dict[int, tuple[float, float]] = {}

            for idx, grid_value in enumerate(row.grid_cells):
                if grid_value == "##":
                    continue
                local_idx = result.header_positions[idx] - row.parent_start
                if not 0 <= local_idx < len(result.chemistry.base_modifications):
                    continue

                x = layout["grid_left"] + idx * layout["pitch"] + layout["radius"]
                modification = result.chemistry.base_modifications[local_idx]
                base_mod = self._bubble_ribose_mod(modification)
                fill = styles["base_colors"].get(base_mod, styles["base_colors"]["DNA"])
                text_color = styles["base_text_colors"].get(base_mod, "#111827")
                if idx == row.penalty_grid_index:
                    fill = "#f6c343"
                    text_color = "#111827"
                self._draw_bubble_canvas_circle(
                    canvas,
                    x,
                    y,
                    layout["radius"],
                    grid_value if show_letters else "",
                    fill,
                    text_color,
                    self._bubble_has_5mec(grid_value, modification),
                )
                if idx == row.penalty_grid_index:
                    canvas.create_oval(
                        x - layout["radius"],
                        y - layout["radius"],
                        x + layout["radius"],
                        y + layout["radius"],
                        outline="#c00000",
                        width=2.6,
                    )
                drawn_centers[local_idx] = (x, y)

            for local_idx, (x, y) in drawn_centers.items():
                if local_idx >= len(result.chemistry.linkages) or local_idx + 1 not in drawn_centers:
                    continue
                linkage = result.chemistry.linkages[local_idx]
                self._draw_linkage_canvas_diamond(
                    canvas,
                    x + layout["radius"],
                    y,
                    layout["diamond_r"],
                    styles["linkage_colors"].get(linkage, "#ffffff"),
                    "#111827",
                )

        legend_bottom = self._draw_bubble_canvas_legend(
            canvas,
            layout["legend_x"],
            layout["legend_y"],
            used_base_mods,
            used_linkages,
            styles,
        )
        legend_bottom = self._draw_penalty_marker_legend(canvas, layout["legend_x"], legend_bottom)
        canvas.configure(scrollregion=(0, 0, layout["width"], max(layout["height"], legend_bottom + 24)))

    def _draw_bubble_canvas_label(
        self,
        canvas,
        x: float,
        y: float,
        label: str,
        *,
        font_size: int = 10,
    ) -> None:
        canvas.create_text(
            x,
            y,
            text=label,
            anchor="e",
            font=("Menlo", font_size, "bold"),
            fill="#1f4e79",
        )

    def _draw_bubble_canvas_circle(
        self,
        canvas,
        x: float,
        y: float,
        radius: float,
        text: str,
        fill: str,
        text_color: str,
        methyl_marker: bool = False,
    ) -> None:
        canvas.create_oval(
            x - radius,
            y - radius,
            x + radius,
            y + radius,
            fill=fill,
            outline="#111827",
            width=1.2,
        )
        if text:
            canvas.create_text(x, y, text=text, font=("Menlo", 9, "bold"), fill=text_color)
        if methyl_marker:
            marker_radius = max(3, radius * 0.34)
            canvas.create_oval(
                x + radius * 0.32,
                y - radius * 0.92,
                x + radius * 0.32 + marker_radius,
                y - radius * 0.92 + marker_radius,
                fill="#f6c343",
                outline="#111827",
                width=1,
            )

    def _draw_linkage_canvas_diamond(
        self,
        canvas,
        x: float,
        y: float,
        radius: float,
        fill: str,
        outline: str,
    ) -> None:
        canvas.create_polygon(
            x,
            y - radius,
            x + radius,
            y,
            x,
            y + radius,
            x - radius,
            y,
            fill=fill,
            outline=outline,
            width=1.1,
        )

    def _draw_bubble_canvas_legend(
        self,
        canvas,
        x: float,
        y: float,
        used_bases: list[str],
        used_linkages: list[str],
        styles: dict,
    ) -> float:
        canvas.create_text(x, y, text="Legend", anchor="w", font=("Menlo", 11, "bold"), fill="#111827")
        cursor_y = y + 28
        line_h = 28

        for mod in used_bases:
            if mod == "5MeC":
                label = "5MeC base modification"
                canvas.create_oval(
                    x + 5,
                    cursor_y - 5,
                    x + 15,
                    cursor_y + 5,
                    fill=styles["methyl_marker"],
                    outline="#111827",
                    width=1.2,
                )
            else:
                label = mod
                canvas.create_oval(
                    x,
                    cursor_y - 9,
                    x + 18,
                    cursor_y + 9,
                    fill=styles["base_colors"].get(mod, styles["base_colors"]["DNA"]),
                    outline="#111827",
                    width=1.2,
                )
            canvas.create_text(
                x + 28,
                cursor_y,
                text=label,
                anchor="w",
                font=("Menlo", 10),
                fill="#111827",
            )
            cursor_y += line_h

        for linkage in used_linkages:
            label = f"{linkage} linkage"
            tri_x = x + 9
            tri_y = cursor_y
            canvas.create_polygon(
                tri_x,
                tri_y - 8,
                tri_x + 8,
                tri_y,
                tri_x,
                tri_y + 8,
                tri_x - 8,
                tri_y,
                fill=styles["linkage_colors"].get(linkage, "#ffffff"),
                outline="#111827",
                width=1.2,
            )
            canvas.create_text(
                x + 28,
                cursor_y,
                text=label,
                anchor="w",
                font=("Menlo", 10),
                fill="#111827",
            )
            cursor_y += line_h

        return cursor_y + 8

    def save_bubble_figure_image(self) -> None:
        if self.last_result is None:
            self.messagebox.showerror(APP_NAME, "Calculate a microwalk before saving a bubble figure.")
            return
        path = self.filedialog.asksaveasfilename(
            title="Save high-resolution bubble figure",
            defaultextension=".png",
            filetypes=[("PNG image", "*.png")],
            initialfile="aso_bubble_figure_600dpi.png",
        )
        if not path:
            return
        try:
            self._save_bubble_figure_png(self.last_result, Path(path))
        except ModuleNotFoundError as exc:
            if exc.name == "PIL":
                self.messagebox.showerror(
                    APP_NAME,
                    "PNG export requires Pillow. Install requirements again, then reopen the app.",
                )
            else:
                self.messagebox.showerror(APP_NAME, str(exc))
            return
        except Exception as exc:
            self.messagebox.showerror(APP_NAME, str(exc))
            return
        self.messagebox.showinfo(APP_NAME, f"Saved high-resolution PNG:\n{path}")

    def _save_bubble_figure_png(self, result, output_path: Path, show_letters: bool | None = None) -> Path:
        from PIL import Image, ImageDraw, ImageFont

        layout = self._bubble_layout(result)
        styles = self._bubble_styles()
        used_base_mods, used_linkages = self._bubble_used_modifications(result)
        if show_letters is None:
            show_letters = self._bubble_base_letters_visible()
        scale = BUBBLE_EXPORT_SCALE
        width = int(layout["width"] * scale)
        height = int(layout["height"] * scale)
        image = Image.new("RGB", (width, height), "white")
        draw = ImageDraw.Draw(image)

        def s(value: float) -> int:
            return int(round(float(value) * scale))

        def font(size: int, bold: bool = False):
            scaled_size = size * scale
            candidates = [
                "/System/Library/Fonts/Supplemental/Arial Bold.ttf" if bold else "/System/Library/Fonts/Supplemental/Arial.ttf",
                "arialbd.ttf" if bold else "arial.ttf",
                "DejaVuSans-Bold.ttf" if bold else "DejaVuSans.ttf",
            ]
            for candidate in candidates:
                try:
                    return ImageFont.truetype(candidate, scaled_size)
                except Exception:
                    continue
            return ImageFont.load_default()

        label_font = font(13, True)
        row_label_font = font(12, True)
        base_font = font(12, True)
        pos_font = font(10, False)
        legend_title_font = font(13, True)
        legend_font = font(12, False)

        def text_right(x: float, y: float, value: str, fnt, fill: str = "#1f4e79") -> None:
            bbox = draw.textbbox((0, 0), value, font=fnt)
            draw.text((s(x) - (bbox[2] - bbox[0]), s(y) - (bbox[3] - bbox[1]) / 2), value, font=fnt, fill=fill)

        def text_center(x: float, y: float, value: str, fnt, fill: str = "#111827") -> None:
            bbox = draw.textbbox((0, 0), value, font=fnt)
            draw.text(
                (s(x) - (bbox[2] - bbox[0]) / 2, s(y) - (bbox[3] - bbox[1]) / 2 - scale),
                value,
                font=fnt,
                fill=fill,
            )

        def circle(
            x: float,
            y: float,
            radius: float,
            value: str,
            fill: str,
            text_fill: str,
            methyl_marker: bool = False,
        ) -> None:
            draw.ellipse(
                (s(x - radius), s(y - radius), s(x + radius), s(y + radius)),
                fill=fill,
                outline="#111827",
                width=s(2),
            )
            if value:
                text_center(x, y, value, base_font, text_fill)
            if methyl_marker:
                marker_radius = max(3, radius * 0.34)
                draw.ellipse(
                    (
                        s(x + radius * 0.32),
                        s(y - radius * 0.92),
                        s(x + radius * 0.32 + marker_radius),
                        s(y - radius * 0.92 + marker_radius),
                    ),
                    fill=styles["methyl_marker"],
                    outline="#111827",
                    width=max(1, s(1)),
                )

        def diamond(x: float, y: float, radius: float, fill: str) -> None:
            points = [(s(x), s(y - radius)), (s(x + radius), s(y)), (s(x), s(y + radius)), (s(x - radius), s(y))]
            draw.polygon(points, fill=fill)
            draw.line(points + [points[0]], fill="#111827", width=max(1, s(1.3)))

        text_right(layout["left_label_x"], layout["rna_y"], "RNA 3' - 5'", label_font)
        for idx, position in enumerate(result.header_positions):
            x = layout["grid_left"] + idx * layout["pitch"] + layout["radius"]
            text_center(x, layout["position_y"], str(position), pos_font)
            circle(
                x,
                layout["rna_y"],
                layout["radius"],
                result.header_bases[idx] if show_letters else "",
                styles["rna_fill"],
                "#111827",
            )

        text_right(layout["left_label_x"], layout["direction_y"], "ASO sequences 5' - 3'", label_font)
        for row in result.rows:
            y = layout["row_top"] + (row.row_number - 1) * layout["row_h"]
            text_right(layout["left_label_x"], y, row.aso_id, row_label_font)
            centers: dict[int, tuple[float, float]] = {}
            for idx, grid_value in enumerate(row.grid_cells):
                if grid_value == "##":
                    continue
                local_idx = idx - row.starting_position
                if not 0 <= local_idx < len(result.chemistry.base_modifications):
                    continue
                x = layout["grid_left"] + idx * layout["pitch"] + layout["radius"]
                modification = result.chemistry.base_modifications[local_idx]
                mod = self._bubble_ribose_mod(modification)
                circle(
                    x,
                    y,
                    layout["radius"],
                    grid_value if show_letters else "",
                    styles["base_colors"].get(mod, styles["base_colors"]["DNA"]),
                    styles["base_text_colors"].get(mod, "#111827"),
                    self._bubble_has_5mec(grid_value, modification),
                )
                centers[local_idx] = (x, y)
            for local_idx, (x, y) in centers.items():
                if local_idx >= len(result.chemistry.linkages) or local_idx + 1 not in centers:
                    continue
                linkage = result.chemistry.linkages[local_idx]
                diamond(x + layout["radius"], y, layout["diamond_r"], styles["linkage_colors"].get(linkage, "#ffffff"))

        self._draw_bubble_png_legend(
            draw,
            layout["legend_x"],
            layout["legend_y"],
            used_base_mods,
            used_linkages,
            styles,
            legend_title_font,
            legend_font,
            scale,
        )

        output_path.parent.mkdir(parents=True, exist_ok=True)
        image.save(output_path, dpi=(BUBBLE_EXPORT_DPI, BUBBLE_EXPORT_DPI))
        return output_path

    def save_penalty_bubble_figure_image(self) -> None:
        if self.last_penalty_result is None:
            self.messagebox.showerror(APP_NAME, "Generate penalty ASOs before saving a visualisation.")
            return
        path = self.filedialog.asksaveasfilename(
            title="Save high-resolution penalty ASO visualisation",
            defaultextension=".png",
            filetypes=[("PNG image", "*.png")],
            initialfile="penalty_aso_visualisation_600dpi.png",
        )
        if not path:
            return
        try:
            self._save_penalty_bubble_figure_png(self.last_penalty_result, Path(path))
        except ModuleNotFoundError as exc:
            if exc.name == "PIL":
                self.messagebox.showerror(
                    APP_NAME,
                    "PNG export requires Pillow. Install requirements again, then reopen the app.",
                )
            else:
                self.messagebox.showerror(APP_NAME, str(exc))
            return
        except Exception as exc:
            self.messagebox.showerror(APP_NAME, str(exc))
            return
        self.messagebox.showinfo(APP_NAME, f"Saved high-resolution PNG:\n{path}")

    def _save_penalty_bubble_figure_png(self, result, output_path: Path, show_letters: bool | None = None) -> Path:
        from PIL import Image, ImageDraw, ImageFont

        layout = self._penalty_bubble_layout(result)
        styles = self._bubble_styles()
        used_base_mods, used_linkages = self._penalty_used_modifications(result)
        if show_letters is None:
            show_letters = self._penalty_bubble_base_letters_visible()
        scale = BUBBLE_EXPORT_SCALE
        width = int(layout["width"] * scale)
        height = int(layout["height"] * scale)
        image = Image.new("RGB", (width, height), "white")
        draw = ImageDraw.Draw(image)

        def s(value: float) -> int:
            return int(round(float(value) * scale))

        def font(size: int, bold: bool = False):
            scaled_size = size * scale
            candidates = [
                "/System/Library/Fonts/Supplemental/Arial Bold.ttf" if bold else "/System/Library/Fonts/Supplemental/Arial.ttf",
                "arialbd.ttf" if bold else "arial.ttf",
                "DejaVuSans-Bold.ttf" if bold else "DejaVuSans.ttf",
            ]
            for candidate in candidates:
                try:
                    return ImageFont.truetype(candidate, scaled_size)
                except Exception:
                    continue
            return ImageFont.load_default()

        label_font = font(13, True)
        row_label_font = font(12, True)
        base_font = font(12, True)
        pos_font = font(10, False)
        legend_title_font = font(13, True)
        legend_font = font(12, False)

        def text_right(x: float, y: float, value: str, fnt, fill: str = "#1f4e79") -> None:
            bbox = draw.textbbox((0, 0), value, font=fnt)
            draw.text((s(x) - (bbox[2] - bbox[0]), s(y) - (bbox[3] - bbox[1]) / 2), value, font=fnt, fill=fill)

        def text_center(x: float, y: float, value: str, fnt, fill: str = "#111827") -> None:
            bbox = draw.textbbox((0, 0), value, font=fnt)
            draw.text(
                (s(x) - (bbox[2] - bbox[0]) / 2, s(y) - (bbox[3] - bbox[1]) / 2 - scale),
                value,
                font=fnt,
                fill=fill,
            )

        def circle(
            x: float,
            y: float,
            radius: float,
            value: str,
            fill: str,
            text_fill: str,
            methyl_marker: bool = False,
            penalty_marker: bool = False,
        ) -> None:
            outline = "#c00000" if penalty_marker else "#111827"
            draw.ellipse(
                (s(x - radius), s(y - radius), s(x + radius), s(y + radius)),
                fill=fill,
                outline=outline,
                width=s(3 if penalty_marker else 2),
            )
            if value:
                text_center(x, y, value, base_font, text_fill)
            if methyl_marker:
                marker_radius = max(3, radius * 0.34)
                draw.ellipse(
                    (
                        s(x + radius * 0.32),
                        s(y - radius * 0.92),
                        s(x + radius * 0.32 + marker_radius),
                        s(y - radius * 0.92 + marker_radius),
                    ),
                    fill=styles["methyl_marker"],
                    outline="#111827",
                    width=max(1, s(1)),
                )

        def diamond(x: float, y: float, radius: float, fill: str) -> None:
            points = [(s(x), s(y - radius)), (s(x + radius), s(y)), (s(x), s(y + radius)), (s(x - radius), s(y))]
            draw.polygon(points, fill=fill)
            draw.line(points + [points[0]], fill="#111827", width=max(1, s(1.3)))

        text_right(layout["left_label_x"], layout["rna_y"], "RNA 3' - 5'", label_font)
        for idx, position in enumerate(result.header_positions):
            x = layout["grid_left"] + idx * layout["pitch"] + layout["radius"]
            text_center(x, layout["position_y"], str(position), pos_font)
            circle(
                x,
                layout["rna_y"],
                layout["radius"],
                result.header_bases[idx] if show_letters else "",
                styles["rna_fill"],
                "#111827",
            )

        text_right(layout["left_label_x"], layout["direction_y"], "Penalty ASOs 5' - 3'", label_font)
        for row in result.rows:
            y = layout["row_top"] + (row.row_number - 1) * layout["row_h"]
            text_right(layout["left_label_x"], y, row.penalty_aso_id, row_label_font)
            centers: dict[int, tuple[float, float]] = {}
            for idx, grid_value in enumerate(row.grid_cells):
                if grid_value == "##":
                    continue
                local_idx = result.header_positions[idx] - row.parent_start
                if not 0 <= local_idx < len(result.chemistry.base_modifications):
                    continue
                x = layout["grid_left"] + idx * layout["pitch"] + layout["radius"]
                modification = result.chemistry.base_modifications[local_idx]
                mod = self._bubble_ribose_mod(modification)
                is_penalty = idx == row.penalty_grid_index
                circle(
                    x,
                    y,
                    layout["radius"],
                    grid_value if show_letters else "",
                    "#f6c343" if is_penalty else styles["base_colors"].get(mod, styles["base_colors"]["DNA"]),
                    "#111827" if is_penalty else styles["base_text_colors"].get(mod, "#111827"),
                    self._bubble_has_5mec(grid_value, modification),
                    is_penalty,
                )
                centers[local_idx] = (x, y)
            for local_idx, (x, y) in centers.items():
                if local_idx >= len(result.chemistry.linkages) or local_idx + 1 not in centers:
                    continue
                linkage = result.chemistry.linkages[local_idx]
                diamond(x + layout["radius"], y, layout["diamond_r"], styles["linkage_colors"].get(linkage, "#ffffff"))

        legend_bottom = self._draw_bubble_png_legend(
            draw,
            layout["legend_x"],
            layout["legend_y"],
            used_base_mods,
            used_linkages,
            styles,
            legend_title_font,
            legend_font,
            scale,
        )
        draw.ellipse(
            (
                s(layout["legend_x"]),
                s(legend_bottom - 9),
                s(layout["legend_x"] + 18),
                s(legend_bottom + 9),
            ),
            fill="#f6c343",
            outline="#c00000",
            width=max(1, s(2)),
        )
        draw.text(
            (s(layout["legend_x"] + 28), s(legend_bottom - 9)),
            "Intentional penalty base",
            font=legend_font,
            fill="#111827",
        )

        output_path.parent.mkdir(parents=True, exist_ok=True)
        image.save(output_path, dpi=(BUBBLE_EXPORT_DPI, BUBBLE_EXPORT_DPI))
        return output_path

    def _draw_bubble_png_legend(
        self,
        draw,
        x: float,
        y: float,
        used_bases: list[str],
        used_linkages: list[str],
        styles: dict,
        title_font,
        legend_font,
        scale: int = 1,
    ) -> float:
        def s(value: float) -> int:
            return int(round(float(value) * scale))

        draw.text((s(x), s(y)), "Legend", font=title_font, fill="#111827")
        cursor_y = y + 32
        line_h = 30
        for mod in used_bases:
            if mod == "5MeC":
                label = "5MeC base modification"
                draw.ellipse(
                    (s(x + 5), s(cursor_y - 5), s(x + 15), s(cursor_y + 5)),
                    fill=styles["methyl_marker"],
                    outline="#111827",
                    width=max(1, s(2)),
                )
            else:
                label = mod
                draw.ellipse(
                    (s(x), s(cursor_y - 9), s(x + 18), s(cursor_y + 9)),
                    fill=styles["base_colors"].get(mod, styles["base_colors"]["DNA"]),
                    outline="#111827",
                    width=max(1, s(2)),
                )
            draw.text((s(x + 28), s(cursor_y - 9)), label, font=legend_font, fill="#111827")
            cursor_y += line_h
        for linkage in used_linkages:
            label = f"{linkage} linkage"
            tri_x = x + 9
            tri_y = cursor_y
            points = [
                (s(tri_x), s(tri_y - 8)),
                (s(tri_x + 8), s(tri_y)),
                (s(tri_x), s(tri_y + 8)),
                (s(tri_x - 8), s(tri_y)),
            ]
            draw.polygon(points, fill=styles["linkage_colors"].get(linkage, "#ffffff"))
            draw.line(points + [points[0]], fill="#111827", width=max(1, s(1.3)))
            draw.text((s(x + 28), s(cursor_y - 9)), label, font=legend_font, fill="#111827")
            cursor_y += line_h
        return cursor_y + 8

    def export_excel(self) -> None:
        self.calculate(show_errors=True)
        if self.last_result is None:
            return
        path = self.filedialog.asksaveasfilename(
            title="Export ASO workbook",
            defaultextension=".xlsx",
            filetypes=[("Excel workbook", "*.xlsx")],
            initialfile="aso_design_output.xlsx",
        )
        if not path:
            return
        try:
            from excel_export import export_result_to_xlsx

            output = export_result_to_xlsx(self.last_result, path)
        except Exception as exc:
            self.messagebox.showerror(APP_NAME, str(exc))
            return
        self.messagebox.showinfo(APP_NAME, f"Exported {output}")


def main(argv: list[str] | None = None) -> int:
    args = parse_args(sys.argv[1:] if argv is None else argv)
    if args.export:
        try:
            from excel_export import export_result_to_xlsx

            result = generate_design(inputs_from_args(args))
            output = export_result_to_xlsx(result, args.export)
        except AsoInputError as exc:
            print(f"error: {exc}", file=sys.stderr)
            return 2
        print(output)
        return 0

    app = AsoDesignerApp()
    app.run()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
