import tkinter as tk
from tkinter import ttk, filedialog, scrolledtext, messagebox
import io

# Biopython imports
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Align import substitution_matrices

# ─────────────────────────────────────────────
#  ALIGNMENT FUNCTIONS (powered by Biopython)
# ─────────────────────────────────────────────

def get_aligner(mode, matrix_name, gap):
    """Build and return a configured Biopython PairwiseAligner."""
    aligner = PairwiseAligner()
    aligner.mode = mode  # 'global' or 'local'

    if matrix_name == "BLOSUM62":
        aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    elif matrix_name == "PAM250":
        aligner.substitution_matrix = substitution_matrices.load("PAM250")
    else:
        # Simple match/mismatch scoring
        aligner.match_score = 1
        aligner.mismatch_score = -1

    aligner.open_gap_score = gap
    aligner.extend_gap_score = gap
    return aligner


def pairwise_align(seq1, seq2, mode, matrix_name, gap):
    """Align two sequences and return (aligned1, aligned2, score)."""
    aligner = get_aligner(mode, matrix_name, gap)
    alignments = aligner.align(seq1, seq2)
    best = next(iter(alignments))  # take the top alignment

    # Extract the two aligned strings (with gaps)
    aligned_str = str(best)
    lines = aligned_str.strip().split("\n")

    # Biopython formats as: seq1_line, middle_line, seq2_line (repeating blocks)
    top_lines = []
    bot_lines = []
    for i, line in enumerate(lines):
        mod = i % 3
        if mod == 0:
            top_lines.append(line)
        elif mod == 2:
            bot_lines.append(line)

    a1 = "".join(top_lines)
    a2 = "".join(bot_lines)
    return a1, a2, best.score


def progressive_msa(sequences, matrix_name, gap):
    """
    ClustalW-style progressive MSA using Biopython for pairwise steps.
    1. Align seq1 vs seq2
    2. Build consensus, align next sequence to it
    3. Expand gaps into the whole alignment
    """
    if len(sequences) < 2:
        return sequences

    a1, a2, _ = pairwise_align(sequences[0], sequences[1], "global", matrix_name, gap)
    aligned = [a1, a2]

    for seq in sequences[2:]:
        # Consensus = most common non-gap base per column
        consensus = ""
        for col in range(len(aligned[0])):
            bases = [s[col] for s in aligned if col < len(s) and s[col] != "-"]
            consensus += bases[0] if bases else "-"

        new_cons, new_seq, _ = pairwise_align(consensus, seq, "global", matrix_name, gap)

        # Re-expand existing aligned rows to match new gaps inserted into consensus
        expanded = [""] * len(aligned)
        ci = 0
        for char in new_cons:
            if char == "-":
                for row in range(len(aligned)):
                    expanded[row] += "-"
            else:
                for row in range(len(aligned)):
                    expanded[row] += aligned[row][ci] if ci < len(aligned[row]) else "-"
                ci += 1

        max_len = max(len(e) for e in expanded)
        for row in range(len(expanded)):
            expanded[row] = expanded[row].ljust(max_len, "-")

        aligned = expanded + [new_seq.ljust(max_len, "-")]

    return aligned


def get_stats(aligned_seqs):
    if len(aligned_seqs) < 2:
        return ""
    a, b = aligned_seqs[0], aligned_seqs[1]
    length = len(a)
    matches = sum(1 for x, y in zip(a, b) if x == y and x != "-")
    gaps = sum(1 for x in a + b if x == "-")
    identity = (matches / length * 100) if length else 0
    return f"Length: {length}  |  Matches: {matches}  |  Identity: {identity:.1f}%  |  Gaps: {gaps}"


# ─────────────────────────────────────────────
#  SEQUENCE PARSER (uses Biopython SeqIO)
# ─────────────────────────────────────────────

def parse_input(text):
    """
    Accepts:
      - FASTA format (detected by '>')
      - Plain text: one sequence per line
    Returns (sequences_list, names_list)
    """
    text = text.strip()
    sequences, names = [], []

    if ">" in text:
        # Use Biopython's SeqIO to parse FASTA
        handle = io.StringIO(text)
        for record in SeqIO.parse(handle, "fasta"):
            sequences.append(str(record.seq).upper())
            names.append(record.id)
    else:
        # Plain text: one sequence per line
        for i, line in enumerate(text.split("\n")):
            seq = line.strip().upper()
            if seq:
                sequences.append(seq)
                names.append(f"Seq{i+1}")

    return sequences, names


# ─────────────────────────────────────────────
#  GUI
# ─────────────────────────────────────────────

class AlignerApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("🧬 Sequence Aligner  (powered by Biopython)")
        self.geometry("1080x700")
        self.configure(bg="#f4f4f4")
        self.resizable(True, True)
        self._build()

    def _build(self):
        # ── TOP BAR ──────────────────────────────
        top = tk.Frame(self, bg="#2c3e50", pady=10)
        top.pack(fill=tk.X)
        tk.Label(top, text="🧬  Sequence Aligner", font=("Helvetica", 18, "bold"),
                 bg="#2c3e50", fg="white").pack()
        tk.Label(top, text="Global · Local · Multiple Sequence Alignment  —  powered by Biopython",
                 font=("Helvetica", 10), bg="#2c3e50", fg="#bdc3c7").pack()

        # ── MAIN FRAME ───────────────────────────
        main = tk.Frame(self, bg="#f4f4f4")
        main.pack(fill=tk.BOTH, expand=True, padx=14, pady=10)

        # ── LEFT PANEL ───────────────────────────
        # Split into scrollable top options + fixed bottom buttons
        left = tk.Frame(main, bg="#f4f4f4", width=310)
        left.pack(side=tk.LEFT, fill=tk.Y, padx=(0, 12))
        left.pack_propagate(False)

        # -- BOTTOM buttons pinned first so they never get pushed off screen
        btn_frame = tk.Frame(left, bg="#f4f4f4")
        btn_frame.pack(side=tk.BOTTOM, fill=tk.X, pady=(6, 0))
        tk.Button(btn_frame, text="▶  Run Alignment", font=("Helvetica", 12, "bold"),
                  bg="#27ae60", fg="white", relief="flat", pady=8,
                  command=self._run).pack(fill=tk.X)
        tk.Button(btn_frame, text="🗑  Clear", bg="#e74c3c", fg="white",
                  relief="flat", pady=4, command=self._clear).pack(fill=tk.X, pady=(4, 0))
        tk.Button(btn_frame, text="💾  Export result", bg="#3498db", fg="white",
                  relief="flat", pady=4, command=self._export).pack(fill=tk.X, pady=(4, 0))

        # -- Scrollable options area above the buttons
        canvas = tk.Canvas(left, bg="#f4f4f4", highlightthickness=0)
        scrollbar = ttk.Scrollbar(left, orient="vertical", command=canvas.yview)
        options = tk.Frame(canvas, bg="#f4f4f4")
        options.bind("<Configure>",
                     lambda e: canvas.configure(scrollregion=canvas.bbox("all")))
        canvas.create_window((0, 0), window=options, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)
        canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

        # Mousewheel scrolling
        canvas.bind_all("<MouseWheel>", lambda e: canvas.yview_scroll(-1*(e.delta//120), "units"))

        # Now all options go into `options` frame instead of `left`
        tk.Label(options, text="Enter sequences:", font=("Helvetica", 11, "bold"),
                 bg="#f4f4f4").pack(anchor="w")
        tk.Label(options, text="FASTA format  or  one sequence per line",
                 font=("Helvetica", 9), bg="#f4f4f4", fg="#666").pack(anchor="w")

        self.seq_input = scrolledtext.ScrolledText(
            options, width=34, height=9, font=("Courier", 11), relief="solid", bd=1)
        self.seq_input.pack(fill=tk.X, pady=(4, 6))

        tk.Button(options, text="Load example sequences", bg="#ecf0f1",
                  relief="flat", command=self._load_example).pack(fill=tk.X, pady=2)
        tk.Button(options, text="📂  Open FASTA / TXT file", bg="#ecf0f1",
                  relief="flat", command=self._load_file).pack(fill=tk.X, pady=2)

        ttk.Separator(options, orient="horizontal").pack(fill=tk.X, pady=10)

        # Algorithm
        tk.Label(options, text="Algorithm:", font=("Helvetica", 10, "bold"),
                 bg="#f4f4f4").pack(anchor="w")
        self.algo_var = tk.StringVar(value="Global (Needleman-Wunsch)")
        for algo in ["Global (Needleman-Wunsch)",
                     "Local (Smith-Waterman)",
                     "Multiple Sequence Alignment"]:
            tk.Radiobutton(options, text=algo, variable=self.algo_var, value=algo,
                           bg="#f4f4f4", font=("Helvetica", 10)).pack(anchor="w")

        ttk.Separator(options, orient="horizontal").pack(fill=tk.X, pady=10)

        # Scoring matrix
        tk.Label(options, text="Scoring matrix:", font=("Helvetica", 10, "bold"),
                 bg="#f4f4f4").pack(anchor="w")
        self.matrix_var = tk.StringVar(value="Simple match/mismatch")
        matrix_menu = ttk.Combobox(options, textvariable=self.matrix_var,
                                    values=["Simple match/mismatch", "BLOSUM62", "PAM250"],
                                    state="readonly")
        matrix_menu.pack(fill=tk.X, pady=(2, 8))

        # Gap penalty
        gap_row = tk.Frame(options, bg="#f4f4f4")
        gap_row.pack(fill=tk.X)
        tk.Label(gap_row, text="Gap penalty:", bg="#f4f4f4",
                 font=("Helvetica", 10)).pack(side=tk.LEFT)
        self.gap_var = tk.StringVar(value="-2")
        tk.Entry(gap_row, textvariable=self.gap_var, width=6,
                 relief="solid").pack(side=tk.LEFT, padx=6)

        # ── RIGHT PANEL ──────────────────────────
        right = tk.Frame(main, bg="#f4f4f4")
        right.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

        tk.Label(right, text="Alignment output:", font=("Helvetica", 11, "bold"),
                 bg="#f4f4f4").pack(anchor="w")

        # Color legend
        legend = tk.Frame(right, bg="#f4f4f4")
        legend.pack(anchor="w", pady=(0, 4))
        tk.Label(legend, text="Color key:", font=("Helvetica", 9),
                 bg="#f4f4f4").pack(side=tk.LEFT, padx=(0, 6))
        for base, color in [("A", "#27ae60"), ("T/U", "#e74c3c"),
                              ("G", "#e67e22"), ("C", "#2980b9"), ("Gap -", "#888")]:
            tk.Label(legend, text=f"  {base}  ", bg=color, fg="white",
                     font=("Courier", 9, "bold"), padx=3).pack(side=tk.LEFT, padx=2)

        # Output text widget
        out_frame = tk.Frame(right)
        out_frame.pack(fill=tk.BOTH, expand=True)

        xscroll = tk.Scrollbar(out_frame, orient=tk.HORIZONTAL)
        yscroll = tk.Scrollbar(out_frame, orient=tk.VERTICAL)
        self.output = tk.Text(
            out_frame, font=("Courier", 13), wrap=tk.NONE,
            xscrollcommand=xscroll.set, yscrollcommand=yscroll.set,
            state=tk.DISABLED, relief="solid", bd=1,
            bg="#1e1e1e", fg="white")
        xscroll.config(command=self.output.xview)
        yscroll.config(command=self.output.yview)
        yscroll.pack(side=tk.RIGHT, fill=tk.Y)
        xscroll.pack(side=tk.BOTTOM, fill=tk.X)
        self.output.pack(fill=tk.BOTH, expand=True)

        # Stats bar
        self.stats_var = tk.StringVar(value="Run an alignment to see stats here.")
        tk.Label(right, textvariable=self.stats_var, font=("Helvetica", 10),
                 bg="#2c3e50", fg="white", anchor="w", pady=4,
                 padx=8).pack(fill=tk.X, pady=(6, 0))

        # Register color tags on the output widget
        BASE_COLORS = {"A": "#27ae60", "T": "#e74c3c", "U": "#e74c3c",
                       "G": "#e67e22", "C": "#2980b9", "-": "#666666"}
        for base, color in BASE_COLORS.items():
            self.output.tag_config(base, foreground=color)
        self.output.tag_config("label", foreground="#aaaaaa",
                               font=("Courier", 13, "bold"))
        self.output.tag_config("cons", foreground="#f1c40f")

    # ── EVENT HANDLERS ───────────────────────

    def _load_example(self):
        self.seq_input.delete("1.0", tk.END)
        self.seq_input.insert(tk.END,
            ">Human\nATGCGTACGTTAGCTAGCTTACG\n"
            ">Mouse\nATGCGTACGTTGGCTAACTTACG\n"
            ">Zebrafish\nATGCGTAGGTTAGCTAGCTTGCG\n"
        )

    def _load_file(self):
        path = filedialog.askopenfilename(
            filetypes=[("FASTA / Text", "*.fasta *.fa *.txt"), ("All files", "*.*")])
        if not path:
            return
        with open(path) as f:
            self.seq_input.delete("1.0", tk.END)
            self.seq_input.insert(tk.END, f.read())

    def _run(self):
        raw = self.seq_input.get("1.0", tk.END)
        seqs, names = parse_input(raw)

        if len(seqs) < 2:
            messagebox.showwarning("Need more sequences",
                                   "Please enter at least 2 sequences.")
            return

        try:
            gap = int(self.gap_var.get())
        except ValueError:
            messagebox.showerror("Bad gap value", "Gap penalty must be an integer (e.g. -2).")
            return

        algo = self.algo_var.get()
        matrix = self.matrix_var.get()
        aligned = []
        score_text = ""

        try:
            if "Multiple" in algo:
                aligned = progressive_msa(seqs, matrix, gap)
                names_out = names[:len(aligned)]
            elif "Local" in algo:
                if len(seqs) > 2:
                    messagebox.showinfo("Note", "Local alignment uses the first 2 sequences.")
                a1, a2, score = pairwise_align(seqs[0], seqs[1], "local", matrix, gap)
                aligned = [a1, a2]
                names_out = names[:2]
                score_text = f"  |  Score: {score:.1f}"
            else:
                if len(seqs) > 2:
                    messagebox.showinfo("Note", "Global pairwise uses the first 2 sequences.")
                a1, a2, score = pairwise_align(seqs[0], seqs[1], "global", matrix, gap)
                aligned = [a1, a2]
                names_out = names[:2]
                score_text = f"  |  Score: {score:.1f}"
        except Exception as e:
            messagebox.showerror("Alignment error", str(e))
            return

        self._display(aligned, names_out)
        self.stats_var.set(get_stats(aligned) + score_text)

    def _display(self, aligned_seqs, names):
        KNOWN_BASES = {"A", "T", "U", "G", "C", "-"}

        self.output.config(state=tk.NORMAL)
        self.output.delete("1.0", tk.END)

        label_w = max((len(n) for n in names), default=4) + 2

        for i, seq in enumerate(aligned_seqs):
            name = names[i] if i < len(names) else f"Seq{i+1}"
            self.output.insert(tk.END, f"{name:<{label_w}} ", "label")
            for char in seq:
                tag = char if char in KNOWN_BASES else "A"
                self.output.insert(tk.END, char, tag)
            self.output.insert(tk.END, "\n")

        # Conservation row for pairwise
        if len(aligned_seqs) == 2:
            cons = ""
            for a, b in zip(aligned_seqs[0], aligned_seqs[1]):
                cons += "|" if a == b and a != "-" else ("." if a != "-" and b != "-" else " ")
            self.output.insert(tk.END, " " * (label_w + 1) + cons + "\n", "cons")

        self.output.config(state=tk.DISABLED)

    def _clear(self):
        self.seq_input.delete("1.0", tk.END)
        self.output.config(state=tk.NORMAL)
        self.output.delete("1.0", tk.END)
        self.output.config(state=tk.DISABLED)
        self.stats_var.set("Run an alignment to see stats here.")

    def _export(self):
        content = self.output.get("1.0", tk.END)
        if not content.strip():
            messagebox.showinfo("Nothing to export", "Run an alignment first.")
            return
        path = filedialog.asksaveasfilename(
            defaultextension=".txt",
            filetypes=[("Text file", "*.txt"), ("All files", "*.*")])
        if path:
            with open(path, "w") as f:
                f.write(content)
            messagebox.showinfo("Saved", f"Result saved to:\n{path}")


# ─────────────────────────────────────────────
if __name__ == "__main__":
    app = AlignerApp()
    app.mainloop()
