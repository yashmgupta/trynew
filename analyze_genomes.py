#!/usr/bin/env python3
"""
Genome analysis tool to compute base composition statistics, sliding-window
entropy/complexity metrics, and produce comparative visualizations.

Features
--------
* Parses FASTA or GenBank flat files without external dependencies.
* Computes genome-wide GC content and AT/GC skew values.
* Aggregates per-position nucleotide frequencies across genomes.
* Calculates sliding-window Shannon entropy and sequence complexity.
* Generates skew trajectory plots, entropy heatmaps, and cross-species
  comparative box plots.

Example
-------
python analyze_genomes.py --output-dir results --window-size 500 --step 100 \
    Chilopoda.gb other_genome.fasta
"""
from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence

try:  # Optional dependency: plotting is skipped if matplotlib is unavailable.
    import matplotlib.pyplot as plt
    HAS_MATPLOTLIB = True
except Exception:
    HAS_MATPLOTLIB = False


@dataclass
class GenomeRecord:
    name: str
    sequence: str

    @property
    def length(self) -> int:
        return len(self.sequence)


BASES = ("A", "T", "G", "C")


def parse_sequence(path: Path) -> GenomeRecord:
    """Parse FASTA or GenBank files without external dependencies.

    The genome name is taken from the FASTA header (first non-empty line after
    ">") or the LOCUS line in a GenBank file. Sequence characters are coerced to
    uppercase and non-ATGC symbols are ignored.
    """

    text = path.read_text()
    name = path.stem

    if "\n>" in text or text.lstrip().startswith(">"):  # FASTA-like
        header, *rest = text.splitlines()
        if header.startswith(">") and header[1:].strip():
            name = header[1:].split()[0]
        sequence = "".join(line.strip() for line in rest if not line.startswith(">"))
    else:  # GenBank flat file
        seq_lines: List[str] = []
        seen_origin = False
        for line in text.splitlines():
            if line.startswith("LOCUS"):
                parts = line.split()
                if len(parts) >= 2:
                    name = parts[1]
            if line.strip().upper() == "ORIGIN":
                seen_origin = True
                continue
            if seen_origin:
                if line.strip() == "//":
                    break
                seq_lines.append("".join(ch for ch in line if ch.isalpha()))
        sequence = "".join(seq_lines)

    sequence = "".join(ch for ch in sequence.upper() if ch in BASES)
    if not sequence:
        raise ValueError(f"No sequence data found in {path}")

    return GenomeRecord(name=name, sequence=sequence)


def base_counts(sequence: str) -> Dict[str, int]:
    return {base: sequence.count(base) for base in BASES}


def compute_global_metrics(genomes: Sequence[GenomeRecord]) -> List[Dict[str, float]]:
    rows: List[Dict[str, float]] = []
    for genome in genomes:
        counts = base_counts(genome.sequence)
        total = sum(counts.values())
        gc = counts["G"] + counts["C"]
        at = counts["A"] + counts["T"]
        gc_content = gc / total if total else 0.0
        at_skew = (counts["A"] - counts["T"]) / at if at else 0.0
        gc_skew = (counts["G"] - counts["C"]) / gc if gc else 0.0
        rows.append(
            {
                "genome": genome.name,
                "length": total,
                "gc_content": gc_content,
                "at_skew": at_skew,
                "gc_skew": gc_skew,
            }
        )
    return rows


def write_tsv(path: Path, headers: Sequence[str], rows: Iterable[Sequence[object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(headers)
        for row in rows:
            writer.writerow(row)


def compute_per_position_frequencies(genomes: Sequence[GenomeRecord]) -> List[List[float]]:
    if not genomes:
        return []
    max_len = max(g.length for g in genomes)
    frequencies: List[List[float]] = []
    for idx in range(max_len):
        base_counts_at_pos = {base: 0 for base in BASES}
        genomes_covering = 0
        for genome in genomes:
            if idx < genome.length:
                base_counts_at_pos[genome.sequence[idx]] += 1
                genomes_covering += 1
        if genomes_covering == 0:
            continue
        frequencies.append(
            [
                idx + 1,
                *(
                    base_counts_at_pos[base] / genomes_covering
                    for base in BASES
                ),
                genomes_covering,
            ]
        )
    return frequencies


def shannon_entropy(counts: Dict[str, int]) -> float:
    total = sum(counts.values())
    if total == 0:
        return 0.0
    entropy = 0.0
    for count in counts.values():
        if count == 0:
            continue
        p = count / total
        entropy -= p * math.log2(p)
    return entropy


def sequence_complexity(counts: Dict[str, int]) -> float:
    """Return a diversity-like complexity metric in [0, 1].

    Defined as 1 - sum(p_i^2), analogous to the Gini-Simpson index.
    """

    total = sum(counts.values())
    if total == 0:
        return 0.0
    return 1.0 - sum((c / total) ** 2 for c in counts.values())


def sliding_window_metrics(genome: GenomeRecord, window_size: int, step: int) -> List[Dict[str, float]]:
    metrics: List[Dict[str, float]] = []
    seq = genome.sequence
    for start in range(0, len(seq) - window_size + 1, step):
        window = seq[start : start + window_size]
        counts = base_counts(window)
        entropy = shannon_entropy(counts)
        complexity = sequence_complexity(counts)
        metrics.append(
            {
                "genome": genome.name,
                "start": start + 1,
                "end": start + window_size,
                "entropy": entropy,
                "complexity": complexity,
            }
        )
    return metrics


def plot_skew(genome: GenomeRecord, out_path: Path) -> None:
    if not HAS_MATPLOTLIB:
        return
    gc_walk_x: List[int] = [0]
    gc_walk_y: List[int] = [0]
    at_walk_x: List[int] = [0]
    at_walk_y: List[int] = [0]

    for base in genome.sequence:
        gc_step = 1 if base == "G" else -1 if base == "C" else 0
        at_step = 1 if base == "A" else -1 if base == "T" else 0
        gc_walk_x.append(gc_walk_x[-1] + 1)
        gc_walk_y.append(gc_walk_y[-1] + gc_step)
        at_walk_x.append(at_walk_x[-1] + 1)
        at_walk_y.append(at_walk_y[-1] + at_step)

    fig, axes = plt.subplots(2, 1, figsize=(8, 8), constrained_layout=True)
    axes[0].plot(gc_walk_x, gc_walk_y, label="Cumulative GC skew", color="#1b9e77")
    axes[0].axhline(0, color="black", linewidth=0.8, linestyle="--")
    axes[0].set_title(f"GC skew random-walk: {genome.name}")
    axes[0].set_xlabel("Position (bp)")
    axes[0].set_ylabel("G - C cumulative")

    axes[1].plot(at_walk_x, at_walk_y, label="Cumulative AT skew", color="#d95f02")
    axes[1].axhline(0, color="black", linewidth=0.8, linestyle="--")
    axes[1].set_title(f"AT skew random-walk: {genome.name}")
    axes[1].set_xlabel("Position (bp)")
    axes[1].set_ylabel("A - T cumulative")

    for ax in axes:
        ax.legend()

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def plot_entropy_heatmap(metrics: List[Dict[str, float]], genome: GenomeRecord, out_path: Path) -> None:
    if not metrics or not HAS_MATPLOTLIB:
        return
    positions = [m["start"] for m in metrics]
    entropy = [m["entropy"] for m in metrics]
    complexity = [m["complexity"] for m in metrics]

    data = [entropy, complexity]
    fig, ax = plt.subplots(figsize=(10, 3), constrained_layout=True)
    cax = ax.imshow(
        [data_i for data_i in data],
        aspect="auto",
        extent=[positions[0], positions[-1], 0, 2],
        cmap="viridis",
        origin="lower",
    )
    ax.set_yticks([0.5, 1.5])
    ax.set_yticklabels(["Entropy", "Complexity"])
    ax.set_xlabel("Genome position (bp)")
    ax.set_title(f"Sliding-window metrics: {genome.name}")
    fig.colorbar(cax, label="Metric value")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def plot_comparative_box(global_metrics: List[Dict[str, float]], out_path: Path) -> None:
    if not global_metrics or not HAS_MATPLOTLIB:
        return
    fig, axes = plt.subplots(1, 3, figsize=(12, 4), constrained_layout=True)

    gc_values = [row["gc_content"] for row in global_metrics]
    at_skews = [row["at_skew"] for row in global_metrics]
    gc_skews = [row["gc_skew"] for row in global_metrics]

    axes[0].boxplot(gc_values, vert=True)
    axes[0].set_title("GC content")
    axes[0].set_ylabel("Fraction")
    axes[0].set_xticklabels(["Genomes"])

    axes[1].violinplot(at_skews, showmeans=True)
    axes[1].set_title("AT skew")
    axes[1].set_ylabel("(A-T)/(A+T)")
    axes[1].set_xticks([1])
    axes[1].set_xticklabels(["Genomes"])

    axes[2].violinplot(gc_skews, showmeans=True)
    axes[2].set_title("GC skew")
    axes[2].set_ylabel("(G-C)/(G+C)")
    axes[2].set_xticks([1])
    axes[2].set_xticklabels(["Genomes"])

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Genome composition and entropy analyzer")
    parser.add_argument("genomes", nargs="+", type=Path, help="Input genome files (FASTA or GenBank)")
    parser.add_argument("--output-dir", type=Path, default=Path("analysis_output"), help="Destination for tables and plots")
    parser.add_argument("--window-size", type=int, default=500, help="Sliding-window size (bp)")
    parser.add_argument("--step", type=int, default=100, help="Sliding-window step size (bp)")
    args = parser.parse_args()

    genomes = [parse_sequence(path) for path in args.genomes]

    if not HAS_MATPLOTLIB:
        print("matplotlib not available; plots will be skipped. Install via 'pip install -r requirements.txt'.")

    global_metrics = compute_global_metrics(genomes)
    write_tsv(
        args.output_dir / "global_metrics.tsv",
        headers=["genome", "length", "gc_content", "at_skew", "gc_skew"],
        rows=[
            [row["genome"], row["length"], f"{row['gc_content']:.6f}", f"{row['at_skew']:.6f}", f"{row['gc_skew']:.6f}"]
            for row in global_metrics
        ],
    )

    per_position = compute_per_position_frequencies(genomes)
    write_tsv(
        args.output_dir / "per_position_frequencies.tsv",
        headers=["position", "A_freq", "T_freq", "G_freq", "C_freq", "genomes_covered"],
        rows=[
            [row[0]] + [f"{value:.6f}" for value in row[1:-1]] + [row[-1]]
            for row in per_position
        ],
    )

    all_window_metrics: List[Dict[str, float]] = []
    for genome in genomes:
        metrics = sliding_window_metrics(genome, args.window_size, args.step)
        all_window_metrics.extend(metrics)
        write_tsv(
            args.output_dir / f"{genome.name}_sliding_window.tsv",
            headers=["genome", "start", "end", "entropy", "complexity"],
            rows=[
                [m["genome"], m["start"], m["end"], f"{m['entropy']:.6f}", f"{m['complexity']:.6f}"]
                for m in metrics
            ],
        )
        plot_skew(genome, args.output_dir / f"{genome.name}_skew.png")
        plot_entropy_heatmap(metrics, genome, args.output_dir / f"{genome.name}_entropy_heatmap.png")

    plot_comparative_box(global_metrics, args.output_dir / "comparative_metrics.png")

    # Combined sliding-window table
    if all_window_metrics:
        write_tsv(
            args.output_dir / "all_sliding_window.tsv",
            headers=["genome", "start", "end", "entropy", "complexity"],
            rows=[
                [m["genome"], m["start"], m["end"], f"{m['entropy']:.6f}", f"{m['complexity']:.6f}"]
                for m in all_window_metrics
            ],
        )


if __name__ == "__main__":
    main()
