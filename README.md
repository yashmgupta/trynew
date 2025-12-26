# Genome entropy and skew analyzer

This repository provides a lightweight Python script to compute genome-wide GC content, AT/GC skew, sliding-window Shannon entropy, and sequence-complexity metrics. It produces tabular outputs plus plots for skew trajectories, entropy heatmaps, and comparative genome-wide summaries.

## Requirements

Install dependencies (matplotlib only) into your environment. If matplotlib is
not available, the script will still emit TSV outputs but skip plotting.

```bash
pip install -r requirements.txt
```

## Usage

Run the analyzer against one or more FASTA or GenBank files. Sliding-window parameters are configurable.

```bash
python analyze_genomes.py --output-dir analysis_output --window-size 500 --step 100 \
    Chilopoda.gb
```

Outputs include:

- `global_metrics.tsv`: length, GC content, AT/GC skew per genome
- `per_position_frequencies.tsv`: nucleotide frequencies at each position across genomes
- `{genome}_sliding_window.tsv`: Shannon entropy and complexity for each window
- Plots: cumulative skew paths, entropy heatmaps, and comparative box/violin plots in the output directory

Sequence complexity is reported as `1 - sum(p_i^2)`, analogous to a Gini-Simpson diversity index, while entropy is Shannon entropy in bits.
