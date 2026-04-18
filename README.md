# STAT 623 Project — Shotgun Sequencing Simulation

This project simulates the shotgun sequencing process, builds an overlap graph for genome assembly, and visualizes the effects of repetitive elements on contig assembly and read mapping quality.

## Dependencies

Install the required Python packages via conda:

```bash
conda install -c conda-forge networkx matplotlib pysam numpy
```

| Package | Purpose |
|---------|---------|
| `pysam` | BAM file creation, read alignment I/O |
| `matplotlib` | All plotting (alignments, contigs, coverage, Lander-Waterman sweeps) |
| `numpy` | Theoretical curve computation for Lander-Waterman plots |
| `networkx` | Overlap graph construction and visualization |

> [!NOTE]
> `random`, `os`, and `math` are used from the Python standard library — no extra install needed.

## Project Files

| File | Description |
|------|-------------|
| `shortgun_simulation.py` | Main entry point — genome generation, read sampling, scenario runners, and Lander-Waterman sweep |
| `alignment_visualization.py` | Realistic alignment pipeline (`create_alignment_bam`), contig/coverage extraction, and all visualization functions |
| `theoretical_stats.py` | Lander-Waterman helper functions (`compute_alpha`, `expected_contigs`, `expected_contig_length`) |
| `overlap.py` | Overlap graph construction (`build_overlap_graph`) and NetworkX visualization (`draw_graph`) |

## Usage

Run all scenarios and Lander-Waterman sweeps:
```bash
python shortgun_simulation.py
```

Run the alignment visualization standalone test:
```bash
python alignment_visualization.py
```

---

## Configuration & Scenarios

### Scenario 1 — Human Mitochondrial Genome Assembly

**Function**: `run_alignment_comparison()`
Simulates a realistic biological scenario based on human mitochondrial DNA (~16.5 kb), which has high unique density and lacks large transpositional elements, focusing instead on hypervariable regions and tandem repeats.

| Parameter | Value |
|-----------|-------|
| Genome length ($G$) | 16 569 bp |
| Read length range ($L$) | 100 – 150 bp |
| Number of reads ($N$) | 2 500 (~18x coverage) |
| SNP rate | 0.01 |
| Repeat sizes | HVR = 150 bp |
| Repeat weights | HVR 3%, Tandem 5%, Unique 92% |

**Visualization** highlights MAPQ scores in simulated realistic conditions:
- **Steelblue** — Unique reads (MAPQ 60)
- **Orange** — Multi-mapped reads (MAPQ 20, ≤ 5 matches)
- **Light Gray** — Ambiguous/repetitive reads (MAPQ 0)

**Outputs**: `alignment_nr_scenario1.png`, `alignment_r_scenario1.png`, `contigs_s1_nr.png`, `contigs_s1_r.png`

---

### Scenario 2 — Contig Structure Comparison (Fragmentation)

**Function**: `run_assembly_comparison()`
Demonstrates how **long, frequent repeats** in the human nuclear genome (where nearly half the genome is repetitive) break assembly continuity and fragment contigs.

| Parameter | Value |
|-----------|-------|
| Genome length ($G$) | 4 000 bp |
| Read length range ($L$) | 50 – 100 bp |
| Number of reads ($N$) | 600 (~11x coverage) |
| SNP rate | 0.005 |
| Repeat sizes | LINE = 1 000 bp, SINE = 300 bp, LTR = 500 bp, DNA = 300 bp |
| Repeat weights | LINE 20%, SINE 13%, LTR 8%, DNA 3%, Unique 55% |

**Visualization**: Shows fragmented **red contig blocks**. Gaps typically appear within or at the boundaries of long repeats where unique overlaps cannot be established.

**Outputs**: `alignment_nr_scenario2.png`, `alignment_r_scenario2.png`, `contigs_s2_nr.png`, `contigs_s2_r.png`

---

### Scenario 3 — Read Length Sweep Analysis

**Function**: `run_read_length_sweep()`
Sweeps read lengths to demonstrate how longer reads bridge repeats and resolve fragmentation, while maintaining a constant coverage depth of 30x on a repeat-rich nuclear genome segment.

| Parameter | Value |
|-----------|-------|
| Genome length ($G$) | 4 000 bp |
| Target Coverage | 30x |
| Read Length Ranges | 50-75, 75-100, 150-200, 250-300, 350-400, 750-800 bp |
| Trials | 5 trials per read length |

**Visualization**: Shows the "stair-step" reduction in contig fragmentation as average read lengths sequentially exceed the sizes of SINEs (300 bp), LTRs (400 bp), and truncated LINEs (800 bp).

**Outputs**: `read_length_sweep_s2.png`

---

## Theoretical Background: Lander-Waterman Model

The simulation validates the **Lander-Waterman Theory** via `run_lander_waterman_sweep()`.

### Key Parameters

| Parameter | Description |
|-----------|-------------|
| Read density ($\alpha$) | $\alpha = N \times L_{\text{mean}} / G$ |
| $L_{\text{mean}}$ | $(L_{\min} + L_{\max}) / 2 = 75$ bp |
| Alpha sweep range | 0.2, 0.5, 0.8, 1.2, 1.5, 1.8, 2.2, 2.6, 3.0, 3.25, 3.6, 4.0, 4.4, 4.8, 5.2, 5.6, 6.0 |
| Multi-trial averaging | S1: 3 trials per $\alpha$,  S2: 5 trials per $\alpha$ |

### Theoretical Curve
$$E[\# \text{contigs}] = N e^{-\alpha} \approx \frac{G}{L} \alpha e^{-\alpha}$$

### LW Sweep Outputs (per scenario)

Each sweep produces a two-panel plot:
- **Left panel (No repeats)**: Empirical contig counts (blue stems) closely follow the theoretical dashed-gray curve.
- **Right panel (With repeats)**: Empirical contig counts (red stems) deviate above the theory, showing how repeat structure violates the Poisson assumption.

**Outputs**: `lander_waterman_sweep_s1.png`, `lander_waterman_sweep_s2.png`

---

## Implementation Details

1. **Genome Generation** — A bio-realistic library including LINEs, SINEs, LTRs, DNA transposons, tandem repeats, and unique segments. Each repeat family uses a single source template; copies are mutated at the specified `snp_rate`.
2. **Alignment Pipeline** — `find_with_mismatches` allows up to 2 mismatches per alignment position to simulate realistic mapper constraints.
3. **Probabilistic Read Dropping** — Highly ambiguous reads (MAPQ 0) are dropped with 90% probability; multi-mapped reads (MAPQ 20) are dropped with 80% probability — simulating real-world assembler filtering.
4. **Overlap Graph** — Built with a minimum overlap of 3 bp; rendered via NetworkX (skipped when > 250 nodes to avoid freezing).
