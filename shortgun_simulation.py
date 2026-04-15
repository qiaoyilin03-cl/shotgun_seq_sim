import random
import os
import pysam
import matplotlib.pyplot as plt
from overlap import build_overlap_graph, draw_graph
from alignment_visualization import (
    visualize_alignment, get_contigs, get_coverage, 
    plot_contigs, create_alignment_bam, plot_lander_waterman_sweep
)
from theoretical_stats import compute_alpha, expected_contigs, expected_contig_length

# Generate random DNA genome
def generate_genome(length):
    """
    Generates a purely random DNA sequence (the 'ground truth' genome).
    
    Args:
        length (int): Total length $G$ of the reference genome.
    """
    bases = ['A', 'C', 'G', 'T']
    return ''.join(random.choice(bases) for _ in range(length))

def generate_genome_with_repeats(length, snp_rate=0.01, repeat_sizes=None, weights=None):
    """
    Generates a bio-realistic genome with various repeat families (LINE, SINE, LTR, DNA).
    
    Args:
        length (int): Total length of the genome.
        snp_rate (float): Mutation rate applied to each repeat copy.
        repeat_sizes (dict, optional): Custom lengths for each repeat family.
        weights (list, optional): Probability weights for each repeat type.
    """
    bases = ['A', 'C', 'G', 'T']
    
    def gen_random(l):
        return ''.join(random.choice(bases) for _ in range(l))
        
    def mutate(seq, rate):
        return ''.join(
            random.choice(bases) if random.random() < rate else b
            for b in seq
        )

    # Use provided repeat sizes or defaults
    if repeat_sizes is None:
        repeat_sizes = {"LINE": 800, "SINE": 150, "LTR": 400, "DNA": 300}

    # Build repeat families (source templates)
    # Using only 1 template per family makes them highly identical across the genome
    repeat_library = {
        name: [gen_random(size) for _ in range(1)]
        for name, size in repeat_sizes.items()
    }

    if weights is None:
        weights = [
            ("LINE", 0.20),
            ("SINE", 0.13),
            ("LTR", 0.08),
            ("DNA", 0.03),
            ("TANDEM", 0.05),
            ("UNIQUE", 0.51)
        ]

    def pick_type():
        r = random.random()
        cumulative = 0
        for t, w in weights:
            cumulative += w
            if r < cumulative:
                return t
        return "UNIQUE" # fallback

    genome = []
    current_len = 0

    while current_len < length:
        t = pick_type()

        if t in repeat_library:
            family = random.choice(repeat_library[t])
            seq = mutate(family, snp_rate)
        elif t == "TANDEM":
            unit = gen_random(random.randint(1, 6))
            seq = unit * random.randint(5, 20)
        else:
            seq = gen_random(100)

        genome.append(seq)
        current_len += len(seq)

    return ''.join(genome)[:length]


# Generate shotgun reads
def generate_reads(genome, num_reads, min_len, max_len):
    """
    Simulates shotgun sequencing by generating random reads from a given genome.
    """
    reads = []
    genome_len = len(genome)
    
    for _ in range(num_reads):
        read_len = random.randint(min_len, max_len)
        start = random.randint(0, genome_len - read_len)
        read = genome[start:start + read_len]
        reads.append(read)
    
    return reads

def display_simulation_results(label, genome, reads, bam_filename, output_contig_png):
    """
    Helper to extract stats from BAM, print theoretical vs empirical, and save contig plot.
    """
    create_alignment_bam(genome, reads, bam_filename)
    reads_data = []
    samfile = pysam.AlignmentFile(bam_filename, "rb")
    for read in samfile.fetch(until_eof=True):
        reads_data.append({
            'start': read.reference_start, 
            'end': read.reference_end,
            'mapq': read.mapping_quality
        })
    samfile.close()

    if os.path.exists(bam_filename):
        os.remove(bam_filename)
        
    contigs = get_contigs(reads_data)
    coverage = get_coverage(len(genome), reads_data)
    
    dropped_count = len(reads) - len(reads_data)
    
    print(f"\n--- {label} ---")
    print(f"Number of reads sampled: {len(reads)}")
    print(f"Number of reads aligned: {len(reads_data)} ({dropped_count} dropped due to mapping ambiguity)")
    print(f"Number of contigs assembled: {len(contigs)}")
    if contigs:
        lengths = [e - s for s, e in contigs]
        print(f"Max contig length: {max(lengths)} bp")
        print(f"Average contig length: {sum(lengths)/len(lengths):.2f} bp")
    
    # Theoretical expectations
    L = max(len(r) for r in reads) if reads else 0
    G = len(genome)
    N = len(reads)
    alpha = compute_alpha(N, L, G)
    exp_contigs = expected_contigs(N, alpha)
    print(f"[Theory] α (read density) = {alpha:.4f}")
    print(f"[Theory] Expected # contigs ≈ {exp_contigs:.2f}")
    
    # Save a dedicated contig plot
    fig, ax = plt.subplots(figsize=(12, 4))
    plot_contigs(ax, contigs, show_labels=True, show_gaps=True)
    ax.set_xlim(0, len(genome))
    ax.set_title(f"Contig Assembly Results ({label})")
    ax.set_xlabel("Genome Position (bp)")
    ax.set_yticks([])
    plt.tight_layout()
    plt.savefig(output_contig_png)
    print(f"Contig-only plot saved to {output_contig_png}")

def run_alignment_comparison():
    """
    Scenario 1: Highlighting Mapping Quality (MAPQ Colors)
    - High coverage (N=500).
    - Identical repeats (snp_rate=0) and short elements (80bp).
    """
    print("\n>>> Scenario 1: Alignment Visualization (MAPQ Colors) <<<")
    G = 4000
    L_min, L_max = 50, 100
    N = 600
    
    # 1. Non-repeat genome
    genome_nr = generate_genome(G)
    reads_nr = generate_reads(genome_nr, N, L_min, L_max)
    visualize_alignment(genome_nr, reads_nr, "alignment_nr_scenario1.png")
    display_simulation_results("S1: Non-Repeat Alignment", genome_nr, reads_nr, 
                               "temp_s1_nr.bam", "contigs_s1_nr.png")
    
    # 2. Repeat genome (realistic repeats)
    repeat_sizes = {"LINE": 1000, "SINE": 300, "LTR": 500, "DNA": 300}
    genome_r = generate_genome_with_repeats(G, snp_rate=0.005, repeat_sizes=repeat_sizes)
    reads_r = generate_reads(genome_r, N, L_min, L_max)
    visualize_alignment(genome_r, reads_r, "alignment_r_scenario1.png")
    display_simulation_results("S1: Repeat Alignment (MAPQ focused)", genome_r, reads_r, 
                               "temp_s1_r.bam", "contigs_s1_r.png")

def run_assembly_comparison():
    """
    Scenario 2: Highlighting Contig Fragmentation (Gaps)
    - Intermediate coverage (N=150).
    - Long repeat elements (400-800bp).
    """
    print("\n>>> Scenario 2: Contig Structure Comparison (Fragmentation) <<<")
    G = 2000
    L_min, L_max = 50, 100
    N = 150
    
    # 1. Non-repeat genome
    genome_nr = generate_genome(G)
    reads_nr = generate_reads(genome_nr, N, L_min, L_max)
    visualize_alignment(genome_nr, reads_nr, "alignment_nr_scenario2.png")
    display_simulation_results("S2: Non-Repeat Assembly", genome_nr, reads_nr, 
                               "temp_s2_nr.bam", "contigs_s2_nr.png")
    
    # 2. Repeat genome (large elements + very high frequency)
    repeat_sizes = {"LINE": 400, "SINE": 100, "LTR": 200, "DNA": 150}
    weights_r = [
        ("LINE", 0.40),
        ("SINE", 0.30),
        ("LTR", 0.10),
        ("DNA", 0.10),
        ("TANDEM", 0.00),
        ("UNIQUE", 0.10)
    ]
    genome_r = generate_genome_with_repeats(G, snp_rate=0.005, repeat_sizes=None, weights=weights_r)
    reads_r = generate_reads(genome_r, N, L_min, L_max)
    visualize_alignment(genome_r, reads_r, "alignment_r_scenario2.png")
    display_simulation_results("S2: Repeat Assembly (Fragmentation focused)", genome_r, reads_r, 
                               "temp_s2_r.bam", "contigs_s2_r.png")

def run_lander_waterman_sweep(genome_r, output_png, label, num_trials=10):
    """
    Performs a sweep of read densities (alpha) and collects contig counts
    averaged over multiple trials for both a non-repeat genome and a 
    provided repeat-rich genome.
    """
    print(f"\n>>> Running Lander-Waterman Sweep Comparison: {label} ({num_trials} trials) <<<")
    G = len(genome_r)
    # L = 100
    L_min, L_max = 50, 100
    L_mean = (L_min + L_max) / 2
    alpha_values = [0.2, 0.5, 0.8, 1.2, 1.5, 1.8, 2.2, 2.6, 3.0, 3.25, 3.6, 4.0, 4.4, 4.8, 5.2, 5.6, 6.0]
    
    sweep_data_nr = []
    sweep_data_r = []
    
    # Non-repeat genome for baseline
    genome_nr = generate_genome(G)
    
    for alpha in alpha_values:
        N = int(alpha * G / L_mean)
        
        trial_counts_nr = []
        trial_counts_r = []
        
        for trial in range(num_trials):
            # Empirical NR
            reads_nr = generate_reads(genome_nr, N, L_min, L_max)
            bam_nr = f"temp_sweep_nr_{alpha}_{trial}.bam"
            create_alignment_bam(genome_nr, reads_nr, bam_nr)
            reads_data_nr = []
            with pysam.AlignmentFile(bam_nr, "rb") as sam:
                for r in sam.fetch(until_eof=True):
                    reads_data_nr.append({'start': r.reference_start, 'end': r.reference_end})
            if os.path.exists(bam_nr): os.remove(bam_nr)
            contigs_nr = get_contigs(reads_data_nr)
            trial_counts_nr.append(len(contigs_nr))
            
            # Empirical R
            reads_r = generate_reads(genome_r, N, L_min, L_max)
            bam_r = f"temp_sweep_r_{alpha}_{trial}.bam"
            create_alignment_bam(genome_r, reads_r, bam_r)
            reads_data_r = []
            with pysam.AlignmentFile(bam_r, "rb") as sam:
                for r in sam.fetch(until_eof=True):
                    reads_data_r.append({'start': r.reference_start, 'end': r.reference_end})
            if os.path.exists(bam_r): os.remove(bam_r)
            contigs_r = get_contigs(reads_data_r)
            trial_counts_r.append(len(contigs_r))
            
        mean_nr = sum(trial_counts_nr) / num_trials
        mean_r = sum(trial_counts_r) / num_trials
        sweep_data_nr.append((alpha, mean_nr))
        sweep_data_r.append((alpha, mean_r))
        
        print(f"Alpha {alpha:.2f}: NR_mean={mean_nr:.2f}, R_mean={mean_r:.2f}")

    plot_lander_waterman_sweep(G, L_mean, sweep_data_nr, sweep_data_r, output_png)

def main():
    """Executes the comparison scenarios and the Lander-Waterman sweeps."""
    run_alignment_comparison()
    run_assembly_comparison()
    
    # 1. Sweep for Scenario 1 (Realistic repeats)
    G_s1 = 4000
    repeat_sizes_s1 = {"LINE": 1000, "SINE": 300, "LTR": 500, "DNA": 300}
    genome_s1 = generate_genome_with_repeats(G_s1, snp_rate=0.005, repeat_sizes=repeat_sizes_s1)
    run_lander_waterman_sweep(genome_s1, "lander_waterman_sweep_s1.png", "Scenario 1 (Realistic Repeats)", num_trials=3)
    
    # 2. Sweep for Scenario 2 (Long slightly mutated repeats)
    G_s2 = 2000
    repeat_sizes_s2 = {"LINE": 400, "SINE": 100, "LTR": 200, "DNA": 150}
    weights_s2 = [
        ("LINE", 0.40), ("SINE", 0.30), ("LTR", 0.10), 
        ("DNA", 0.10), ("TANDEM", 0.00), ("UNIQUE", 0.10)
    ]
    genome_s2 = generate_genome_with_repeats(G_s2, snp_rate=0.005, repeat_sizes=repeat_sizes_s2, weights=weights_s2)
    run_lander_waterman_sweep(genome_s2, "lander_waterman_sweep_s2.png", "Scenario 2 (Long Repeats)", num_trials=5)

if __name__ == "__main__":
    main()
