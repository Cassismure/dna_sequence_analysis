"""
Reads directly from dna_analysis.db (built by analysis.py)
"""
import sqlite3
import matplotlib.pyplot as plt
import os

os.makedirs("charts", exist_ok=True)
conn = sqlite3.connect("dna_analysis.db")


# 1. Average GC content by organism (bar chart)
def chart_avg_gc_by_organism():
    cur = conn.execute("""
        SELECT organism, AVG(gc_content) AS avg_gc
        FROM sequences
        GROUP BY organism
        ORDER BY avg_gc DESC
    """)
    rows = cur.fetchall()
    organisms = [r[0] for r in rows]
    avg_gc = [r[1] for r in rows]

    plt.figure(figsize=(7, 5))
    plt.bar(organisms, avg_gc, color="#B04C86")
    plt.xlabel("Organism")
    plt.ylabel("Average GC Content (%)")
    plt.title("Average GC Content by Organism")
    plt.ylim(0, max(avg_gc) + 10)
    plt.tight_layout()
    plt.savefig("charts/avg_gc_by_organism.png", dpi=150)
    plt.close()


# 2. Nucleotide composition per organism (stacked bar)
def chart_nucleotide_composition():
    cur = conn.execute("""
        SELECT organism,
               SUM(a_count) AS A, SUM(t_count) AS T,
               SUM(g_count) AS G, SUM(c_count) AS C
        FROM sequences
        GROUP BY organism
    """)
    rows = cur.fetchall()
    organisms = [r[0] for r in rows]
    a = [r[1] for r in rows]
    t = [r[2] for r in rows]
    g = [r[3] for r in rows]
    c = [r[4] for r in rows]

    plt.figure(figsize=(7, 5))
    plt.bar(organisms, a, label="A", color="#55A868")
    plt.bar(organisms, t, bottom=a, label="T", color="#C44E52")
    bottom_at = [a[i] + t[i] for i in range(len(a))]
    plt.bar(organisms, g, bottom=bottom_at, label="G", color="#4C72B0")
    bottom_atg = [bottom_at[i] + g[i] for i in range(len(a))]
    plt.bar(organisms, c, bottom=bottom_atg, label="C", color="#DD8452")

    plt.xlabel("Organism")
    plt.ylabel("Total Nucleotide Count")
    plt.title("Nucleotide Composition by Organism")
    plt.legend()
    plt.tight_layout()
    plt.savefig("charts/nucleotide_composition.png", dpi=150)
    plt.close()


# 3. Sequence length distribution (histogram)
def chart_length_distribution():
    cur = conn.execute("SELECT length FROM sequences")
    lengths = [r[0] for r in cur.fetchall()]

    plt.figure(figsize=(7, 5))
    plt.hist(lengths, bins=10, color="#8172B2", edgecolor="black")
    plt.xlabel("Sequence Length (bp)")
    plt.ylabel("Number of Sequences")
    plt.title("Sequence Length Distribution")
    plt.tight_layout()
    plt.savefig("charts/length_distribution.png", dpi=150)
    plt.close()


# 4. GC content vs sequence length (scatter, colored by organism)
def chart_gc_vs_length():
    cur = conn.execute("SELECT organism, length, gc_content FROM sequences")
    rows = cur.fetchall()

    organisms = sorted(set(r[0] for r in rows))
    colors = plt.cm.tab10.colors

    plt.figure(figsize=(7, 5))
    for i, org in enumerate(organisms):
        xs = [r[1] for r in rows if r[0] == org]
        ys = [r[2] for r in rows if r[0] == org]
        plt.scatter(xs, ys, label=org, color=colors[i % len(colors)], alpha=0.8)

    plt.xlabel("Sequence Length (bp)")
    plt.ylabel("GC Content (%)")
    plt.title("GC Content vs Sequence Length")
    plt.legend()
    plt.tight_layout()
    plt.savefig("charts/gc_vs_length.png", dpi=150)
    plt.close()


# 5. Melting temperature by organism (box-style comparison via scatter + mean line)
def chart_melting_temp():
    cur = conn.execute("SELECT organism, melting_temp FROM sequences ORDER BY organism")
    rows = cur.fetchall()
    organisms = sorted(set(r[0] for r in rows))

    plt.figure(figsize=(7, 5))
    data_by_org = [[r[1] for r in rows if r[0] == org] for org in organisms]
    plt.boxplot(data_by_org, tick_labels=organisms)
    plt.xlabel("Organism")
    plt.ylabel("Melting Temperature (°C)")
    plt.title("Melting Temperature Distribution by Organism")
    plt.tight_layout()
    plt.savefig("charts/melting_temp_by_organism.png", dpi=150)
    plt.close()


if __name__ == "__main__":
    chart_avg_gc_by_organism()
    chart_nucleotide_composition()
    chart_length_distribution()
    chart_gc_vs_length()
    chart_melting_temp()
    print("Saved 5 charts to charts/")

conn.close()
