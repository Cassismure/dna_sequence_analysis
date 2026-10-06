import csv
import sqlite3
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
from Bio.SeqUtils import MeltingTemp as mt


def analyze_sequence(sequence: str) -> dict:
    seq_obj = Seq(sequence.upper())
    gc_pct = round(gc_fraction(seq_obj) * 100, 2)

    # Tm_Wallace (2*(A+T) + 4*(G+C)) only works for short primers (<14bp) —
    # Tm_GC is the correct empirical formula for longer sequences,
    tm = round(mt.Tm_GC(seq_obj), 2)

    counts = {base: sequence.upper().count(base) for base in "ATGC"}

    return {
        "length": len(sequence),
        "gc_content": gc_pct,
        "at_content": round(100 - gc_pct, 2),
        "melting_temp": tm,
        "a_count": counts["A"],
        "t_count": counts["T"],
        "g_count": counts["G"],
        "c_count": counts["C"],
    }


def load_sequences(csv_path: str) -> list[dict]:
    with open(csv_path, newline="") as f:
        return list(csv.DictReader(f))


def analyze(sequences: list[dict]) -> list[dict]:
    results = []
    for row in sequences:
        metrics = analyze_sequence(row["sequence"])
        results.append({
            "sequence_id": row["sequence_id"],
            "organism": row["organism"],
            "sequence": row["sequence"],
            **metrics,
        })
    return results


def build_database(results: list[dict], db_path: str = "dna_analysis.db"):
    conn = sqlite3.connect(db_path)
    cur = conn.cursor()

    cur.executescript("""
    DROP TABLE IF EXISTS sequences;
    CREATE TABLE sequences (
        sequence_id   TEXT PRIMARY KEY,
        organism      TEXT NOT NULL,
        sequence      TEXT NOT NULL,
        length        INTEGER NOT NULL,
        gc_content    REAL NOT NULL,
        at_content    REAL NOT NULL,
        melting_temp  REAL NOT NULL,
        a_count       INTEGER NOT NULL,
        t_count       INTEGER NOT NULL,
        g_count       INTEGER NOT NULL,
        c_count       INTEGER NOT NULL
    );
    """)

    cur.executemany("""
        INSERT INTO sequences
        (sequence_id, organism, sequence, length, gc_content, at_content,
         melting_temp, a_count, t_count, g_count, c_count)
        VALUES (:sequence_id, :organism, :sequence, :length, :gc_content,
                :at_content, :melting_temp, :a_count, :t_count, :g_count, :c_count)
    """, results)

    conn.commit()
    conn.close()
    print(f"Loaded {len(results)} sequences into {db_path}")


if __name__ == "__main__":
    raw = load_sequences("data/sequences.csv")
    results = analyze(raw)
    build_database(results)

    for r in results[:3]:
        print(r["sequence_id"], r["organism"], f"GC={r['gc_content']}%",
              f"Tm={r['melting_temp']}°C", f"len={r['length']}")
