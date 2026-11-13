"""Dominant tandem-repeat period of each unique giant sequence.

For lags p = 1..2000 we measure the fraction of positions i (sampled every 7th)
with s[i] == s[i+p]. A translated tandem DNA array of unit U bp gives a protein
period of U aa (U not divisible by 3) or U/3 aa, with identity near 1. Natural
repeat proteins (titin Ig/FN3, spidroins, adhesins) usually show diverged
repeats (identity well below 0.5 at their period) - spidroins are an exception
and are flagged by the same rule, so this is a flag, not a verdict.
Output: data/period_flags.tsv
"""
import os
import sys
from multiprocessing import Pool

sys.path.insert(0, os.path.dirname(__file__))
from period import periods, read_fasta  # noqa: E402

ROOT = os.path.join(os.path.dirname(__file__), "..")
D = os.path.join(ROOT, "data")


def work(item):
    name, s = item
    top = periods(s, maxp=2000, top=3)
    return name, len(s), top


if __name__ == "__main__":
    items = [(n, s) for n, _, s in read_fasta(os.path.join(D, "fasta", "unique_ge35992.fasta"))]
    with Pool(8) as pool:
        res = pool.map(work, items, chunksize=4)
    with open(os.path.join(D, "period_flags.tsv"), "w") as f:
        f.write("md5\tlength\ttop_period_aa\ttop_period_identity\tsecond_period_aa\tsecond_identity\n")
        for name, n, top in res:
            f.write(f"{name}\t{n}\t{top[0][1]}\t{top[0][0]:.3f}\t{top[1][1]}\t{top[1][0]:.3f}\n")
    print("done", len(res))
