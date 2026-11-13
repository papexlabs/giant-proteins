"""Dominant internal repeat period of a protein sequence: for each lag p (1..maxp)
the fraction of positions i with s[i] == s[i+p]; reports the best lags."""
import sys

def read_fasta(p):
    name, hdr, buf = None, "", []
    for line in open(p):
        if line.startswith(">"):
            if name:
                yield name, hdr, "".join(buf)
            name, hdr, buf = line[1:].split()[0], line, []
        else:
            buf.append(line.strip())
    if name:
        yield name, hdr, "".join(buf)

def periods(s, maxp=2000, top=5):
    n = len(s)
    res = []
    for p in range(1, min(maxp, n // 2)):
        m = sum(1 for i in range(0, n - p, 7) if s[i] == s[i + p])  # sampled every 7th position
        res.append((m / len(range(0, n - p, 7)), p))
    res.sort(reverse=True)
    return res[:top]

if __name__ == "__main__":
    fasta, ids = sys.argv[1], set(sys.argv[2:])
    for name, hdr, s in read_fasta(fasta):
        if name in ids or (len(hdr.split()) > 1 and hdr.split()[1] in ids):
            print(hdr.split()[1] if len(hdr.split()) > 1 else name, len(s), [(p, round(f, 3)) for f, p in periods(s)])
