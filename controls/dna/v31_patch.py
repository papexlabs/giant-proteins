"""Apply rules v3.1 to cloud/dna_check.py: a copy-number test added AFTER the v3 controls missed a
collapsed and/or an expanded repeat (controls/dna/summary.json, 'v3'). No v3 threshold changes.

v3.1 copy-number test (Illumina, from the stored 1-kb window arrays, so stored results can be re-judged):
  M = median all-MAPQ depth of the windows outside repeat or unplaceable blocks
  a copy-number block is a run of >= 3 consecutive windows whose all-MAPQ depth is each > sqrt(2) M
  (closer to a doubled copy number than to one copy) or each < M / sqrt(2) (closer to half)
  any such block -> 'against: copy-number anomaly' (collapsed repeat if high, expanded or absent if low)
The log-midpoints sqrt(2) and 1/sqrt(2) separate 1 copy from 2 copies and from half a copy; the
3-window minimum ignores single-window GC dips. A real copy-number difference shorter than 3 kb is
not detected, and in a metagenome a strain-variable region can also give a low run (stated caveat).
Usage: python3 controls/dna/v31_patch.py   (idempotent)
"""
import os

P = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "cloud", "dna_check.py")
s = open(P).read()
if "v3.1" in s.split("RULES = ")[1].split("\n")[0]:
    raise SystemExit("already v3.1")
s = s.replace('RULES = "v3 2026-11-12: break test only where the library can place pairs (simulated from the assembly)"',
              'RULES = ("v3.1 2026-11-13: v3 (break test only where the library can place pairs) + copy-number test "\n'
              '         "(>= 3 consecutive 1-kb windows at > sqrt(2) or < 1/sqrt(2) of the gene\'s all-MAPQ depth)")')
old = '''    res["library_unplaceable_blocks"] = blocks[:50]'''
new = '''    # v3.1 copy-number test (added 13 Nov after the controls in controls/dna): runs of >= 3 windows whose
    # all-MAPQ depth is closer to a doubled or a halved copy number than to the gene's single-copy depth
    if kind == "illumina":
        excl = set()
        for b_ in res.get("repeat_blocks", []) + blocks:
            for k in range(max(0, (b_["start"] - lo) // WIN), min(len(wa), (b_["end"] - lo) // WIN + 1)):
                excl.add(k)
        base = [wa[k] for k in range(len(wa)) if k not in excl and wa[k] > 0]
        M = sorted(base)[len(base) // 2] if base else 0
        cn, k = [], 0
        while M and k < len(wa):
            hi_, lo_ = wa[k] > M * 2 ** 0.5, wa[k] < M / 2 ** 0.5
            if hi_ or lo_:
                j = k
                while j + 1 < len(wa) and ((wa[j + 1] > M * 2 ** 0.5) if hi_ else (wa[j + 1] < M / 2 ** 0.5)):
                    j += 1
                if j - k + 1 >= 3:
                    cn.append({"start": lo + k * WIN, "end": lo + (j + 1) * WIN - 1, "windows": j - k + 1,
                               "depth_ratio": round(sum(wa[k:j + 1]) / (j - k + 1) / M, 2),
                               "kind": "excess copies (collapsed repeat?)" if hi_ else "deficit (expanded or absent?)"})
                k = j + 1
            else:
                k += 1
        res["copy_number_reference_depth"] = M
        res["copy_number_blocks"] = cn
        # the depth ratio of each repeat or unplaceable block (all-MAPQ, windows overlapping it)
        for b_ in res.get("repeat_blocks", []) + blocks:
            ks = list(range(max(0, (b_["start"] - lo) // WIN), min(len(wa), (b_["end"] - lo) // WIN + 1)))
            b_["all_mapq_depth_ratio"] = round(sum(wa[x] for x in ks) / len(ks) / M, 2) if ks and M else None
    res["library_unplaceable_blocks"] = blocks[:50]'''
assert old in s
s = s.replace(old, new, 1)
old2 = '''    nc = [b for b in m.get("library_unplaceable_blocks", []) if not b.get("crossed")]'''
new2 = '''    if kind == "illumina" and adequate and m.get("copy_number_blocks"):
        why_against.append("copy-number anomaly (v3.1): " + "; ".join(
            f"{(b['end'] - b['start'] + 1) // 1000} kb at {b['depth_ratio']}x, {b['kind']}" for b in m["copy_number_blocks"]))
    nc = [b for b in m.get("library_unplaceable_blocks", []) if not b.get("crossed")]'''
assert old2 in s
s = s.replace(old2, new2, 1)
open(P, "w").write(s)
print("cloud/dna_check.py now applies rules v3.1")
