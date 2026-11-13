"""Collect the control results (controls/dna/results/<version>/ctl_*.json) into controls/dna/summary.json
and print one line per control: final verdict, which rules fired, and whether that is the right answer.
Usage: python3 controls/dna/summarize_controls.py v3 [v3.1 ...]
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
mods = json.load(open(os.path.join(HERE, "refs", "modifications.json")))
summary_path = os.path.join(HERE, "summary.json")
summary = json.load(open(summary_path)) if os.path.exists(summary_path) else {}
EXPECT = {"a_deletion": "against", "b_collapse": "against", "c_expansion": "against",
          "d_rrna_correct": "supported", "e_segdup_correct": "supported"}
for ver in sys.argv[1:] or ["v3"]:
    res = {}
    d = os.path.join(HERE, "results", ver)
    for f in sorted(os.listdir(d)):
        if not f.endswith(".json"):
            continue
        j = json.load(open(os.path.join(d, f)))
        tag, mode = j["job_id"][4:].rsplit("_", 1)
        il, on = j.get("illumina", {}), j.get("ont", {})
        fv = j["final_verdict"]
        want = EXPECT[tag]
        if want == "against":
            right = fv.startswith("against")
            partial = (not right) and (fv.startswith("conflicting") or "inconclusive" in fv or "outside" in fv)
            outcome = "caught" if right else ("not caught as a break; grade left unchanged" if partial else "MISSED (passes)")
        else:
            right = fv.startswith("supported")
            outcome = ("correct" if fv.startswith("supported by the reads") else
                       "correct, repeat left unresolved (grade unchanged)" if right else
                       "FALSE ALARM (against)" if fv.startswith("against") else "not confirmed (inconclusive; grade unchanged)")
        res[f"{tag}/{mode}"] = {
            "truth": mods["regions"][tag]["truth"], "final_verdict": fv, "outcome": outcome,
            "illumina_verdict": il.get("verdict"), "ont_verdict": on.get("verdict"),
            "break_boundaries": il.get("break_boundaries"), "clip_pileups_continuity_lost": len(il.get("clip_pileups_with_continuity_lost", [])),
            "max_discordant_frac_window": il.get("max_discordant_frac_window"), "depth_ratio": il.get("depth_ratio_gene_vs_contig"),
            "min_spanning_testable": il.get("min_spanning_unique_boundaries"), "unplaceable_blocks": [
                (b["start"], b["end"], b["crossing"], b["crossing_how"]) for b in il.get("library_unplaceable_blocks", [])],
            "window_repeat_blocks": il.get("repeat_blocks"), "copy_number_blocks": il.get("copy_number_blocks"),
            "ont_blocks_spanned": on.get("short_read_repeat_blocks_spanned_by_long_reads"),
            "ont_clip_pileups_continuity_lost": len(on.get("clip_pileups_with_continuity_lost", [])) if on else None,
            "reads": (j.get("illumina_reads_processed"), j.get("ont_reads_processed"))}
        print(f"{ver:5s} {tag:17s} {mode:6s} {outcome:45s} | {fv[:150]}")
    summary[ver] = res
json.dump(summary, open(summary_path, "w"), indent=1)
