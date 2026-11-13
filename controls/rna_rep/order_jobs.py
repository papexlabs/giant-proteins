#!/usr/bin/env python3
"""Order cloud/bundle_rep_rna/jobs.tsv so that the species whose models the RNA test can decide
run first (highest share of introns with a unique junction, controls/rna_rep/junction_mappability.json),
in case the box stops early. Usage: python3 controls/rna_rep/order_jobs.py"""
import json
import os

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")
J = os.path.join(ROOT, "cloud", "bundle_rep_rna", "jobs.tsv")
mp = json.load(open(os.path.join(ROOT, "controls", "rna_rep", "junction_mappability.json")))
best = {}
for x in mp:
    if x.get("n_introns"):
        best[x["species_key"]] = max(best.get(x["species_key"], 0), x["share_unique_junction"])
lines = open(J).read().splitlines()
head, rows = lines[0], lines[1:]
rows.sort(key=lambda l: -best.get(l.split("\t")[0], 0))
open(J, "w").write("\n".join([head] + rows) + "\n")
for l in rows:
    print(l.split("\t")[0], best.get(l.split("\t")[0]))
