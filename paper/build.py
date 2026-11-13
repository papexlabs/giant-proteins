"""Build the manuscript and the supplement: numbers -> figures -> tables -> PDF (tectonic).

One command, from anywhere:  .venv/bin/python paper/build.py
Needs the project venv (matplotlib, numpy) and tectonic (/opt/homebrew/bin/tectonic, or TECTONIC=...).
tectonic fetches LaTeX packages over the network on first use and caches them; nothing else goes online.
Sources of truth: paper/manuscript.tex, paper/supplement.tex, paper/preamble.tex, paper/refs.bib.
Outputs: paper/numbers.json, paper/numbers.tex, paper/NUMBERS_generated.tsv, paper/figures/*, paper/tables/*,
paper/supplementary/*, paper/manuscript.pdf, paper/supplement.pdf, and word counts on stdout.
"""
import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PY = sys.executable
TECTONIC = os.environ.get("TECTONIC") or shutil.which("tectonic") or "/opt/homebrew/bin/tectonic"

for step in ("make_numbers.py", "make_figures.py", "make_tables.py"):
    subprocess.run([PY, os.path.join(HERE, step)], check=True, stdout=subprocess.DEVNULL)

for doc in ("manuscript", "supplement"):
    r = subprocess.run([TECTONIC, "-X", "compile", "--keep-logs", f"{doc}.tex"], cwd=HERE, capture_output=True, text=True)
    log = r.stdout + r.stderr
    if r.returncode != 0:
        print(log[-3000:])
        sys.exit(f"tectonic failed on {doc}.tex")
    bad = [l for l in log.splitlines() if re.search(r"undefined|unknown number key|Citation .* undefined|multiply defined", l)]
    texlog = open(os.path.join(HERE, f"{doc}.log")).read() if os.path.exists(os.path.join(HERE, f"{doc}.log")) else ""
    bad += [l for l in texlog.splitlines() if re.search(r"Reference .* undefined|Citation .* undefined|unknown number key", l)]
    if bad:
        print("\n".join(sorted(set(bad))))
        sys.exit(f"{doc}: unresolved references or number keys")
    for ext in (".log", ".blg"):
        p = os.path.join(HERE, doc + ext)
        if os.path.exists(p):
            os.remove(p)


def plain_words(tex):
    """Approximate word count of LaTeX: drop comments, floats (figures, tables), commands with their
    file/label arguments; keep the text inside \\nm{...} as one token per number."""
    t = re.sub(r"(?<!\\)%.*", "", tex)
    t = re.sub(r"\\begin\{(figure|table)\}.*?\\end\{\1\}", " ", t, flags=re.S)
    t = re.sub(r"\\nm\{[^}]*\}", " 0 ", t)
    t = re.sub(r"\\(cite[pt]?|ref|label|file|input|includegraphics|bibliography|bibliographystyle)\*?(\[[^]]*\])?\{[^}]*\}", " ", t)
    t = re.sub(r"\\[a-zA-Z]+\*?", " ", t)
    t = re.sub(r"[{}\[\]$\\~&]", " ", t)
    return len(re.findall(r"[A-Za-z0-9][\w'.,%-]*", t))


tex = open(os.path.join(HERE, "manuscript.tex")).read()


def section(a, b):
    i, j = tex.index(a), tex.index(b)
    return tex[i:j]


wc = {"abstract": plain_words(section("\\begin{abstract}", "\\end{abstract}")),
      "introduction to discussion and limitations": plain_words(section("\\section{Introduction}", "\\section{Methods}")),
      "methods": plain_words(section("\\section{Methods}", "\\section*{Data and code availability}"))}
wc["main text (introduction to methods)"] = wc["introduction to discussion and limitations"] + wc["methods"]
sup = open(os.path.join(HERE, "supplement.tex")).read()
wc["supplement"] = plain_words(sup[sup.index("\\begin{document}"):])
print("PDF:", os.path.join(HERE, "manuscript.pdf"), "and", os.path.join(HERE, "supplement.pdf"))
print("words:", "; ".join(f"{k} {v}" for k, v in wc.items()))
