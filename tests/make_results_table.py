# tests/make_results_table.py
# Owner: Member A. Runs every Step 10 test case through the real pipeline and writes
# results_table.md for the report.   Run from the project root:  python tests/make_results_table.py
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from core import *  # noqa: E402,F403

REF = clean((ROOT / "data" / "hbb_cds.txt").read_text())
DB = load_db(str(ROOT / "data" / "hbb_variants.json"))

CASES = [
    ("Reference vs itself", "No variants", lambda: REF),
    ("HbS c.20A>T", "1 substitution, missense, p.Glu7Val, HbS", lambda: substitute(REF, 20, "A", "T")),
    ("HbC c.19G>A", "1 substitution, missense, p.Glu7Lys, HbC", lambda: substitute(REF, 19, "G", "A")),
    ("HbE c.79G>A", "1 substitution, missense, p.Glu27Lys, HbE", lambda: substitute(REF, 79, "G", "A")),
    ("Codon 39 c.118C>T", "nonsense, p.Gln40Ter", lambda: substitute(REF, 118, "C", "T")),
    ("Codon 15 c.47G>A", "nonsense, p.Trp16Ter", lambda: substitute(REF, 47, "G", "A")),
    ("Delete bases 126-129", "c.126_129del, frameshift, Codon 41/42", lambda: delete(REF, 126, 129)),
    ("Duplicate base 27", "c.27dup, frameshift, Codon 8/9", lambda: insert_after(REF, 27, REF[26])),
    ("Synonymous c.21G>A", "synonymous, p.Glu7=", lambda: substitute(REF, 21, "G", "A")),
    ("Not in table c.4G>A", "missense, Not in local database", lambda: substitute(REF, 4, "G", "A")),
    ("Two substitutions (c.20 + c.79)", "2 separate events", lambda: substitute(substitute(REF, 20, "A", "T"), 79, "G", "A")),
    ("Start codon c.3G>A", "start-loss, p.Met1?", lambda: substitute(REF, 3, "G", "A")),
    ("In-frame deletion 100-102", "in-frame deletion", lambda: delete(REF, 100, 102)),
    ("Invalid characters 'ACGTX'", "Clear error, no crash", "ACGTX"),
    ("Empty input", "Clear error, no crash", ""),
]


def run(case):
    label, expected, build = case
    if isinstance(build, str):
        try:
            clean(build)
            return label, expected, "Accepted (unexpected)", "FAIL"
        except ValueError as err:
            return label, expected, f"Error: {err}", "PASS"
    events = call_variants(REF, build())
    if not events:
        return label, expected, "No variants", "PASS" if expected == "No variants" else "FAIL"
    parts = []
    for e in events:
        name = to_hgvs(e, REF)
        effect, protein = consequence(REF, e)
        hit = lookup(name, DB)
        parts.append(f"{name}, {effect}, {protein}, " + (hit["name"] if hit else "Not in local database"))
    got = "; ".join(parts)
    # simple automatic check: every expected keyword appears in the output
    keys = [k.strip() for k in expected.replace(",", " ").split() if k.startswith(("p.", "c.", "Codon", "Hb", "missense", "start-loss",
            "nonsense", "frameshift", "synonymous", "in-frame", "Not"))]
    ok = all(k in got for k in keys) and (("2 separate" not in expected) or len(events) == 2)
    return label, expected, got, "PASS" if ok else "CHECK"


if __name__ == "__main__":
    lines = ["| # | Test case | Expected | Actual output | Result |", "|---|---|---|---|---|"]
    for i, c in enumerate(CASES, 1):
        label, exp, got, res = run(c)
        lines.append(f"| {i} | {label} | {exp} | {got.replace('|', '/')} | {res} |")
    out = "\n".join(lines) + "\n"
    (ROOT / "results_table.md").write_text(out)
    print(out)
