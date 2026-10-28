"""
core.py - VariantLens analysis engine (merged team module).

Pipeline: validate -> compare/align -> variant events -> 3' normalize -> HGVS name
          -> consequence -> database match

Conventions shared by every function:
  * Positions are 1-based coding-DNA (c.) positions on the HBB CDS (NM_000518.5);
    position 1 is the A of the ATG start codon.
  * A variant event is a dict {"type", "start", "end", "ref", "alt"} where type is
    substitution | deletion | insertion | delins. For an insertion, start and end are
    the two reference bases on either side of the insertion point.

Authors: Member B (validation, comparison, normalization, HGVS),
         Member C (genetic code, consequence, legacy names, lookup ideas),
         Member A (data, test helpers, merge).
"""
import json
from itertools import product
from pathlib import Path
from urllib.parse import quote

from Bio.Align import PairwiseAligner

DATA_DIR = Path(__file__).resolve().parent / "data"


# =====================================================================
# Step 4: input handling and validation  (Member B)
# =====================================================================

def clean(text):
    """Return a validated DNA sequence from raw text or FASTA input.

    FASTA header lines (first non-space character ``>``), blank lines, spaces and
    line-wrapping are ignored. Only A/C/G/T are accepted; ambiguous bases such as
    N are rejected deliberately.
    """
    if not isinstance(text, str):
        raise TypeError("Sequence input must be a string.")
    sequence_lines = []
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith(">"):
            continue
        sequence_lines.append("".join(line.split()))
    seq = "".join(sequence_lines).upper()
    if not seq:
        raise ValueError("Sequence is empty.")
    bad = sorted(set(seq) - set("ACGT"))
    if bad:
        raise ValueError(f"Invalid characters found: {bad}. Only A, C, G and T are allowed.")
    return seq


def length_warning(ref, test, limit=0.20):
    """Warning text if the test length differs from the reference by more than 20%."""
    diff = abs(len(test) - len(ref)) / len(ref)
    if diff > limit:
        return (f"Test sequence is {len(test)} bases but the reference is {len(ref)} "
                f"({diff:.0%} different). Is this the right gene?")
    return None


# =====================================================================
# Step 5: comparison engine  (Member B)
# =====================================================================

def make_aligner():
    """Global pairwise aligner used to compare two DNA sequences."""
    aligner = PairwiseAligner()
    aligner.mode = "global"          # compare the whole length of both sequences
    aligner.match_score = 2
    aligner.mismatch_score = -3
    aligner.open_gap_score = -5      # opening a gap is costly...
    aligner.extend_gap_score = -1    # ...extending it is cheap, so indels stay in one block
    return aligner


def _gap_event(events, ref, test, r0, r1, t0, t1):
    """Convert one unaligned block (0-based, half-open) to an event."""
    deleted, inserted = ref[r0:r1], test[t0:t1]
    if deleted and inserted:
        events.append({"type": "delins", "start": r0 + 1, "end": r1, "ref": deleted, "alt": inserted})
    elif deleted:
        events.append({"type": "deletion", "start": r0 + 1, "end": r1, "ref": deleted, "alt": ""})
    elif inserted:
        events.append({"type": "insertion", "start": r0, "end": r0 + 1, "ref": "", "alt": inserted})


def call_variants(ref, test, normalize=True):
    """Return sequence differences as event dicts with 1-based positions.

    Equal-length sequences are compared base by base; sequences of different lengths
    are globally aligned so insertions and deletions are found. Indels are shifted to
    their 3'-most equivalent position (HGVS rule) unless normalize=False.
    Inputs should already be cleaned with clean().
    """
    if not isinstance(ref, str) or not isinstance(test, str):
        raise TypeError("Reference and test sequences must be strings.")
    if not ref or not test:
        raise ValueError("Reference and test sequences must not be empty.")
    if set(ref) - set("ACGT") or set(test) - set("ACGT"):
        raise ValueError("Sequences must contain only A/C/G/T; call clean() first.")

    events = []
    if len(ref) == len(test):
        for i, (rbase, tbase) in enumerate(zip(ref, test)):
            if rbase != tbase:
                events.append({"type": "substitution", "start": i + 1, "end": i + 1,
                               "ref": rbase, "alt": tbase})
        return events

    alignment = make_aligner().align(ref, test)[0]
    ref_pos = test_pos = 0
    ref_blocks, test_blocks = alignment.aligned
    for (rs, re), (ts, te) in zip(ref_blocks, test_blocks):
        rs, re, ts, te = int(rs), int(re), int(ts), int(te)
        _gap_event(events, ref, test, ref_pos, rs, test_pos, ts)
        for offset in range(min(re - rs, te - ts)):      # aligned blocks can still hold mismatches
            rbase, tbase = ref[rs + offset], test[ts + offset]
            if rbase != tbase:
                pos = rs + offset + 1
                events.append({"type": "substitution", "start": pos, "end": pos,
                               "ref": rbase, "alt": tbase})
        ref_pos, test_pos = re, te
    _gap_event(events, ref, test, ref_pos, len(ref), test_pos, len(test))
    if normalize:
        events = [normalize_indel(e, ref) for e in events]
    return sorted(events, key=lambda e: (e["start"], e["end"], e["type"]))


# =====================================================================
# Step 6: 3' normalization and HGVS naming  (Member B; insertion fix by Member A)
# =====================================================================

def normalize_indel(e, ref):
    """Return the 3'-most equivalent insertion/deletion event (HGVS 3' rule).

    In a repeat such as TTCT, deleting or inserting at several positions gives the
    same sequence; HGVS reports the one furthest towards the 3' end.
    Substitutions and delins are returned unchanged.
    """
    event = dict(e)
    if event["type"] == "deletion":
        start, end = event["start"], event["end"]
        while end < len(ref) and ref[start - 1] == ref[end]:
            start += 1
            end += 1
        event.update(start=start, end=end, ref=ref[start - 1:end])
    elif event["type"] == "insertion":
        pos, inserted = event["start"], event["alt"]
        # Shift right while the next reference base equals the first inserted base.
        # The inserted bases rotate as they move (e.g. "TC" before T C... becomes "CT").
        while pos < len(ref) and inserted and ref[pos] == inserted[0]:
            inserted = inserted[1:] + inserted[0]
            pos += 1
        event.update(start=pos, end=pos + 1, alt=inserted)
    return event


def to_hgvs(e, ref):
    """Format an event as a coding-DNA HGVS expression, e.g. c.20A>T, c.27dup."""
    event = normalize_indel(e, ref)
    kind, start, end = event["type"], event["start"], event["end"]
    span = str(start) if start == end else f"{start}_{end}"
    if kind == "substitution":
        return f"c.{start}{event['ref']}>{event['alt']}"
    if kind == "deletion":
        return f"c.{span}del"
    if kind == "insertion":
        inserted, n = event["alt"], len(event["alt"])
        # Duplication takes precedence when the insertion repeats the bases just before it.
        if start >= n and ref[start - n:start] == inserted:
            return f"c.{start - n + 1}dup" if n == 1 else f"c.{start - n + 1}_{start}dup"
        return f"c.{start}_{end}ins{inserted}"
    return f"c.{span}delins{event['alt']}"


# =====================================================================
# Step 7: consequence prediction  (Member C)
# =====================================================================

_AAS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON = {"".join(c): a for c, a in zip(product("TCAG", repeat=3), _AAS)}   # standard genetic code

THREE = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys", "Q": "Gln",
    "E": "Glu", "G": "Gly", "H": "His", "I": "Ile", "L": "Leu", "K": "Lys",
    "M": "Met", "F": "Phe", "P": "Pro", "S": "Ser", "T": "Thr", "W": "Trp",
    "Y": "Tyr", "V": "Val", "*": "Ter",
}


def translate(seq):
    """Translate until the first stop (stop included). Ignores a trailing partial codon."""
    prot = []
    for i in range(0, len(seq) - 2, 3):
        aa = CODON.get(seq[i:i + 3], "X")
        prot.append(aa)
        if aa == "*":
            break
    return "".join(prot)


def apply_event(ref, e):
    """Return the mutated sequence that an event describes."""
    if e["type"] == "insertion":
        return ref[:e["start"]] + e["alt"] + ref[e["start"]:]
    return ref[:e["start"] - 1] + e["alt"] + ref[e["end"]:]


def annotate(ref, e):
    """Classify one variant's effect on the protein.

    Returns {"effect", "hgvs_p", "aa_pos", "legacy_name", "codon"}.
    effect is one of: synonymous, missense, nonsense, start-loss, stop-loss,
    frameshift, in-frame deletion, in-frame insertion, in-frame delins.
    legacy_name uses haemoglobin numbering (initiator Met not counted),
    e.g. Glu6Val for HbS, because that is what the clinical literature uses.
    """
    mut = apply_event(ref, e)
    ref_p, alt_p = translate(ref), translate(mut)
    codon = ""
    if e["type"] == "substitution":
        i, off = (e["start"] - 1) // 3, (e["start"] - 1) % 3
        old = ref[i * 3:i * 3 + 3]
        codon = f"{old} > {old[:off] + e['alt'] + old[off + 1:]}"

    if ref_p == alt_p:
        if e["type"] == "substitution":
            n = (e["start"] - 1) // 3 + 1
            return {"effect": "synonymous", "hgvs_p": f"p.{THREE[ref_p[n - 1]]}{n}=",
                    "aa_pos": n, "legacy_name": None, "codon": codon}
        return {"effect": "synonymous", "hgvs_p": "p.=", "aa_pos": None,
                "legacy_name": None, "codon": codon}

    # first amino acid that differs
    i = next((k for k in range(min(len(ref_p), len(alt_p))) if ref_p[k] != alt_p[k]),
             min(len(ref_p), len(alt_p)))
    n = i + 1
    r = ref_p[i] if i < len(ref_p) else None
    a = alt_p[i] if i < len(alt_p) else None
    r3 = THREE.get(r, "?") if r else "?"
    a3 = THREE.get(a, "?") if a else "?"

    def result(effect, hgvs_p, legacy=None):
        return {"effect": effect, "hgvs_p": hgvs_p, "aa_pos": n,
                "legacy_name": legacy, "codon": codon}

    shift = (len(mut) - len(ref)) % 3 != 0
    if i == 0:
        return result("start-loss", "p.Met1?")
    if shift:
        return result("frameshift", f"p.{r3}{n}{a3}fs")
    if r == "*":
        return result("stop-loss", f"p.Ter{n}{a3}ext")
    if e["type"] == "substitution":
        if a == "*":
            return result("nonsense", f"p.{r3}{n}Ter", f"{r3}{n - 1}Ter")
        return result("missense", f"p.{r3}{n}{a3}", f"{r3}{n - 1}{a3}")
    if e["type"] == "deletion":
        k = len(e["ref"]) // 3
        end = n + k - 1
        txt = f"p.{r3}{n}del" if k == 1 else f"p.{r3}{n}_{THREE.get(ref_p[end - 1], '?')}{end}del"
        return result("in-frame deletion", txt)
    if e["type"] == "insertion":
        return result("in-frame insertion", f"p.{r3}{n}ins")       # simplified naming
    return result("in-frame delins", "p.?")


def consequence(ref, e):
    """Return (effect, protein_change), e.g. ("missense", "p.Glu7Val")."""
    res = annotate(ref, e)
    return res["effect"], res["hgvs_p"]


# =====================================================================
# Step 8: database match and links  (Members C and A)
# =====================================================================

def load_db(path=DATA_DIR / "hbb_variants.json"):
    """Load the known-variant table (a list of entries)."""
    with open(path, encoding="utf-8") as f:
        data = json.load(f)
    if isinstance(data, dict):                    # tolerate {"variants": [...]}
        data = data.get("variants", list(data.values()))
    return data


def load_reference(path=DATA_DIR / "hbb_cds.txt"):
    with open(path, encoding="utf-8") as f:
        return clean(f.read())


def _norm(hgvs):
    """'NM_000518.5:c.20A>T ' -> 'c.20A>T'"""
    hgvs = hgvs.strip()
    return hgvs.split(":", 1)[1] if ":" in hgvs else hgvs


def lookup(hgvs, db):
    """Return the database entry for an HGVS c. name (or one of its aliases), else None."""
    hgvs = _norm(hgvs)
    for entry in db:
        if hgvs == entry["hgvs"] or hgvs in entry.get("aliases", []):
            return entry
    return None                                   # not in database: uncharacterized, NOT benign


def clinvar_link(hgvs, entry=None):
    """Direct ClinVar record if verified in the database, otherwise a ClinVar search link."""
    if entry and entry.get("clinvar_url"):
        return entry["clinvar_url"]
    return "https://www.ncbi.nlm.nih.gov/clinvar/?term=" + quote(f"HBB[gene] AND {_norm(hgvs)}")


# =====================================================================
# Step 9b: test-sequence generators  (used by the app examples and tests)
# =====================================================================

def substitute(ref, pos, expected, alt):
    """Change the base at 1-based pos to alt; `expected` guards against wrong numbering."""
    assert ref[pos - 1] == expected, f"Position {pos} is {ref[pos - 1]}, not {expected}"
    return ref[:pos - 1] + alt + ref[pos:]


def delete(ref, start, end):
    """Delete bases start..end (1-based, inclusive)."""
    return ref[:start - 1] + ref[end:]


def insert_after(ref, pos, bases):
    """Insert bases after 1-based pos (pos=0 inserts at the very start)."""
    return ref[:pos] + bases + ref[pos:]


# =====================================================================
# Full analysis
# =====================================================================

def analyze(ref, test, db):
    """Run the whole pipeline. Returns (events, rows) with one row per variant."""
    events = call_variants(ref, test)
    rows = []
    for e in events:
        name = to_hgvs(e, ref)
        info = annotate(ref, e)
        hit = lookup(name, db)
        rows.append({
            "event": e, "hgvs": name, "effect": info["effect"], "protein": info["hgvs_p"],
            "codon": info["codon"], "legacy": hit["legacy"] if hit else info["legacy_name"],
            "hit": hit, "clinvar": clinvar_link(name, hit),
        })
    return events, rows
