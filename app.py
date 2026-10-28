"""VariantLens: Streamlit interface (Member D, merged with the team's core.py).

Run:  streamlit run app.py
Needs core.py and the data/ folder (hbb_cds.txt, hbb_variants.json) next to this file.
"""
from html import escape
from pathlib import Path

import pandas as pd
import streamlit as st

from core import (analyze, apply_event, clean, delete, insert_after, length_warning, load_db,
                  load_reference, substitute, translate)

HERE = Path(__file__).parent
# Development fallback only (first 14 codons of HBB). The real hbb_cds.txt replaces it.
DEV_REF = "ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCC"

COLORS = {
    "substitution": "#e5383b",  # red
    "deletion": "#7b2cbf",      # purple
    "insertion": "#f77f00",     # orange
    "delins": "#2a9d8f",        # teal
}

EFFECT_HELP = {
    "synonymous": "Codon changes but the amino acid stays the same.",
    "missense": "One amino acid is replaced by another.",
    "nonsense": "Codon becomes a stop; the protein is cut short.",
    "start-loss": "The ATG start codon is lost; the protein may not be made.",
    "stop-loss": "The stop codon is lost; the protein runs on.",
    "frameshift": "Indel not a multiple of 3; every codon after it is misread.",
    "in-frame deletion / insertion": "Multiple of 3 bases; amino acids removed or added.",
}


# ------------------------------------------------------------------ highlighting
def highlight_html(ref, events, line_width=60):
    """HTML for the reference sequence with variant positions coloured.

    Event format (1-based): {"type", "start", "end", "ref", "alt"}.
    For an insertion, start and end are the reference bases on either side of
    the insertion point; both flank bases are highlighted.
    """
    marks = {}
    for e in events:
        for pos in range(e["start"], e["end"] + 1):
            if 1 <= pos <= len(ref):
                marks[pos] = e["type"]

    lines = []
    for i in range(0, len(ref), line_width):
        chunk = ref[i:i + line_width]
        parts = []
        for j, base in enumerate(chunk):
            pos = i + j + 1
            kind = marks.get(pos)
            if kind:
                colour = COLORS.get(kind, "#555555")
                parts.append(
                    f'<span style="background:{colour};color:#fff;border-radius:3px;'
                    f'padding:0 1px;font-weight:700" title="{kind} at position {pos}">{base}</span>'
                )
            else:
                parts.append(base)
            if (j + 1) % 10 == 0 and j + 1 < len(chunk):
                parts.append(" ")
        label = (
            '<span style="color:#888;display:inline-block;width:3.5em;'
            f'text-align:right;margin-right:1em">{i + 1}</span>'
        )
        lines.append(f"<div>{label}{''.join(parts)}</div>")

    legend = " ".join(
        f'<span style="background:{c};color:#fff;border-radius:3px;padding:0 6px">{k}</span>'
        for k, c in COLORS.items()
    )
    return (
        '<div style="font-family:Consolas,Menlo,monospace;font-size:15px;'
        'line-height:1.8;overflow-x:auto;white-space:pre">'
        + "".join(lines)
        + f'</div><div style="margin-top:8px;font-size:13px">{legend}</div>'
    )


# ------------------------------------------------------------------ data loading
@st.cache_resource
def read_reference():
    """Return (sequence, loaded_from_file)."""
    try:
        return load_reference(), True
    except (OSError, ValueError):
        return DEV_REF, False


@st.cache_resource
def get_db():
    try:
        return load_db()
    except (OSError, ValueError):
        return []


def make_examples(ref):
    """Example test sequences built from the reference (needs the full 444-base CDS)."""
    if len(ref) < 444:
        return {"Normal (no change)": ref}
    return {
        "Normal (no change)": ref,
        "HbS (sickle cell): c.20A>T": substitute(ref, 20, "A", "T"),
        "HbC: c.19G>A": substitute(ref, 19, "G", "A"),
        "HbE: c.79G>A": substitute(ref, 79, "G", "A"),
        "Codon 39 nonsense: c.118C>T": substitute(ref, 118, "C", "T"),
        "Codon 15 nonsense: c.47G>A": substitute(ref, 47, "G", "A"),
        "Codon 41/42 deletion: c.126_129del": delete(ref, 126, 129),
        "Codon 8/9 insertion: c.27dup": insert_after(ref, 27, ref[26]),
        "Synonymous: c.21G>A": substitute(ref, 21, "G", "A"),
        "Start-loss: c.3G>A": substitute(ref, 3, "G", "A"),
        "Unknown variant: c.4G>A": substitute(ref, 4, "G", "A"),
        "Two mutations: HbS + HbE": substitute(substitute(ref, 20, "A", "T"), 79, "G", "A"),
        "Invalid input (contains X)": ref[:30] + "X" + ref[31:],
    }


# ------------------------------------------------------------------ result rows
def change_text(e):
    if e["type"] == "substitution":
        return f"{e['ref']} → {e['alt']}"
    if e["type"] == "deletion":
        return f"del {e['ref']}"
    if e["type"] == "insertion":
        return f"ins {e['alt']}"
    return f"{e['ref']} → {e['alt']}"


# ------------------------------------------------------------------ plain-language conclusions
AA_NAMES = {
    "Ala": "Alanine", "Arg": "Arginine", "Asn": "Asparagine", "Asp": "Aspartic acid",
    "Cys": "Cysteine", "Gln": "Glutamine", "Glu": "Glutamic acid", "Gly": "Glycine",
    "His": "Histidine", "Ile": "Isoleucine", "Leu": "Leucine", "Lys": "Lysine",
    "Met": "Methionine", "Phe": "Phenylalanine", "Pro": "Proline", "Ser": "Serine",
    "Thr": "Threonine", "Trp": "Tryptophan", "Tyr": "Tyrosine", "Val": "Valine",
}

# predicted impact on the protein (not a clinical verdict)
IMPACT = {
    "nonsense": "high", "frameshift": "high", "start-loss": "high", "stop-loss": "high",
    "missense": "moderate", "in-frame deletion": "moderate", "in-frame insertion": "moderate",
    "in-frame delins": "moderate", "synonymous": "low",
}
IMPACT_STYLE = {
    "high": ("#d62828", "High predicted impact"),
    "moderate": ("#e09f3e", "Moderate predicted impact"),
    "low": ("#2a9d8f", "Low predicted impact"),
}
UNKNOWN = "Not in database (uncharacterized)"


def aa(code):
    """'Glu' -> 'Glutamic acid (Glu)'."""
    return f"{AA_NAMES[code]} ({code})" if code in AA_NAMES else code


def dna_text(e):
    if e["type"] == "substitution":
        return f"{e['ref']} → {e['alt']} at position {e['start']}"
    if e["type"] == "deletion":
        n = len(e["ref"])
        where = f"position {e['start']}" if n == 1 else f"positions {e['start']}–{e['end']}"
        return f"{n} base{'s' if n > 1 else ''} deleted ({e['ref']}) at {where}"
    if e["type"] == "insertion":
        n = len(e["alt"])
        return f"{n} base{'s' if n > 1 else ''} inserted ({e['alt']}) after position {e['start']}"
    return f"{e['ref']} replaced by {e['alt']} at positions {e['start']}–{e['end']}"


def protein_story(ref, r):
    """One plain-English sentence on what the change does to the beta-globin protein."""
    e, eff, prot = r["event"], r["effect"], r["protein"]
    normal_len = len(translate(ref).rstrip("*"))
    new_len = len(translate(apply_event(ref, e)).rstrip("*"))
    p = prot[2:] if prot.startswith("p.") else prot
    if eff == "synonymous":
        return "The codon changes but still codes for the same amino acid, so the protein is unchanged."
    if eff == "missense":
        old, pos, new = p[:3], p[3:-3], p[-3:]
        return f"Amino acid {pos} changes from {aa(old)} to {aa(new)}. One building block of the protein is swapped."
    if eff == "nonsense":
        return (f"Amino acid {p[3:-3]}, {aa(p[:3])}, becomes a STOP signal. The protein is cut short "
                f"to {new_len} of its normal {normal_len} amino acids.")
    if eff == "frameshift":
        pos = "".join(ch for ch in p[3:] if ch.isdigit())
        return (f"The reading frame shifts at amino acid {pos}, "
                f"so every codon after it is misread. The altered protein is {new_len} amino acids "
                f"long instead of {normal_len}.")
    if eff == "start-loss":
        return "The ATG start codon is destroyed, so the protein may not be made at all."
    if eff == "stop-loss":
        return f"The stop codon is lost, so the protein runs on past its normal end ({normal_len} amino acids)."
    if eff == "in-frame deletion":
        k = len(e["ref"]) // 3
        return f"{k} amino acid{'s are' if k > 1 else ' is'} removed; the rest of the protein is read correctly."
    if eff == "in-frame insertion":
        k = len(e["alt"]) // 3
        return f"{k} amino acid{'s are' if k > 1 else ' is'} added; the rest of the protein is read correctly."
    return "Several bases are replaced; the protein changes at this point."


def sentence_case(cond):
    """'Sickle cell disease' -> 'sickle cell disease', but keep 'HbC disease'."""
    return cond if cond.startswith("Hb") else cond[0].lower() + cond[1:]


def overall_conclusion(rows):
    """Headline sentence and colour for the whole analysis."""
    n = len(rows)
    known = [r for r in rows if r["hit"]]
    unknown = [r for r in rows if not r["hit"]]
    levels = [IMPACT.get(r["effect"], "moderate") for r in rows]
    top = "high" if "high" in levels else "moderate" if "moderate" in levels else "low"
    parts = [f"<b>{n} variant{'s' if n > 1 else ''} found.</b>"]
    for r in known:
        cond = sentence_case(r["hit"]["condition"].split(";")[0].strip())
        parts.append(f"{escape(r['hgvs'])} is the known variant <b>{escape(r['hit']['name'])}</b>, "
                     f"associated with {escape(cond)}.")
    if unknown:
        names = ", ".join(escape(r["hgvs"]) for r in unknown)
        verb = "is" if len(unknown) == 1 else "are"
        parts.append(f"{names} {verb} not in the database, so {'it is' if len(unknown) == 1 else 'they are'} "
                     "<b>uncharacterized</b>.")
    return " ".join(parts), top, len(known), len(unknown)


def conclusion_banner_html(rows):
    text, top, n_known, n_unknown = overall_conclusion(rows)
    colour, label = IMPACT_STYLE[top]
    return f"""
<div style="border:2px solid {colour};border-radius:12px;padding:18px 22px;margin:8px 0 18px;
            background:{colour}14">
  <div style="font-size:13px;letter-spacing:.06em;text-transform:uppercase;opacity:.75">Conclusion</div>
  <div style="font-size:20px;line-height:1.5;margin:6px 0 12px">{text}</div>
  <span style="background:{colour};color:#fff;border-radius:20px;padding:3px 12px;font-size:13px;
               font-weight:600">Highest: {label.lower()}</span>
  <span style="border:1px solid currentColor;border-radius:20px;padding:2px 12px;font-size:13px;
               margin-left:6px;opacity:.8">{n_known} known · {n_unknown} uncharacterized</span>
  <div style="font-size:12.5px;opacity:.7;margin-top:12px">Predicted effect on the beta-globin
    protein only. Educational result, not a diagnosis.</div>
</div>"""


def variant_card_html(ref, r, i):
    e, hit = r["event"], r["hit"]
    colour, label = IMPACT_STYLE[IMPACT.get(r["effect"], "moderate")]
    if hit:
        status = (f'<span style="background:#1d3557;color:#fff;border-radius:20px;padding:3px 12px;'
                  f'font-size:13px;font-weight:600">Known: {escape(hit["name"])}</span>')
        known_line = (f"<b>Known variant.</b> Associated with {escape(sentence_case(hit['condition']))}. "
                      f"Legacy name: {escape(hit['legacy'])}.")
    else:
        status = ('<span style="border:1px solid currentColor;border-radius:20px;padding:2px 12px;'
                  'font-size:13px;opacity:.8">Uncharacterized</span>')
        known_line = ("<b>Not in the local database.</b> This means it is uncharacterized here, "
                      "not that it is harmless. Check ClinVar.")
    codon = escape(r["codon"]) if r["codon"] else "—"
    step = ('<div style="flex:1 1 180px;border-radius:8px;padding:10px 12px;'
            'background:rgba(128,128,128,.10)">'
            '<div style="font-size:11.5px;text-transform:uppercase;letter-spacing:.05em;opacity:.7">{t}</div>'
            '<div style="font-size:15px;font-weight:600;margin-top:2px">{v}</div></div>')
    arrow = '<div style="align-self:center;font-size:20px;opacity:.5">→</div>'
    steps = arrow.join([
        step.format(t="1 · DNA change", v=escape(dna_text(e))),
        step.format(t="2 · Codon", v=codon),
        step.format(t="3 · Protein", v=escape(r["protein"])),
        step.format(t="4 · Effect", v=escape(r["effect"].capitalize())),
    ])
    return f"""
<div style="border:1px solid rgba(128,128,128,.35);border-left:6px solid {colour};border-radius:10px;
            padding:16px 20px;margin-bottom:14px">
  <div style="display:flex;flex-wrap:wrap;gap:8px;align-items:center;margin-bottom:10px">
    <span style="font-size:13px;opacity:.6">Variant {i}</span>
    <span style="font-family:Consolas,Menlo,monospace;font-size:20px;font-weight:700">{escape(r['hgvs'])}</span>
    {status}
    <span style="background:{colour};color:#fff;border-radius:20px;padding:3px 12px;font-size:13px;
                 font-weight:600">{label}</span>
  </div>
  <div style="display:flex;flex-wrap:wrap;gap:8px;margin-bottom:12px">{steps}</div>
  <div style="font-size:15.5px;line-height:1.55;margin-bottom:6px">{escape(protein_story(ref, r))}</div>
  <div style="font-size:15px;line-height:1.55">{known_line}
    <a href="{r['clinvar']}" target="_blank">Open in ClinVar ↗</a></div>
</div>"""


def build_rows(rows):
    table = []
    for r in rows:
        e, hit = r["event"], r["hit"]
        table.append({
            "Position": e["start"] if e["start"] == e["end"] else f"{e['start']}–{e['end']}",
            "Change": change_text(e),
            "Type": e["type"].capitalize(),
            "HGVS": r["hgvs"],
            "Codon": r["codon"],
            "Effect": r["effect"],
            "Protein": r["protein"],
            "Legacy name": r["legacy"] or "",
            "Known variant": hit["name"] if hit else UNKNOWN,
            "Condition": hit["condition"] if hit else "",
            "ClinVar": r["clinvar"],
        })
    return table


# ------------------------------------------------------------------ interface callbacks
def load_example():
    ref = read_reference()[0]
    st.session_state.ref_text = ref
    st.session_state.test_text = make_examples(ref)[st.session_state.example]


def load_reference_text():
    st.session_state.ref_text = read_reference()[0]


def clear_all():
    st.session_state.ref_text = ""
    st.session_state.test_text = ""


def main():
    st.set_page_config(page_title="VariantLens", page_icon="🧬", layout="wide")
    ref_default, from_file = read_reference()
    st.session_state.setdefault("ref_text", ref_default)
    st.session_state.setdefault("test_text", "")

    with st.sidebar:
        st.header("About")
        st.write("VariantLens compares a test DNA sequence with the **HBB** (beta-globin) "
                 "reference, finds every mutation, names it in HGVS notation, predicts its "
                 "effect on the protein and checks it against known variants.")
        st.markdown("**Reference:** NM_000518.5 CDS, 444 bases")
        st.subheader("Effect types")
        for k, v in EFFECT_HELP.items():
            st.markdown(f"**{k}**: {v}")
        st.warning("Educational tool. Not for clinical or diagnostic use.")

    st.title("🧬 VariantLens")
    st.caption("Compare a reference and a test DNA sequence and label known HBB variants. "
               "Educational use only, not for diagnosis.")

    if not from_file:
        st.info("hbb_cds.txt not found, so a short development reference is used for now.")
    db = get_db()
    if not db:
        st.warning("Known-variant database not found or empty, so every variant will show "
                   "'Not in database'.")

    st.selectbox("Load an example", list(make_examples(ref_default)), key="example",
                 index=None, placeholder="Choose an example test sequence...",
                 on_change=load_example)
    st.text_area("Reference sequence", key="ref_text", height=120,
                 placeholder="Paste reference DNA (A, C, G, T)")
    st.text_area("Test sequence", key="test_text", height=120,
                 placeholder="Paste test DNA (A, C, G, T) or FASTA, or choose an example above")

    c1, c2 = st.columns(2)
    c1.button("Load HBB reference", on_click=load_reference_text, use_container_width=True)
    c2.button("Clear", on_click=clear_all, use_container_width=True)

    if st.button("Analyse", type="primary", use_container_width=True):
        try:
            ref = clean(st.session_state.ref_text)
            test = clean(st.session_state.test_text)
        except (ValueError, TypeError) as err:
            st.error(f"Invalid input: {err}")
            st.stop()

        warn = length_warning(ref, test)
        if warn:
            st.warning(warn)

        events, rows = analyze(ref, test, db)
        if not events:
            st.success("The sequences are identical. No variants found.")
            st.stop()

        st.subheader("Results")
        st.markdown(conclusion_banner_html(rows), unsafe_allow_html=True)

        st.markdown("#### What each variant means")
        for i, r in enumerate(rows, 1):
            st.markdown(variant_card_html(ref, r, i), unsafe_allow_html=True)

        st.markdown("#### Where the variants are")
        st.markdown(highlight_html(ref, events), unsafe_allow_html=True)

        st.markdown("#### Detailed table")
        st.dataframe(
            pd.DataFrame(build_rows(rows)),
            hide_index=True,
            use_container_width=True,
            column_config={"ClinVar": st.column_config.LinkColumn("ClinVar", display_text="Open")},
        )

        with st.expander("Test sequence (cleaned)"):
            st.code(test, language=None)

    st.divider()
    st.caption("Limits: coding sequence only, small local HBB database, "
               "no splice-site, large deletion or structural variants.")


if __name__ == "__main__":
    main()
