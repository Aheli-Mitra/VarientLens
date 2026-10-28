# tests/test_core.py  (run from the project root with: pytest -v)
# Owner: Member A (Data and testing)
import json
import sys
from pathlib import Path

import pytest
from Bio.Seq import Seq

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from core import *  # noqa: E402,F401,F403

REF = clean((ROOT / "data" / "hbb_cds.txt").read_text())
DB = load_db(str(ROOT / "data" / "hbb_variants.json"))

BETA_GLOBIN = ("MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDN"
               "LKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH")


def analyze(test_seq):
    """Run the full pipeline and return one row per variant (what the app shows)."""
    rows = []
    for e in call_variants(REF, test_seq):
        name = to_hgvs(e, REF)
        effect, protein = consequence(REF, e)
        hit = lookup(name, DB)
        rows.append({"event": e, "hgvs": name, "effect": effect, "protein": protein, "hit": hit})
    return rows


# ---------- Part 1: reference sanity (Step 2) ----------

def test_reference_length():
    assert len(REF) == 444 and len(REF) % 3 == 0


def test_reference_start_and_stop():
    assert REF.startswith("ATG")
    assert REF[-3:] in ("TAA", "TAG", "TGA")


def test_reference_translates_to_beta_globin():
    protein = str(Seq(REF).translate())
    assert protein.startswith("MVHLTPEEKS")
    assert protein == BETA_GLOBIN + "*"        # 147 aa, single stop at the end
    assert protein.count("*") == 1


# ---------- Part 2: database integrity (Step 3) ----------

def test_db_entries_have_same_keys():
    keys = set(DB[0])
    assert all(set(e) == keys for e in DB)
    assert {"hgvs", "aliases", "name", "legacy", "protein", "effect", "condition", "source"} <= keys


def test_db_names_are_unique():
    names = [e["hgvs"] for e in DB] + [a for e in DB for a in e["aliases"]]
    assert len(names) == len(set(names))


# How to build each database variant from the reference
BUILDERS = {
    "c.20A>T": lambda: substitute(REF, 20, "A", "T"),
    "c.19G>A": lambda: substitute(REF, 19, "G", "A"),
    "c.79G>A": lambda: substitute(REF, 79, "G", "A"),
    "c.118C>T": lambda: substitute(REF, 118, "C", "T"),
    "c.47G>A": lambda: substitute(REF, 47, "G", "A"),
    "c.126_129del": lambda: delete(REF, 126, 129),
    "c.27dup": lambda: insert_after(REF, 27, REF[26]),
}


def test_every_db_entry_has_a_builder():
    assert {e["hgvs"] for e in DB} == set(BUILDERS)


@pytest.mark.parametrize("entry", DB, ids=[e["hgvs"] for e in DB])
def test_db_entry_round_trip(entry):
    """Build the variant, run the pipeline, and check it finds itself in the database."""
    rows = analyze(BUILDERS[entry["hgvs"]]())
    assert len(rows) == 1
    row = rows[0]
    assert row["hit"] is entry, f"pipeline named it {row['hgvs']}, not found in database"
    assert row["effect"] == entry["effect"]
    if entry["protein"].startswith("p."):
        assert row["protein"] == entry["protein"]


@pytest.mark.parametrize("entry", [e for e in DB if e["aliases"]], ids=lambda e: e["hgvs"])
def test_aliases_are_truly_equivalent(entry):
    """Every alias spelling must give exactly the same mutated sequence as the main name."""
    target = BUILDERS[entry["hgvs"]]()
    for alias in entry["aliases"]:
        body = alias[2:]
        if body.endswith("del"):
            s, e = map(int, body[:-3].split("_"))
            assert delete(REF, s, e) == target, alias
        elif "ins" in body:
            left, bases = body.split("ins")
            assert insert_after(REF, int(left.split("_")[0]), bases) == target, alias


# ---------- Part 3: Step 10 test cases ----------

def test_identical():
    assert call_variants(REF, REF) == []


def test_hbs():
    rows = analyze(substitute(REF, 20, "A", "T"))
    assert len(rows) == 1
    assert rows[0]["event"]["type"] == "substitution"
    assert rows[0]["hgvs"] == "c.20A>T"
    assert (rows[0]["effect"], rows[0]["protein"]) == ("missense", "p.Glu7Val")
    assert rows[0]["hit"]["name"].startswith("HbS")


def test_hbc():
    rows = analyze(substitute(REF, 19, "G", "A"))
    assert len(rows) == 1 and rows[0]["hgvs"] == "c.19G>A"
    assert (rows[0]["effect"], rows[0]["protein"]) == ("missense", "p.Glu7Lys")
    assert rows[0]["hit"]["name"] == "HbC"


def test_hbe():
    rows = analyze(substitute(REF, 79, "G", "A"))
    assert len(rows) == 1 and rows[0]["hgvs"] == "c.79G>A"
    assert (rows[0]["effect"], rows[0]["protein"]) == ("missense", "p.Glu27Lys")
    assert rows[0]["hit"]["name"] == "HbE"


def test_codon39_nonsense():
    rows = analyze(substitute(REF, 118, "C", "T"))
    assert len(rows) == 1
    assert (rows[0]["effect"], rows[0]["protein"]) == ("nonsense", "p.Gln40Ter")


def test_deletion():
    ev = call_variants(REF, delete(REF, 126, 129))
    # Assert type and length, not exact position: the 4 bases sit in a repeat (Step 6).
    assert len(ev) == 1 and ev[0]["type"] == "deletion" and len(ev[0]["ref"]) == 4
    assert consequence(REF, ev[0])[0] == "frameshift"


def test_duplication():
    ev = call_variants(REF, insert_after(REF, 27, REF[26]))
    assert len(ev) == 1 and ev[0]["type"] == "insertion" and len(ev[0]["alt"]) == 1
    assert consequence(REF, ev[0])[0] == "frameshift"


def test_synonymous():
    # c.21G>A: codon 7 GAG > GAA, still Glu
    rows = analyze(substitute(REF, 21, "G", "A"))
    assert len(rows) == 1
    assert (rows[0]["effect"], rows[0]["protein"]) == ("synonymous", "p.Glu7=")


def test_not_in_database():
    # c.4G>A: codon 2 GTG > ATG (Val2Met), not in our table
    rows = analyze(substitute(REF, 4, "G", "A"))
    assert len(rows) == 1 and rows[0]["hgvs"] == "c.4G>A"
    assert rows[0]["effect"] == "missense"
    assert rows[0]["hit"] is None          # app must show "Not in local database (uncharacterized)"


def test_two_substitutions():
    test = substitute(substitute(REF, 20, "A", "T"), 79, "G", "A")
    rows = analyze(test)
    assert [r["hgvs"] for r in rows] == ["c.20A>T", "c.79G>A"]


@pytest.mark.parametrize("bad", ["ACGTX", "ACG T1 23", "ACGTN"])
def test_invalid_characters(bad):
    with pytest.raises(ValueError, match="Invalid characters"):
        clean(bad)


@pytest.mark.parametrize("empty", ["", "   \n  ", ">just a header\n"])
def test_empty_input(empty):
    with pytest.raises(ValueError, match="empty"):
        clean(empty)


def test_fasta_and_lowercase_are_accepted():
    assert clean(">seq1\natg gtg\ncat\n") == "ATGGTGCAT"


# ---------- Part 4: extra edge cases ----------

def test_inframe_deletion():
    ev = call_variants(REF, delete(REF, 100, 102))
    assert len(ev) == 1 and len(ev[0]["ref"]) == 3
    assert consequence(REF, ev[0])[0] == "in-frame deletion"


def test_stop_loss():
    rows = analyze(substitute(REF, 443, "A", "C"))   # TAA > TCA
    assert rows[0]["effect"] == "stop-loss"


def test_clinvar_link_is_search_url():
    link = clinvar_link("c.20A>T")
    assert link.startswith("https://www.ncbi.nlm.nih.gov/clinvar/?term=")
    assert "HBB" in link and " " not in link


def test_substitute_guard_catches_wrong_base():
    with pytest.raises(AssertionError):
        substitute(REF, 20, "G", "T")    # position 20 is A, not G


# ---------- Part 5: merged-code features (3' normalization, Member C's annotation) ----------

def test_deletion_is_3prime_normalized():
    # deleting 124-127, 125-128 or 126-129 gives the same sequence; HGVS names the 3'-most
    for s in (124, 125, 126):
        ev = call_variants(REF, delete(REF, s, s + 3))
        assert to_hgvs(ev[0], REF) == "c.126_129del"


def test_duplication_is_3prime_normalized():
    for p in (26, 27):
        ev = call_variants(REF, insert_after(REF, p, "G"))
        assert to_hgvs(ev[0], REF) == "c.27dup"


def test_normalized_events_describe_the_same_sequence():
    """For every 1-4 base deletion and 1-3 base insertion, the event must rebuild the test sequence."""
    for p in range(1, 441, 7):
        cases = [delete(REF, p, p + k) for k in range(4)]
        cases += [insert_after(REF, p, b) for b in ("G", REF[p:p + 2], REF[p:p + 3], "TA")]
        for test in cases:
            ev = call_variants(REF, test)
            if len(ev) == 1:
                assert apply_event(REF, ev[0]) == test


def test_frameshift_protein_name():
    ev = call_variants(REF, delete(REF, 126, 129))
    assert consequence(REF, ev[0]) == ("frameshift", "p.Phe42Leufs")   # ClinVar: p.Phe42LeufsTer19


def test_start_loss():
    ev = call_variants(REF, substitute(REF, 3, "G", "A"))       # ATG > ATA
    assert consequence(REF, ev[0]) == ("start-loss", "p.Met1?")


def test_legacy_name_uses_hb_numbering():
    ev = call_variants(REF, substitute(REF, 20, "A", "T"))
    assert annotate(REF, ev[0])["legacy_name"] == "Glu6Val"
    assert annotate(REF, ev[0])["codon"] == "GAG > GTG"


def test_lookup_accepts_transcript_prefix():
    assert lookup("NM_000518.5:c.20A>T", DB)["name"].startswith("HbS")


def test_clinvar_direct_link_when_verified():
    entry = lookup("c.20A>T", DB)
    assert clinvar_link("c.20A>T", entry).endswith("/variation/15333/")
