# tests/test_app.py - chooses every example in the Streamlit app, clicks Analyse, checks the result.
from pathlib import Path

import pytest
from streamlit.testing.v1 import AppTest

APP = str(Path(__file__).resolve().parent.parent / "app.py")

EXPECT = {
    "Normal (no change)": ("success", None),
    "HbS (sickle cell): c.20A>T": ("table", "c.20A>T"),
    "HbC: c.19G>A": ("table", "c.19G>A"),
    "HbE: c.79G>A": ("table", "c.79G>A"),
    "Codon 39 nonsense: c.118C>T": ("table", "c.118C>T"),
    "Codon 15 nonsense: c.47G>A": ("table", "c.47G>A"),
    "Codon 41/42 deletion: c.126_129del": ("table", "c.126_129del"),
    "Codon 8/9 insertion: c.27dup": ("table", "c.27dup"),
    "Synonymous: c.21G>A": ("table", "c.21G>A"),
    "Start-loss: c.3G>A": ("table", "c.3G>A"),
    "Unknown variant: c.4G>A": ("table", "c.4G>A"),
    "Two mutations: HbS + HbE": ("table", "c.79G>A"),
    "Invalid input (contains X)": ("error", None),
}


def run_example(example):
    at = AppTest.from_file(APP, default_timeout=30).run()
    at.selectbox(key="example").select(example).run()
    next(b for b in at.button if b.label == "Analyse").click().run()
    assert not at.exception
    return at


@pytest.mark.parametrize("example", list(EXPECT))
def test_example(example):
    at = run_example(example)
    kind, variant = EXPECT[example]
    if kind == "success":
        assert at.success
    elif kind == "error":
        assert at.error and "Invalid characters" in at.error[0].value
    else:
        assert variant in list(at.dataframe[0].value["HGVS"])


def test_known_variant_label():
    at = run_example("HbS (sickle cell): c.20A>T")
    df = at.dataframe[0].value
    assert df["Known variant"][0].startswith("HbS")
    assert df["Legacy name"][0].startswith("Glu6Val")


def test_clear_button():
    at = AppTest.from_file(APP, default_timeout=30).run()
    next(b for b in at.button if b.label == "Clear").click().run()
    assert at.text_area(key="ref_text").value == ""


def page_text(at):
    return " ".join(m.value for m in at.markdown)


def test_conclusion_names_known_variant():
    text = page_text(run_example("HbS (sickle cell): c.20A>T"))
    assert "Conclusion" in text and "HbS" in text and "sickle cell disease" in text
    assert "Glutamic acid (Glu) to Valine (Val)" in text


def test_conclusion_flags_uncharacterized():
    text = page_text(run_example("Unknown variant: c.4G>A"))
    assert "uncharacterized" in text and "not that it is harmless" in text


def test_nonsense_card_reports_truncated_length():
    text = page_text(run_example("Codon 39 nonsense: c.118C>T"))
    assert "High predicted impact" in text and "39 of its normal 147" in text


def test_frameshift_card_reports_new_length():
    text = page_text(run_example("Codon 41/42 deletion: c.126_129del"))
    assert "59 amino acids long instead of 147" in text
