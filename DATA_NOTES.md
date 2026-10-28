# VariantLens: Data and Testing Notes (Member A)

## 1. Reference: `data/hbb_cds.txt`

- **Source:** HBB transcript **NM_000518.5**, CDS feature (c.1 to c.444). The older NM_000518.4 has the same CDS.
- **Format:** FASTA (one header line, 60 bases per line). `clean()` strips the header automatically.

| Sanity check (Step 2) | Result |
|---|---|
| Length 444, divisible by 3 | Pass |
| Starts with ATG | Pass |
| Ends with stop codon (TAA) | Pass |
| Translates to MVHLTPEEKS... | Pass |
| 147 amino acids, identical to beta-globin protein, single stop at end | Pass |

## 2. Database: `data/hbb_variants.json`

Every entry was checked in three ways on **2026-10-27**:

1. **Reference base matches:** the base at the c. position is the stated reference base.
2. **Protein change matches:** translating the mutated codon gives the stated p. change.
3. **External record exists:** ClinVar, with LOVD / IthaGenes / HbVar as cross-checks.

| HGVS | Known as | Protein | Ref base OK | ClinVar |
|---|---|---|---|---|
| c.20A>T | HbS | p.Glu7Val | A | Pathogenic (VCV 15333) |
| c.19G>A | HbC | p.Glu7Lys | G | Pathogenic (VCV 15126) |
| c.79G>A | HbE | p.Glu27Lys | G | Record found |
| c.118C>T | Codon 39 | p.Gln40Ter | C | Record found (also LOVD) |
| c.47G>A | Codon 15 | p.Trp16Ter | G | Record found |
| c.126_129del | Codon 41/42 (-TTCT) | frameshift (ClinVar: p.Phe42fs) | CTTT | Pathogenic (VCV 15417) |
| c.27dup | Codon 8/9 (+G) | frameshift | G | Record found (also LOVD, HbVar) |

**Additions to the guide's JSON format:**

- **`clinvar_url`** is the direct ClinVar record link where the variation ID was confirmed. It is empty otherwise, and the app falls back to `clinvar_link()` search.
- **Extra aliases for the two indels.** See section 3.

## 3. Key finding: indel aliases

Both indels sit in repeated bases, so several positions give the **same mutated sequence**. I checked this by brute force over every position:

- Deleting 4 bases at **124-127 (TTCT), 125-128 (TCTT) or 126-129 (CTTT)** gives an identical sequence. HGVS uses the 3'-most, `c.126_129del`. LOVD lists all three spellings for the same variant.
- Inserting G after **base 26 or base 27** gives an identical sequence. HGVS calls it `c.27dup`.

**Our aligner picks the 5'-most position:** `c.124_127del` and `c.26_27insG`. The guide's JSON had only `c.27_28insG` as an alias, so without the extra aliases both famous thalassemia variants would show as *"Not in local database"* in the demo.

**Fix applied (two layers):**

1. Member B's `normalize_indel()` now shifts every indel to its 3'-most position, so the app reports the standard names `c.126_129del` and `c.27dup` directly.
2. The aliases `c.124_127del`, `c.125_128del` and `c.26_27insG`, `c.27_28insG` stay in the database as a backup, for example for names typed by hand. The test `test_aliases_are_truly_equivalent` proves every alias gives the same sequence.

**Bug found during merging:** the first version of the normalizer moved multi-base insertions without rotating the inserted bases, which gave a wrong name in 1,551 of 7,902 test insertions. After the fix, all 10,616 tested insertions and deletions are described correctly, checked by `test_normalized_events_describe_the_same_sequence`.

> **For the report and viva:** this is the 3-prime normalization issue from Step 6, shown with a real example, and the reason the stretch goal was implemented.

## 4. Test suite: `tests/test_core.py`

Run from the project root: `pytest -v` (44 core and data tests plus 15 app tests in `tests/test_app.py`, all passing).

- **Part 1, Reference sanity:** length, start/stop, exact protein.
- **Part 2, Database integrity:**
  - Same keys in every entry, and no duplicate names.
  - **Round-trip:** each entry is built from the reference, run through the full pipeline, and must find itself in the database with the right effect and protein.
  - Aliases are truly equivalent.
- **Part 3, Step 10 cases:** all 12 rows of the guide's table.
- **Part 4, Extras:** in-frame deletion, stop-loss, ClinVar link format, and the `substitute()` guard catching a wrong base.

`python tests/make_results_table.py` regenerates `results_table.md` from real runs. Rerun it after any change to `core.py`.

## 5. How the teammates' code was merged

- **B (`core.py`):** validation, comparison and 3' normalization kept as written; insertion rotation fixed.
- **C (`annotation.py`):**
    - The genetic code, `translate()` and the consequence logic moved into `core.py` as `annotate()`.
    - It was adapted to B's event format (`type/start/end/ref/alt` instead of `SNV/DEL/INS` with `pos`).
    - `consequence()` still returns `(effect, protein)` for the tests; `annotate()` also gives the legacy name and codon change.
    - C's idea of stripping `NM_000518.5:` prefixes was added to `lookup()`.
    - C's `substitute/delete/insert_after` had different arguments from the guide's, so the guide's versions were kept (the app and tests use them).
- **D (`app.py`):**
    - Layout, colour highlighting with position numbers, buttons and the results table kept.
    - An example dropdown, Codon / Protein / Legacy name / Condition columns, the length warning and the "uncharacterized" note were added.
    - The fallbacks for an unknown `core.py` API were removed, now that the API is fixed.
    - The page now uses the wide layout so the 11-column table fits.
- **Database:** C's loader expected `hgvs_c` and `clinvar_id` fields, so the merged `lookup()` uses the existing `hgvs`, `aliases` and `clinvar_url` fields instead.
