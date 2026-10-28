# VariantLens: Five-Minute Demo Script

**Before you start:** run `streamlit run app.py`, open the browser tab, and close other windows. Run `pytest` once beforehand to confirm everything passes.

## Minute 1: The problem

> "A single changed letter in DNA can cause a serious disease. In sickle cell disease, one A in the HBB gene becomes a T. Finding such changes by eye in a 444-letter sequence is slow and error-prone. VariantLens finds them automatically, names them in the international standard format, and explains what they do to the protein."

## Minute 2: Sickle cell (HbS)

1. Select **HbS (sickle cell): c.20A>T** and click **Analyse**.
2. Point at the single red base in the first line (position 20).
3. Read the **Conclusion** box aloud, then point out on the variant card:
   - **c.20A>T**: position 20, counted from the A of ATG, changes A to T.
   - **Codon GAG > GTG** and **p.Glu7Val**: glutamic acid becomes valine.
   - **Known as HbS**, with the legacy name Glu6Val.

> "Textbooks say Glu6Val because they skip the first methionine. The modern HGVS standard counts it, so it is Glu7Val. Same mutation, two names."

## Minute 3: Deletion and alignment

1. Select **Codon 41/42 deletion: c.126_129del** and click **Analyse**.
2. Point at the four purple bases.

> "Four bases are missing, so the test sequence is shorter. If we compared letter by letter, every base after the deletion would look different. That would be hundreds of false mutations. So we use global sequence alignment, which lines the sequences up and finds the gap as one event. Four is not a multiple of three, so the reading frame shifts: a frameshift, which causes beta-thalassemia."

**If asked about the position:** "The deleted bases sit in a repeat, so deleting at 124, 125 or 126 gives the same sequence. The aligner picks the leftmost, so we shift it to the 3'-most position, which is the HGVS rule. That gives c.126_129del, the name ClinVar uses. The protein name p.Phe42Leufs means phenylalanine 42 becomes leucine and the frame shifts."

## Minute 4: Honest reporting

1. Select **Unknown variant: c.4G>A** and click **Analyse**.

> "This mutation is not in our database, so the app says uncharacterized. It never says benign, because absence from a small table is not evidence of safety. The ClinVar link lets you check the full database."

2. Optionally select **Start-loss: c.3G>A**: the ATG start codon is destroyed, so the protein may not be made at all.
3. Optionally select **Invalid input (contains X)** to show the clear error message.

## Minute 5: Limitations and future work

> "VariantLens is educational, not diagnostic. It only covers the coding sequence, so splice mutations common in India, like IVS-I-5, are not detected, and its database has seven variants. Future work includes splice-site detection with the genomic sequence and live ClinVar lookup. All 63 automated tests pass."

## Backup plan

If the app fails to start, show the screenshots in `screenshots/` and run `pytest -v` in the terminal.
