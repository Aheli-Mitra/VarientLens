| # | Test case | Expected | Actual output | Result |
|---|---|---|---|---|
| 1 | Reference vs itself | No variants | No variants | PASS |
| 2 | HbS c.20A>T | 1 substitution, missense, p.Glu7Val, HbS | c.20A>T, missense, p.Glu7Val, HbS (sickle cell variant) | PASS |
| 3 | HbC c.19G>A | 1 substitution, missense, p.Glu7Lys, HbC | c.19G>A, missense, p.Glu7Lys, HbC | PASS |
| 4 | HbE c.79G>A | 1 substitution, missense, p.Glu27Lys, HbE | c.79G>A, missense, p.Glu27Lys, HbE | PASS |
| 5 | Codon 39 c.118C>T | nonsense, p.Gln40Ter | c.118C>T, nonsense, p.Gln40Ter, Codon 39 (C>T) | PASS |
| 6 | Codon 15 c.47G>A | nonsense, p.Trp16Ter | c.47G>A, nonsense, p.Trp16Ter, Codon 15 (G>A) | PASS |
| 7 | Delete bases 126-129 | c.126_129del, frameshift, Codon 41/42 | c.126_129del, frameshift, p.Phe42Leufs, Codon 41/42 (-TTCT) | PASS |
| 8 | Duplicate base 27 | c.27dup, frameshift, Codon 8/9 | c.27dup, frameshift, p.Ser10Valfs, Codon 8/9 (+G) | PASS |
| 9 | Synonymous c.21G>A | synonymous, p.Glu7= | c.21G>A, synonymous, p.Glu7=, Not in local database | PASS |
| 10 | Not in table c.4G>A | missense, Not in local database | c.4G>A, missense, p.Val2Met, Not in local database | PASS |
| 11 | Two substitutions (c.20 + c.79) | 2 separate events | c.20A>T, missense, p.Glu7Val, HbS (sickle cell variant); c.79G>A, missense, p.Glu27Lys, HbE | PASS |
| 12 | Start codon c.3G>A | start-loss, p.Met1? | c.3G>A, start-loss, p.Met1?, Not in local database | PASS |
| 13 | In-frame deletion 100-102 | in-frame deletion | c.102_104del, in-frame deletion, p.Val35del, Not in local database | PASS |
| 14 | Invalid characters 'ACGTX' | Clear error, no crash | Error: Invalid characters found: ['X']. Only A, C, G and T are allowed. | PASS |
| 15 | Empty input | Clear error, no crash | Error: Sequence is empty. | PASS |
