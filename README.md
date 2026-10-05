# Spurious-pockets-dark-proteome

Analysis code for *Most predicted drug pockets in the apicomplexan dark
proteome are artefacts of non-globular structure*.

Every numeric value in the manuscript, figures and table is computed here and
read from `numbers_p2.json`. Nothing is transcribed by hand, so a disagreement
between the text and the data surfaces as a failed assertion rather than as a
silent error.

Data: Zenodo, doi to follow. Companion methods paper:
https://github.com/ktosmanprotein/structural-novelty-bias

## Reproducing

```
python3 add_fold_class.py            # radius of gyration, fold class, lining topology
python3 check_annotation_sources.py  # UniProt against VEuPathDB -> the strict dark set
python3 cross_species_foldclass.py   # the same fold test across 15 genomes
python3 audit_numbers_p2.py          # everything above -> numbers_p2.json
```

Run them in that order; `audit_numbers_p2.py` reads what the first three write.
It also re-derives the figures in the companion methods paper's record and
asserts agreement, so it fails if the two papers drift apart.

## What each script does

| Script | Purpose |
| --- | --- |
| `add_fold_class.py` | Radius of gyration over the confident core (pLDDT >= 70), the ratio to the compact-globule expectation 2.2N^0.38, fold class, and the number of distinct sequence segments lining each pocket |
| `check_annotation_sources.py` | Compares the UniProt name with the VEuPathDB product for all 128 candidates and applies the strict and lenient annotation rules of section 4.3 |
| `cross_species_foldclass.py` | Applies the fold test to the dark proteome of each of the 15 panel genomes and correlates the extended fraction with low-complexity content |
| `audit_numbers_p2.py` | Regenerates every number quoted in the manuscript and writes `numbers_p2.json` |
| `references_p2.py` | The reference list, with the PMID or Crossref record each entry was verified against |
| `run_fpocket_dark128.command` | Runs fpocket over the dark set |
| `research_316_3di_aa.command` | Re-searches the entry set in Foldseek's combined 3Di and amino-acid mode (section 4.6) |

## Not in this repository

The figure-rendering scripts. They restyle values `audit_numbers_p2.py` has
already computed, so they reproduce the pictures rather than the results. The
manuscript builder is also absent: the submitted manuscript was edited by hand
after it was generated, so re-running the builder would not reproduce it.

## Inputs

AlphaFold DB v6 models for *P. falciparum*, the fpocket output, and the
novelty and verification tables of the companion methods paper
(doi:10.5281/zenodo.23045417).
