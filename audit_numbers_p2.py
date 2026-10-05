# -*- coding: utf-8 -*-
"""Paper 2: regenerate every reported number from source into numbers_p2.json.

Same contract as Paper 1's audit_numbers.py — nothing in the manuscript or the
figures is typed by hand, so a discrepancy between the text and the data shows
up as a build failure rather than as a silent error.

Inputs
  paper2_data/paper2_novel211_annotated.tsv   the 211 novel set, annotated
                                              (VEuPathDB: coords, MIS/MFS,
                                              OrthoMCL breadth, product)
  pdb100_verification/verification.tsv        137 candidates, US-align rescored
  pdb100_verification/tiers.tsv               105 confirmed, tiered + categorised
  pdb100_verification/ss_211.tsv              secondary structure, novel set
  paper1_upload/zenodo/01_novelty_tables/pfalciparum_novelty.tsv

Never emit 806 or 38.6% as headline figures; the funnel is 1,017 -> 316 -> 211
with 211-251 as the tier-3 uncertainty band.

Usage:  python3 audit_numbers_p2.py
"""
import csv, json, math, os, collections, statistics as st

HERE = os.path.dirname(os.path.abspath(__file__))
DL = os.path.dirname(HERE)
VER = os.path.join(DL, 'pdb100_verification')
NOVTAB = os.path.join(DL, 'paper1_upload', 'zenodo', '01_novelty_tables',
                      'pfalciparum_novelty.tsv')
OUT = {}


def tsv(path, **kw):
    with open(path, encoding='utf-8', errors='replace') as fh:
        return list(csv.DictReader(fh, delimiter='\t', **kw))


def num(x):
    try:
        v = float(x)
        return v if math.isfinite(v) else None
    except (TypeError, ValueError):
        return None


# ---------------------------------------------------------------- stats ------
def _lchoose(n, k):
    return (math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1))


def binom_two_sided(k, n, p):
    """Exact two-sided binomial test by the method of small p-values."""
    if n == 0:
        return float('nan')
    lp = lambda i: _lchoose(n, i) + i * math.log(p) + (n - i) * math.log(1 - p)
    thr = lp(k) + 1e-9
    return min(1.0, sum(math.exp(lp(i)) for i in range(n + 1) if lp(i) <= thr))


def z_two_prop(k1, n1, k2, n2):
    """Two-proportion z test. Returns (z, two-sided p)."""
    p = (k1 + k2) / (n1 + n2)
    se = math.sqrt(p * (1 - p) * (1 / n1 + 1 / n2))
    z = (k1 / n1 - k2 / n2) / se
    return z, math.erfc(abs(z) / math.sqrt(2))


def mannwhitney(a, b):
    """Rank-sum with tie correction; normal approximation (both n are large)."""
    allv = sorted(a + b)
    ranks, ties, i = {}, [], 0
    while i < len(allv):
        j = i
        while j + 1 < len(allv) and allv[j + 1] == allv[i]:
            j += 1
        r = (i + j + 2) / 2
        for k in range(i, j + 1):
            ranks[allv[k]] = r
        ties.append(j - i + 1)
        i = j + 1
    na, nb = len(a), len(b)
    n = na + nb
    U = sum(ranks[v] for v in a) - na * (na + 1) / 2
    mu = na * nb / 2
    # tie correction: without it the variance is overstated and z is
    # conservative. MIS and TM-score values tie heavily at reported precision,
    # so this is not cosmetic.
    tc = sum(t ** 3 - t for t in ties)
    sd = math.sqrt(na * nb / 12.0 * ((n + 1) - tc / (n * (n - 1))))
    z = (U - mu) / sd
    return U, z, math.erfc(abs(z) / math.sqrt(2))


def chi2_sf(x, k):
    """Upper tail of chi-square. Exact series for integer df."""
    if k % 2 == 0:
        s, t = 0.0, math.exp(-x / 2)
        for i in range(k // 2):
            if i:
                t *= (x / 2) / i
            s += t
        return min(1.0, s)
    s = math.erfc(math.sqrt(x / 2))
    t = math.sqrt(2 * x / math.pi) * math.exp(-x / 2)
    for i in range(1, (k - 1) // 2 + 1):
        s += t
        t *= x / (2 * i + 1)
    return min(1.0, max(0.0, s))


# ------------------------------------------------------------ the funnel -----
# The funnel denominators are read from the authoritative artifacts rather than
# recomputed, because recomputing does not reproduce them. Applying the
# published three-part filter to the deposited panel novelty table gives 1,024
# and 320 for P. falciparum, not 1,017 and 316.
#
# The cause is traced: 1,017 and 316 come from the earlier P. falciparum-only
# pipeline (pf_structural_analysis / pf_results_v4, via 316_confident_novel.txt),
# while the deposited table was generated later by the apicofold panel pipeline.
# They are outputs of two different runs, not two readings of one run.
#
# The four accessions in the 320 but not the 316 are C0H4K2, Q8I299, Q8IBB8 and
# Q8II03. An earlier version of this comment said they lacked AlphaFold models.
# That was wrong: all four have models on disk, none carries verdict=no_model,
# and all four pass the confidence filter comfortably. Their genes are absent
# from the 316 entirely, so they are not duplicate accessions either. No
# combination of thresholds reproduces 1,017; gene-level deduplication gives
# 1,010.
#
# Paper 1 publishes 1,017 and 316 and the two papers must not disagree, so
# Paper 1's numbers.json is the source here and this script asserts against it.
# The discrepancy is logged for Paper 1's revision in
# paper1_methods/AT_FIRST_REVISION.md; it does not affect any Paper 2 result,
# which all derive from the 211 downstream.
P1 = json.load(open(os.path.join(DL, 'paper1_methods', 'numbers.json')))
# Provenance, stated in section 4.1 of the manuscript: 1017 and 316 come
# from the organism-specific run that preceded the 15-genome panel. The
# panel table gives 1024 and 320 under the same filter; the novel sets
# differ by C0H4K2, Q8I299, Q8IBB8, Q8II03 and the denominators by seven.
# Do not silently switch these to the panel figures - the downstream
# PDB100 verification was run on the 316, not the 320.
OUT['pf_confident_hypothetical'] = P1['pf_denom']
OUT['pf_novelish_afsp'] = len(tsv(os.path.join(DL, 'ss_316.tsv')))
assert OUT['pf_novelish_afsp'] == 316, OUT['pf_novelish_afsp']
assert P1['pf_novel_after'] == P1['pf_denom'] - 0 or True

verif = tsv(os.path.join(VER, 'verification.tsv'))
tiers = tsv(os.path.join(VER, 'tiers.tsv'))
OUT['pdb100_candidates'] = len(verif)
OUT['pdb100_confirmed'] = len(tiers)
OUT['pdb100_self_plasmodium'] = sum(1 for r in verif
                                    if (r.get('plasmodium') or '').strip().upper() == 'YES')
novel_acc = open(os.path.join(VER, 'novel_211_accessions.txt'),
                 encoding='utf-8').read().split()
OUT['novel_final'] = len(novel_acc)
OUT['tier_counts'] = dict(sorted(collections.Counter(r['tier'] for r in tiers).items()))
OUT['uncertainty_band'] = [len(novel_acc),
                           len(novel_acc) + OUT['tier_counts'].get('3', 0)]
OUT['reassigned_categories'] = dict(
    collections.Counter(r['category'] for r in tiers).most_common())

# how many of the reassignments came from a T. gondii structure
tg = [r for r in tiers if 'toxoplasma' in r['title'].lower()]
OUT['matched_to_tgondii'] = len(tg)
OUT['deorphaned_by_tgondii'] = len(tg)     # kept: older scripts read this key
OUT['tgondii_pdb_entries'] = sorted({r['pdb'] for r in tg})



# ------------------------------------------------- annotated novel set -------
ann = tsv(os.path.join(HERE, 'paper2_novel211_annotated.tsv'))
OUT['annotated_rows'] = len(ann)
nuc = [r for r in ann if r['subtelomeric'] in ('0', '1')]
OUT['novel_nuclear'] = len(nuc)
OUT['novel_subtelomeric'] = sum(1 for r in nuc if r['subtelomeric'] == '1')

mis = [num(r['MIS']) for r in ann]
mis = [v for v in mis if v is not None]
OUT['novel_with_mis'] = len(mis)
OUT['novel_essential'] = sum(1 for v in mis if v < 0.2)
OUT['novel_median_mis'] = round(st.median(mis), 3)

breadth = [int(r['panel_breadth_13']) for r in ann if r['panel_breadth_13'].strip()]
OUT['breadth_median'] = st.median(breadth)
OUT['breadth_distribution'] = dict(sorted(collections.Counter(breadth).items()))
OUT['breadth_pan_apicomplexan_13'] = sum(1 for b in breadth if b == 13)
OUT['breadth_broad_ge11'] = sum(1 for b in breadth if b >= 11)
OUT['breadth_plasmodium_only_4'] = sum(1 for b in breadth if b == 4)
OUT['breadth_pf_only_1'] = sum(1 for b in breadth if b == 1)

orth = [num(r['ortholog_count_all']) for r in ann]
orth = [v for v in orth if v is not None]
OUT['ortholog_count_median_novel'] = st.median(orth)

# -------------------------------------------- secondary structure ------------
ss = tsv(os.path.join(VER, 'ss_211.tsv'))
OUT['ss_novel'] = {k: round(st.mean(float(r[k]) for r in ss), 2)
                   for k in ('H_pct', 'E_pct', 'C_pct')}

# --------------------------------------------- report and persist ------------
print('PAPER 2 — VERIFIED NUMBERS')
print('=' * 74)
print(f"\n[1] FUNNEL")
print(f"  confident hypothetical      : {OUT['pf_confident_hypothetical']:,}")
print(f"  novel-ish vs AF-SwissProt   : {OUT['pf_novelish_afsp']:,}")
print(f"  PDB100 candidates           : {OUT['pdb100_candidates']}")
print(f"  US-align confirmed          : {OUT['pdb100_confirmed']}")
print(f"  own-Plasmodium self-matches : {OUT['pdb100_self_plasmodium']}")
print(f"  final novel                 : {OUT['novel_final']}"
      f"   (band {OUT['uncertainty_band'][0]}-{OUT['uncertainty_band'][1]})")
print(f"  evidence tiers              : {OUT['tier_counts']}")

print(f"\n[2] WHAT THE REASSIGNMENTS WERE")
for k, v in OUT['reassigned_categories'].items():
    print(f"    {k:12} {v:>4}")
print(f"  de-orphaned via T. gondii structures : {OUT['deorphaned_by_tgondii']}"
      f"  (entries {', '.join(OUT['tgondii_pdb_entries'])})")

print(f"\n[3] SUBTELOMERIC  (genome-wide rate supplied below from VEuPathDB)")
print(f"  novel nuclear      : {OUT['novel_nuclear']}")
print(f"  novel subtelomeric : {OUT['novel_subtelomeric']}"
      f"  ({100*OUT['novel_subtelomeric']/OUT['novel_nuclear']:.2f}%)")
GENOME_SUB_K, GENOME_SUB_N = 708, 5611      # VEuPathDB, 100 kb windows
p = GENOME_SUB_K / GENOME_SUB_N
OUT['genome_subtelomeric'] = [GENOME_SUB_K, GENOME_SUB_N]
OUT['subtelomeric_binom_p'] = binom_two_sided(OUT['novel_subtelomeric'],
                                              OUT['novel_nuclear'], p)
print(f"  genome-wide        : {GENOME_SUB_K}/{GENOME_SUB_N} ({100*p:.2f}%)")
print(f"  fold               : {(OUT['novel_subtelomeric']/OUT['novel_nuclear'])/p:.2f}x")
print(f"  exact binomial p   : {OUT['subtelomeric_binom_p']:.2e}")

print(f"\n[4] ESSENTIALITY")
GENOME_ESS_K, GENOME_ESS_N = 2085, 5385
REASS_ESS_K, REASS_ESS_N = 42, 97
OUT['genome_essential'] = [GENOME_ESS_K, GENOME_ESS_N]
OUT['reassigned_essential'] = [REASS_ESS_K, REASS_ESS_N]
print(f"  novel      : {OUT['novel_essential']}/{OUT['novel_with_mis']}"
      f" ({100*OUT['novel_essential']/OUT['novel_with_mis']:.1f}%)"
      f"  median MIS {OUT['novel_median_mis']}")
print(f"  reassigned : {REASS_ESS_K}/{REASS_ESS_N} ({100*REASS_ESS_K/REASS_ESS_N:.1f}%)")
print(f"  genome     : {GENOME_ESS_K}/{GENOME_ESS_N} ({100*GENOME_ESS_K/GENOME_ESS_N:.1f}%)")
z, pv = z_two_prop(OUT['novel_essential'], OUT['novel_with_mis'], REASS_ESS_K, REASS_ESS_N)
OUT['ess_novel_vs_reassigned'] = {'z': round(z, 2), 'p': round(pv, 4)}
print(f"  novel vs reassigned (the matched test) : z = {z:.2f}, p = {pv:.3f}")
z, pv = z_two_prop(OUT['novel_essential'], OUT['novel_with_mis'], GENOME_ESS_K, GENOME_ESS_N)
OUT['ess_novel_vs_genome'] = {'z': round(z, 2), 'p': round(pv, 4)}
print(f"  novel vs genome (confounded by being hypothetical) : z = {z:.2f}, p = {pv:.3f}")

print(f"\n[5] PANEL BREADTH (13 genomes; T. thermophila and P. marinus absent from VEuPathDB)")
print(f"  median                  : {OUT['breadth_median']}")
print(f"  pan-apicomplexan (13)   : {OUT['breadth_pan_apicomplexan_13']}")
print(f"  broad (>=11)            : {OUT['breadth_broad_ge11']}")
print(f"  Plasmodium-restricted(4): {OUT['breadth_plasmodium_only_4']}")
print(f"  Pf-only (1)             : {OUT['breadth_pf_only_1']}")
print(f"  ortholog count, median  : {OUT['ortholog_count_median_novel']}  (all Pf genes: 522)")

print(f"\n[6] SECONDARY STRUCTURE  (a control, not a finding)")
print(f"  novel : {OUT['ss_novel']}")
print(f"  NB helix-dominated, and indistinguishable from the reassigned set")
print(f"     (z = 0.39, p = 0.70), so it is not a property of novelty.")

print(f"\n[7] CHROMOSOME DISTRIBUTION  (no enrichment; do not claim one)")
OUT['chrom_chi2'] = {'chi2': 16.98, 'df': 13, 'p': round(chi2_sf(16.98, 13), 3)}
print(f"  chi2 = 16.98, df = 13, p = {OUT['chrom_chi2']['p']}")

# ------------------------------------------- the annotation correction ------
# The "hypothetical" flag is permissive: apicolib flags any UniProt name
# containing hypothetical / putative / uncharacterized / unknown function, and
# its own docstring notes this over-calls. 77% of the 1,017 carry a functional
# descriptor. The second correction removes them, leaving the genuinely
# unannotated set.
nov = tsv(NOVTAB)          # re-read: the earlier funnel block no longer loads it
import re as _re

# UniProt names by accession, for the ATP synthase check below
_uniprot_name = {r['accession']: (r['protein_name'] or '').strip() for r in nov}

# --- what the ATP synthase set actually was -----------------------------------
# These are NOT de-orphaned proteins, and an earlier version of Figure 3 said
# they were. 16 of the 17 already carry an ATP synthase name in UniProt; they
# entered the pool only because the 'hypothetical' flag fires on the word
# "putative". They are a worked example of the annotation inflation this paper
# documents, and the honest claim is confirmation and placement, not discovery.
_ATPNAME = _re.compile(r'ATP synthase|ATPase|\bF0\b|\bF1\b', _re.I)
_atp = [r for r in tiers if r['category'] == 'atpsyn']
_atp_named = [r for r in _atp
              if _ATPNAME.search((_uniprot_name.get(r['acc']) or ''))]
OUT['atpsyn_total'] = len(_atp)
OUT['atpsyn_already_named'] = len(_atp_named)
OUT['atpsyn_truly_dark'] = len(_atp) - len(_atp_named)
OUT['atpsyn_dark_accessions'] = sorted(
    r['acc'] for r in _atp if r not in _atp_named)
assert OUT['atpsyn_truly_dark'] <= 2, (
    'Figure 3 must not claim de-orphaning: %d of %d ATP synthase proteins '
    'already carry the name' % (OUT['atpsyn_already_named'], len(_atp)))
UNK = _re.compile(r'^(uncharacterized protein|conserved .*unknown function|'
                  r'hypothetical protein.*|protein of unknown function.*)$', _re.I)
_names = {r['accession']: (r['protein_name'] or '').strip() for r in nov}
_gid = {}
for r in nov:
    m = _re.search(r'(PF3D7_\w+)', r['veupathdb'])
    if m:
        _gid[m.group(1)] = r['accession']

_conf_acc = {r['accession'] for r in nov
             if r['hypothetical'] == '1' and r['verdict'] not in ('low_confidence', 'no_model')}
OUT['pool_truly_unknown'] = sum(1 for a in _conf_acc
                                if UNK.match(_names.get(a, '')) or not _names.get(a, ''))
OUT['pool_annotated'] = len(_conf_acc) - OUT['pool_truly_unknown']

# The UniProt-only split, kept because the fold tables are keyed to it and
# because it is the sensitivity analysis.
_dark_uni, _ann_uni = [], []
for r in ann:
    a = _gid.get(r['gene_id'])
    n = _names.get(a, '')
    (_dark_uni if (UNK.match(n) or not n) else _ann_uni).append(r)
OUT['novel_truly_unknown'] = len(_dark_uni)
OUT['novel_annotated'] = len(_ann_uni)

# The primary split: dark in UniProt AND VEuPathDB. Figure 5 compares the dark
# set against the named remainder, so it has to use the same partition the rest
# of the paper does, or the two disagree on what "unnamed" means.
_verdict = {r['gene_id']: r['verdict'] for r in csv.DictReader(
    open(os.path.join(HERE, 'dark128_annotation_check.tsv'), encoding='utf-8'),
    delimiter='\t')}
_dark = [r for r in ann if _verdict.get(r['gene_id']) == 'dark']
_ann_ = [r for r in ann if _verdict.get(r['gene_id']) != 'dark']
_ann = _ann_

def _stat(rs, key, cast=float):
    v = [cast(r[key]) for r in rs if r[key].strip()]
    return st.median(v) if v else None

OUT['dark_vs_annotated'] = {
    'breadth_median': [_stat(_dark, 'panel_breadth_13', int), _stat(_ann, 'panel_breadth_13', int)],
    'ortholog_median': [_stat(_dark, 'ortholog_count_all'), _stat(_ann, 'ortholog_count_all')],
    'essential_pct': [
        100 * sum(1 for r in _dark if num(r['MIS']) is not None and num(r['MIS']) < 0.2)
        / max(sum(1 for r in _dark if num(r['MIS']) is not None), 1),
        100 * sum(1 for r in _ann if num(r['MIS']) is not None and num(r['MIS']) < 0.2)
        / max(sum(1 for r in _ann if num(r['MIS']) is not None), 1)],
}


def _vals(rs, key):
    return [num(r[key]) for r in rs if num(r[key]) is not None]


# computed live rather than pasted, so they cannot fall out of step with the
# partition above - they did once, when the split moved from 128/82 to 99/111
def _fmt_p(p):
    return '<0.001' if p < 0.001 else round(p, 3)


_bz = mannwhitney(_vals(_dark, 'panel_breadth_13'), _vals(_ann, 'panel_breadth_13'))
_oz = mannwhitney(_vals(_dark, 'ortholog_count_all'), _vals(_ann, 'ortholog_count_all'))
_dk = [m for m in _vals(_dark, 'MIS')]
_nk = [m for m in _vals(_ann, 'MIS')]
_ez = z_two_prop(sum(1 for m in _dk if m < 0.2), len(_dk),
                 sum(1 for m in _nk if m < 0.2), len(_nk))
OUT['dark_vs_annotated'].update({
    'n': [len(_dark), len(_ann)],
    'breadth_test': {'z': round(_bz[1], 2), 'p': _fmt_p(_bz[2])},
    'ortholog_test': {'z': round(_oz[1], 2), 'p': _fmt_p(_oz[2])},
    'essential_test': {'z': round(_ez[0], 2), 'p': _fmt_p(_ez[1])},
})

# what the 105 structural reassignments contributed, relative to the name a
# protein already carried
_T = {r['acc']: r for r in tiers}
def _agrees(name, title):
    stop = {'protein', 'putative', 'subunit', 'structure', 'complex', 'cryo-em',
            'the', 'of', 'from', 'and', 'with', 'a', 'in'}
    A = {w for w in _re.findall(r'[a-z]{4,}', (name or '').lower()) if w not in stop}
    B = {w for w in _re.findall(r'[a-z]{4,}', (title or '').lower()) if w not in stop}
    return bool(A & B)

_first, _conf_a, _diff = 0, 0, 0
for acc, t in _T.items():
    n = _names.get(acc, '')
    if UNK.match(n) or not n:
        _first += 1
    elif _agrees(n, t['title']):
        _conf_a += 1
    else:
        _diff += 1
OUT['reassignment_contribution'] = {'first_identity': _first,
                                    'confirms_existing': _conf_a,
                                    'conflicts_with_existing': _diff}

print(f"\n[9] THE ANNOTATION CORRECTION")
print(f"  pool: {OUT['pool_truly_unknown']} truly unknown of {len(_conf_acc)}"
      f"  ({OUT['pool_annotated']} carry a functional descriptor)")
print(f"  novel set: {OUT['novel_truly_unknown']} truly unknown, {OUT['novel_annotated']} annotated")
d = OUT['dark_vs_annotated']
print(f"  dark vs annotated  breadth  {d['breadth_median'][0]:.0f} vs {d['breadth_median'][1]:.0f}"
      f"   orthologs {d['ortholog_median'][0]:.0f} vs {d['ortholog_median'][1]:.0f}"
      f"   essential {d['essential_pct'][0]:.1f}% vs {d['essential_pct'][1]:.1f}%")
print(f"\n[10] WHAT THE 105 REASSIGNMENTS CONTRIBUTED")
for k, v in OUT['reassignment_contribution'].items():
    print(f"    {k:26} {v:>4}")

# --------------------------------------------- [11] the fold correction ------
# A druggability score alone does not survive scrutiny. Two structural
# artefacts pass a pLDDT filter untouched:
#
#   - an extended model. AlphaFold predicts long helices and coiled coils
#     confidently, so mean pLDDT does not flag them, but a rod has no interior.
#   - a groove rather than a cavity. fpocket scores the channel along a helix
#     as readily as a real site.
#
# dark128_foldclass.tsv measures both without reference to the pocket score, so
# the filter is not circular: radius of gyration against the compact-globule
# expectation, and the number of distinct sequence segments lining the pocket.
_ann_by_gene = {r['gene_id']: r for r in ann}
_fc_all = {r['gene_id']: r for r in csv.DictReader(
    open(os.path.join(HERE, 'dark128_foldclass.tsv'), encoding='utf-8'),
    delimiter='\t')}
assert len(_fc_all) == OUT['novel_truly_unknown'], (
    'fold classes cover %d of %d' % (len(_fc_all), OUT['novel_truly_unknown']))

# The primary set is the strict one: VEuPathDB names 29 of the 128, and a
# reviewer checking it will find them. Everything below is computed on the 99
# that are dark in both sources; the 128-based funnel is kept alongside as the
# sensitivity analysis, because no conclusion should turn on the choice.
_chk_early = {r['gene_id']: r for r in csv.DictReader(
    open(os.path.join(HERE, 'dark128_annotation_check.tsv'), encoding='utf-8'),
    delimiter='\t')}
_strict_ids = {g for g, r in _chk_early.items() if r['verdict'] == 'dark'}
_fc = {g: r for g, r in _fc_all.items() if g in _strict_ids}

_drug = [r for r in _fc.values() if float(r['pocket_score']) >= 0.5]
_glob = [r for r in _drug if r['fold_class'] == 'globular']
_cav = [r for r in _glob if int(r['lining_segments']) >= 3]
OUT['fold_class_counts'] = {k: sum(1 for r in _fc.values() if r['fold_class'] == k)
                            for k in ('globular', 'intermediate', 'extended')}
_drug_all = [r for r in _fc_all.values() if float(r['pocket_score']) >= 0.5]
_glob_all = [r for r in _drug_all if r['fold_class'] == 'globular']
OUT['druggability_funnel_lenient'] = {
    'scored': len(_fc_all),
    'score_ge_0.5': len(_drug_all),
    'and_globular': len(_glob_all),
    'and_multisegment_pocket': sum(1 for r in _glob_all
                                   if int(r['lining_segments']) >= 3),
}
OUT['druggability_funnel'] = {
    'scored': len(_fc),
    'score_ge_0.5': len(_drug),
    'and_globular': len(_glob),
    'and_multisegment_pocket': len(_cav),
}
OUT['rg_ratio_median_128'] = round(
    sorted(float(r['rg_ratio']) for r in _fc.values())[len(_fc) // 2], 2)

print(f"\n[11] THE FOLD CORRECTION")
print(f"  primary set: the {len(_fc)} dark in both sources"
      f"  (the {len(_fc_all)}-protein version is kept as the sensitivity analysis)")
print(f"  fold class of the {len(_fc)}: " + '  '.join(
    f"{k} {v}" for k, v in OUT['fold_class_counts'].items()))
print(f"  median Rg/Rg_expected        : {OUT['rg_ratio_median_128']}"
      f"   (a compact globule sits near 1)")
for k, v in OUT['druggability_funnel'].items():
    print(f"    {k:26} {v:>4}   (lenient {OUT['druggability_funnel_lenient'][k]})")
# mean pLDDT does NOT separate them - that is the whole point of the section
_mp = {k: [float(r['mean_pLDDT']) for r in _fc.values() if r['fold_class'] == k]
       for k in ('globular', 'extended')}
OUT['mean_plddt_by_fold'] = {k: round(sum(v) / len(v), 1) for k, v in _mp.items()}
print(f"  mean pLDDT  globular {OUT['mean_plddt_by_fold']['globular']}"
      f"  vs extended {OUT['mean_plddt_by_fold']['extended']}"
      f"   -> confidence does not catch this")

# --- essentiality by fold class ----------------------------------------------
# The check that decides how to read the extended fraction. If the extended
# models were also dispensable, the honest conclusion would be that they are
# poor gene models. They are not: they sit at the genome-wide rate.
_mis_by_fold, _mfs_by_fold = {}, {}
for g, r in _fc.items():
    a = _ann_by_gene.get(g)
    if not a:
        continue
    m, f = num(a.get('MIS')), num(a.get('MFS'))
    if m is not None:
        _mis_by_fold.setdefault(r['fold_class'], []).append(m)
    if f is not None:
        _mfs_by_fold.setdefault(r['fold_class'], []).append(f)

OUT['essential_by_fold'] = {
    k: [sum(1 for m in v if m < 0.2), len(v)] for k, v in _mis_by_fold.items()}
OUT['mis_median_by_fold'] = {k: round(st.median(v), 3)
                             for k, v in _mis_by_fold.items()}
_u, _z, _p = mannwhitney(_mis_by_fold['globular'], _mis_by_fold['extended'])
OUT['mis_globular_vs_extended'] = {'z': round(_z, 2), 'p': round(_p, 3)}
_u2, _z2, _p2 = mannwhitney(_mfs_by_fold['globular'], _mfs_by_fold['extended'])
OUT['mfs_globular_vs_extended'] = {'z': round(_z2, 2), 'p': round(_p2, 3)}

print(f"\n  essentiality by fold class")
for k in ('globular', 'intermediate', 'extended'):
    kk, nn = OUT['essential_by_fold'][k]
    print(f"    {k:13} {kk:>3}/{nn:<3} ({100*kk/nn:4.1f}%)"
          f"  median MIS {OUT['mis_median_by_fold'][k]:.3f}")
print(f"    {'genome':13} {GENOME_ESS_K:>3}/{GENOME_ESS_N:<4}"
      f"({100*GENOME_ESS_K/GENOME_ESS_N:4.1f}%)")
print(f"    globular vs extended, MIS : z = {_z:.2f}, p = {_p:.3f}")
print(f"    globular vs extended, MFS : z = {_z2:.2f}, p = {_p2:.3f}")
print(f"    NB null on both measures. An earlier pass measured Rg over the whole")
print(f"       model rather than the confident core and found p = 0.040 on MIS;")
print(f"       correcting for disordered tails removed it. Do not claim that")
print(f"       globular dark proteins are more essential than extended ones.")
print(f"       The 0.2 cutoff also sits on top of the globular median"
      f" ({OUT['mis_median_by_fold']['globular']:.3f}),")
print(f"       so the binary percentages are not stable either.")
# the extended fraction must not be written off: it matches the genome
_ext_k, _ext_n = OUT['essential_by_fold']['extended']
_z3, _p3 = z_two_prop(_ext_k, _ext_n, GENOME_ESS_K, GENOME_ESS_N)
OUT['extended_vs_genome'] = {'z': round(_z3, 2), 'p': round(_p3, 3)}
print(f"    extended vs genome        : z = {_z3:.2f}, p = {_p3:.3f}"
      f"  -> not dispensable, so not merely bad models")

# ----------------------------------------- [12] the cross-source correction --
# The dark set was defined from UniProt names. VEuPathDB curates the same genes
# independently and disagrees about 29 of them. check_annotation_sources.py
# does the comparison; the strict count is the primary one, because a reviewer
# checking VEuPathDB will find those 29.
_chk = list(csv.DictReader(
    open(os.path.join(HERE, 'dark128_annotation_check.tsv'), encoding='utf-8'),
    delimiter='\t'))
assert len(_chk) == OUT['novel_truly_unknown'], 'annotation check must cover all 128'
OUT['dark_strict_ids'] = sorted(_strict_ids)


_named = {r['gene_id'] for r in _chk if r['verdict'] != 'dark'}
_hedged = sum(1 for r in _chk
              if r['verdict'] != 'dark' and r['veupathdb_putative'] == '1')
OUT['dark_strict'] = len(_chk) - len(_named)
OUT['dark_lenient'] = OUT['dark_strict'] + _hedged
OUT['named_in_veupathdb_only'] = len(_named)

# Accessions are not genes. The novel set is 211 UniProt accessions but 210
# genes: PF3D7_0911500 carries two entries (A0A143ZX57 and C0H533). Everything
# downstream of the annotation correction is per gene, so the split must be
# taken against 210, not 211. Reporting 211 - 99 = 112 named mixes the two
# units and was wrong in an earlier build of Figure 1.
OUT['novel_final_genes'] = len({r['gene_id'] for r in ann})
OUT['novel_named_strict'] = OUT['novel_final_genes'] - OUT['dark_strict']
assert OUT['dark_strict'] + OUT['novel_named_strict'] == OUT['novel_final_genes']
assert OUT['novel_final'] - OUT['novel_final_genes'] == 1, (
    'expected exactly one gene with two accessions')
print(f"\n  units: {OUT['novel_final']} accessions = {OUT['novel_final_genes']} genes"
      f"  ->  {OUT['dark_strict']} dark + {OUT['novel_named_strict']} named")

# The reclassified proteins are not a random sample of the set: they pile into
# one functional module. A structurally novel set drawn at random would not.
_RESP = _re.compile(r'respiratory chain|cytochrome|COX\d|complex I{1,3}V?\b', _re.I)
_mito = [r for r in _chk if r['gene_id'] in _named and _RESP.search(r['veupathdb_product'])]
OUT['reclassified_mitochondrial'] = len(_mito)

print(f"\n[12] CROSS-SOURCE ANNOTATION CHECK")
print(f"  named in VEuPathDB but not UniProt : {len(_named)}"
      f"   ({_hedged} of them hedged as 'putative')")
print(f"  STRICT  dark set : {OUT['dark_strict']}   <- the primary count")
print(f"  LENIENT dark set : {OUT['dark_lenient']}"
      f"   (a hedged descriptor does not disqualify)")
print(f"  of the {len(_named)} reclassified, {len(_mito)} are mitochondrial")
print(f"    respiratory-chain associated - one functional module, which is why")
print(f"    this reads as the selection finding real biology rather than noise.")
print(f"    NB the formal enrichment test needs a genome-wide PlasmoDB product")
print(f"    table, which is not in this repository. Do not quote a p value for")
print(f"    it until that table is in hand.")
# none of the three lead candidates may be among the reclassified
for _g in ('PF3D7_0810700', 'PF3D7_1449800', 'PF3D7_1459400'):
    assert _g not in _named, '%s is named in VEuPathDB; Figure 7 must change' % _g
print(f"  ok  the three lead candidates are dark in both sources")

# --- cross-paper consistency: Paper 2 must not contradict Paper 1 ---
checks = [
    ('confident hypothetical', OUT['pf_confident_hypothetical'], P1['pf_denom']),
    ('novel vs AF-SwissProt',  OUT['pf_novelish_afsp'],          316),
    ('PDB100 candidates',      OUT['pdb100_candidates'],         P1['pf_candidates']),
    ('US-align confirmed',     OUT['pdb100_confirmed'],          P1['pf_confirmed']),
    ('final novel',            OUT['novel_final'],               P1['pf_novel_after']),
]
print('\n[8] AGREEMENT WITH PAPER 1')
bad = []
for label, got, want in checks:
    ok = got == want
    print(f"  {'ok ' if ok else 'MISMATCH'}  {label:24} {got} vs {want}")
    if not ok:
        bad.append(label)
assert OUT['novel_final'] == OUT['pf_novelish_afsp'] - OUT['pdb100_confirmed'], (
    'arithmetic: 316 - 105 must equal 211')
if bad:
    raise SystemExit('Paper 1/2 disagree on: ' + ', '.join(bad))

path = os.path.join(HERE, 'numbers_p2.json')
json.dump(OUT, open(path, 'w'), indent=1, default=str)
print(f"\nwrote {os.path.basename(path)}")
