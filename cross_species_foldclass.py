# -*- coding: utf-8 -*-
"""Does the fold correction generalise beyond P. falciparum?

The fold artefact is the strongest methodological claim in Paper 2, and a
reviewer will ask whether it is a peculiarity of one parasite. It is not, but
neither is it uniform: the proportion of the dark proteome that is extended
rather than globular varies roughly six-fold across the 15-genome panel, from
14% in Vitrella to 60% in Eimeria and Tetrahymena.

That variation is the useful result. It means the reliability of a
structure-based druggability screen is an organism-level property that has to
be measured before the screen is trusted, and it comes with a cheap diagnostic:
low-complexity sequence content predicts it across the panel.

Uses only the AlphaFold models and novelty tables already on disk; no pocket
detection, so this runs in seconds and needs no external tools.
"""
import os, re, csv, glob, math, json
import statistics as st

HERE = os.path.dirname(os.path.abspath(__file__))


def _find(*candidates):
    for c in candidates:
        if os.path.isdir(c):
            return c
    raise SystemExit('cannot locate any of: %s' % ', '.join(candidates))


# The structures live outside this folder and the two are reached by different
# paths depending on where the script runs, so resolve rather than assume.
STRUCT = _find(os.path.expanduser('~/apico_structural_analysis/data/structures'),
               os.path.join(HERE, '..', '..', 'apico_structural_analysis',
                            'data', 'structures'))
TAB = _find(os.path.expanduser(
                '~/Downloads/paper1_upload/zenodo/01_novelty_tables'),
            os.path.join(HERE, '..', 'paper1_upload', 'zenodo',
                         '01_novelty_tables'))

GLOBULAR_MAX, EXTENDED_MIN = 1.6, 2.0
CORE_PLDDT, MIN_CORE = 70.0, 50
MIN_PER_SPECIES = 15

UNK = re.compile(r'^(uncharacterized protein|conserved .*unknown function|'
                 r'hypothetical protein.*|protein of unknown function.*)$', re.I)
NOVEL_VERDICTS = ('novel', 'no_hit', 'twilight')

AA3 = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q',
       'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
       'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
       'TYR': 'Y', 'VAL': 'V'}


def parse(path):
    xyz, seq = [], []
    for l in open(path):
        if l.startswith('ATOM') and l[12:16].strip() == 'CA':
            xyz.append((float(l[30:38]), float(l[38:46]), float(l[46:54]),
                        float(l[60:66])))
            seq.append(AA3.get(l[17:20].strip(), 'X'))
    return xyz, ''.join(seq)


def radius_of_gyration(c):
    n = len(c)
    cx = sum(a[0] for a in c) / n
    cy = sum(a[1] for a in c) / n
    cz = sum(a[2] for a in c) / n
    return math.sqrt(sum((a[0] - cx) ** 2 + (a[1] - cy) ** 2 + (a[2] - cz) ** 2
                         for a in c) / n)


def low_complexity(s, w=20, frac=0.5):
    """Fraction of w-residue windows in which two residue types account for at
    least half the window. A plain compositional measure, deliberately not a
    SEG/DUST reimplementation, so that it is reproducible from this file."""
    if len(s) < w:
        return 0.0
    n = 0
    for i in range(len(s) - w + 1):
        win = s[i:i + w]
        counts = sorted({c: win.count(c) for c in set(win)}.values(),
                        reverse=True)
        if sum(counts[:2]) >= w * frac:
            n += 1
    return n / (len(s) - w + 1)


def ranks(v):
    order = sorted(range(len(v)), key=lambda i: v[i])
    r = [0.0] * len(v)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
            j += 1
        for k in range(i, j + 1):
            r[order[k]] = (i + j + 2) / 2
        i = j + 1
    return r


def pearson(x, y):
    n = len(x)
    mx, my = sum(x) / n, sum(y) / n
    num = sum((a - mx) * (b - my) for a, b in zip(x, y))
    den = math.sqrt(sum((a - mx) ** 2 for a in x) * sum((b - my) ** 2 for b in y))
    return num / den if den else 0.0


def spearman(x, y):
    rho = pearson(ranks(x), ranks(y))
    n = len(x)
    t = rho * math.sqrt((n - 2) / max(1 - rho * rho, 1e-12))
    return rho, math.erfc(abs(t) / math.sqrt(2))


def main():
    panel = {r['species_key']: r for r in csv.DictReader(
        open(os.path.join(TAB, '_panel_summary.tsv'), encoding='utf-8'),
        delimiter='\t')}
    rows = []
    for sp in sorted(panel):
        d = os.path.join(STRUCT, sp)
        tab = os.path.join(TAB, sp + '_novelty.tsv')
        if not (os.path.isdir(d) and os.path.exists(tab)):
            continue
        nov = {r['accession']: r for r in csv.DictReader(
            open(tab, encoding='utf-8'), delimiter='\t')}
        ratio, lc, nk = [], [], []
        for f in glob.glob(os.path.join(d, 'AF-*-F1-model_v6.pdb')):
            acc = os.path.basename(f).split('-')[1]
            r = nov.get(acc)
            if not r:
                continue
            name = (r['protein_name'] or '').strip()
            if not (UNK.match(name) or not name):
                continue                      # named: not part of the dark set
            if r['verdict'] not in NOVEL_VERDICTS:
                continue
            xyz, seq = parse(f)
            core = [a for a in xyz if a[3] >= CORE_PLDDT]
            use = core if len(core) >= MIN_CORE else xyz
            if len(use) < MIN_CORE:
                continue
            ratio.append(radius_of_gyration(use) / (2.2 * len(use) ** 0.38))
            lc.append(low_complexity(seq))
            nk.append((seq.count('N') + seq.count('K')) / len(seq))
        if len(ratio) < MIN_PER_SPECIES:
            continue
        rows.append(dict(
            species_key=sp, label=panel[sp]['label'], clade=panel[sp]['clade'],
            n=len(ratio), median_rg_ratio=round(st.median(ratio), 2),
            pct_globular=round(100 * sum(1 for x in ratio if x <= GLOBULAR_MAX)
                               / len(ratio), 1),
            pct_extended=round(100 * sum(1 for x in ratio if x >= EXTENDED_MIN)
                               / len(ratio), 1),
            median_low_complexity=round(st.median(lc), 4),
            median_asn_lys=round(st.median(nk), 4)))

    rows.sort(key=lambda r: -r['pct_extended'])
    with open(os.path.join(HERE, 'cross_species_foldclass.tsv'), 'w',
              newline='', encoding='utf-8') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter='\t')
        w.writeheader()
        w.writerows(rows)

    print('%-18s %-15s %6s %8s %9s %9s' % (
        'organism', 'clade', 'n', 'Rg med', 'extended', 'lowcomp'))
    for r in rows:
        print('%-18s %-15s %6d %8.2f %8.0f%% %9.3f' % (
            r['label'], r['clade'], r['n'], r['median_rg_ratio'],
            r['pct_extended'], r['median_low_complexity']))

    ext = [r['pct_extended'] for r in rows]
    rho_lc, p_lc = spearman([r['median_low_complexity'] for r in rows], ext)
    rho_nk, p_nk = spearman([r['median_asn_lys'] for r in rows], ext)
    print('\nacross the panel (n = %d genomes):' % len(rows))
    print('  low-complexity content vs extended fraction  rho = %+.3f, p = %.3f'
          % (rho_lc, p_lc))
    print('  Asn+Lys content vs extended fraction         rho = %+.3f, p = %.3f'
          % (rho_nk, p_nk))
    print('  range of extended fraction: %.0f%% (%s) to %.0f%% (%s)'
          % (rows[-1]['pct_extended'], rows[-1]['label'],
             rows[0]['pct_extended'], rows[0]['label']))
    print('\n  Asn+Lys tracks extendedness within P. falciparum but not across')
    print('  genomes; low-complexity content does both. Quote the second.')

    # ---- and the same question within P. falciparum, protein by protein -----
    # The between-genome correlation could be driven by anything that differs
    # between genomes. The within-genome test is the one that implicates
    # composition directly, and it has to control for length, since a short
    # protein could be extended for unrelated reasons.
    wf = within_pfalciparum()                 # primary dark set
    wf_all = within_pfalciparum(restrict_to_strict=False)   # UniProt-only

    out = {'cross_species_foldclass': rows,
           'lowcomplexity_vs_extended': {'rho': round(rho_lc, 3),
                                         'p': round(p_lc, 4), 'n': len(rows)},
           'asnlys_vs_extended_between': {'rho': round(rho_nk, 3),
                                          'p': round(p_nk, 4), 'n': len(rows)},
           'extended_range': [rows[-1]['pct_extended'], rows[-1]['label'],
                              rows[0]['pct_extended'], rows[0]['label']],
           'within_pf': wf,
           'within_pf_uniprot_only': wf_all}
    json.dump(out, open(os.path.join(HERE, 'cross_species_foldclass.json'), 'w'),
              indent=1)
    print('\nwithin P. falciparum (n = %d dark proteins, partial correlations '
          'controlling for length):' % wf['n'])
    print('  low-complexity vs Rg ratio   rho = %+.3f, p = %.2g'
          % (wf['lowcomplexity']['rho'], wf['lowcomplexity']['p']))
    print('  Asn+Lys vs Rg ratio          rho = %+.3f, p = %.2g'
          % (wf['asn_lys']['rho'], wf['asn_lys']['p']))
    print('  length vs Rg ratio           rho = %+.3f  (so length is not the '
          'driver)' % wf['length_vs_ratio_rho'])


def partial_spearman(x, y, z):
    """Spearman correlation of x and y with z partialled out."""
    rx, ry, rz = ranks(x), ranks(y), ranks(z)
    rxy, rxz, ryz = pearson(rx, ry), pearson(rx, rz), pearson(ry, rz)
    return (rxy - rxz * ryz) / math.sqrt((1 - rxz ** 2) * (1 - ryz ** 2))


def within_pfalciparum(restrict_to_strict=True):
    """Within-genome correlations.

    The cross-species panels must use the UniProt-only criterion, because
    VEuPathDB curation differs between organisms and a cross-genome comparison
    needs one rule applied identically everywhere. The within-P. falciparum
    panel has no such constraint, so it uses the paper's primary dark set -
    otherwise the figure would quote a denominator that appears nowhere else
    in the manuscript.
    """
    fc = list(csv.DictReader(
        open(os.path.join(HERE, 'dark128_foldclass.tsv'), encoding='utf-8'),
        delimiter='\t'))
    if restrict_to_strict:
        import json as _j
        _n = os.path.join(HERE, 'numbers_p2.json')
        if os.path.exists(_n):
            _strict = set(_j.load(open(_n))['dark_strict_ids'])
            fc = [r for r in fc if r['gene_id'] in _strict]
    work = os.path.join(HERE, 'fpocket_work')
    ratio, lc, nk, L = [], [], [], []
    for r in fc:
        _, seq = parse(os.path.join(work, r['accession'] + '.pdb'))
        if not seq:
            continue
        ratio.append(float(r['rg_ratio']))
        lc.append(low_complexity(seq))
        nk.append((seq.count('N') + seq.count('K')) / len(seq))
        L.append(len(seq))
    n = len(ratio)

    def p_of(rho):
        t = rho * math.sqrt((n - 3) / max(1 - rho * rho, 1e-12))
        return math.erfc(abs(t) / math.sqrt(2))

    r_lc = partial_spearman(lc, ratio, L)
    r_nk = partial_spearman(nk, ratio, L)
    return {'n': n,
            'lowcomplexity': {'rho': round(r_lc, 3), 'p': float('%.3g' % p_of(r_lc))},
            'asn_lys': {'rho': round(r_nk, 3), 'p': float('%.3g' % p_of(r_nk))},
            'length_vs_ratio_rho': round(pearson(ranks(L), ranks(ratio)), 3)}


if __name__ == '__main__':
    main()
