# -*- coding: utf-8 -*-
"""Fold class and pocket topology for the 128, and the columns they add to Table 1.

Two artefacts survive the pLDDT filter and have to be caught structurally.

  1. An extended model. AlphaFold predicts long single helices and coiled coils
     with high confidence, so mean pLDDT does not flag them, but a rod has no
     interior. Radius of gyration against the compact-globule expectation
     Rg ~ 2.2 * N**0.38 separates them: a globular domain sits near 1, an
     extended one at 2 and above.

  2. A groove rather than a cavity. fpocket scores the channel running along a
     helix as readily as a real site. A binding cavity is walled by residues
     far apart in sequence; a groove is lined by one contiguous stretch. The
     number of distinct sequence segments contributing the lining separates
     them without any assumption about shape.

Neither test uses the druggability score, so neither is circular.
"""
import os, csv, json, math, re

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, 'fpocket_work')

GLOBULAR_MAX = 1.6      # Rg / Rg_expected
EXTENDED_MIN = 2.0
MIN_SEGMENTS = 3        # distinct sequence stretches lining the pocket
SEG_GAP = 4             # residues; a larger gap starts a new stretch


CORE_PLDDT = 70.0       # a residue belongs to the confident core above this
MIN_CORE = 50           # below this there is no core to measure


def ca_coords(path):
    out = []
    for l in open(path):
        if l.startswith('ATOM') and l[12:16].strip() == 'CA':
            out.append((float(l[30:38]), float(l[38:46]), float(l[46:54]),
                        float(l[60:66])))
    return out


def radius_of_gyration(c):
    n = len(c)
    cx = sum(a[0] for a in c) / n
    cy = sum(a[1] for a in c) / n
    cz = sum(a[2] for a in c) / n
    return math.sqrt(sum((a[0]-cx)**2 + (a[1]-cy)**2 + (a[2]-cz)**2
                         for a in c) / n)


def best_pocket(info):
    """Highest druggability score in an fpocket _info.txt. Pocket 1 is not
    always the best-scoring one, so this searches rather than assumes."""
    cur, best = None, (None, -1.0)
    for line in open(info, errors='replace'):
        line = line.strip()
        m = re.match(r'Pocket\s+(\d+)', line)
        if m:
            cur = int(m.group(1))
            continue
        m = re.match(r'Druggability Score\s*:\s*([\d.eE+-]+)', line)
        if m and cur is not None:
            s = float(m.group(1))
            if s > best[1]:
                best = (cur, s)
    return best


def fold_class(ratio):
    if ratio <= GLOBULAR_MAX:
        return 'globular'
    if ratio >= EXTENDED_MIN:
        return 'extended'
    return 'intermediate'


def analyse(acc):
    pdb = os.path.join(WORK, acc + '.pdb')
    info = os.path.join(WORK, acc + '_out', acc + '_info.txt')
    if not (os.path.exists(pdb) and os.path.exists(info)):
        return None
    c = ca_coords(pdb)
    n = len(c)
    # Measure the confident core, not the whole model. A compact domain with a
    # long disordered tail (PF3D7_1449800 is the clear case) has its Rg
    # inflated by the tail and can be called extended when the fold is not.
    # Two proteins move each way on this; the median over the 128 barely
    # shifts, 2.21 to 2.17, so the headline does not rest on the choice.
    core = [a for a in c if a[3] >= CORE_PLDDT]
    use = core if len(core) >= MIN_CORE else c
    rg_full = radius_of_gyration(c)
    rg = radius_of_gyration(use)
    ratio = rg / (2.2 * len(use) ** 0.38)
    pnum, pscore = best_pocket(info)
    patm = os.path.join(WORK, acc + '_out', 'pockets', 'pocket%s_atm.pdb' % pnum)
    res = sorted({int(l[22:26]) for l in open(patm)
                  if l.startswith(('ATOM', 'HETATM'))}) if os.path.exists(patm) else []
    segs = 1 if res else 0
    for a, b in zip(res, res[1:]):
        if b - a > SEG_GAP:
            segs += 1
    return dict(accession=acc, length=n, core_residues=len(use),
                rg=round(rg, 1), rg_full_model=round(rg_full, 1),
                rg_ratio=round(ratio, 2), fold_class=fold_class(ratio),
                best_pocket=pnum, pocket_score=pscore,
                lining_residues=len(res), lining_segments=segs,
                lining_span=(res[-1] - res[0]) if res else 0,
                mean_pLDDT=round(sum(a[3] for a in c) / n, 1))


def main():
    ids = list(csv.DictReader(open(os.path.join(HERE, 'dark128_ids.tsv'),
                                   encoding='utf-8'), delimiter='\t'))
    rows = []
    for r in ids:
        a = analyse(r['accession'])
        if a is None:
            print('  !! no fpocket output for %s' % r['gene_id'])
            continue
        a['gene_id'] = r['gene_id']
        rows.append(a)

    cols = ['gene_id', 'accession', 'length', 'core_residues', 'mean_pLDDT',
            'rg', 'rg_full_model', 'rg_ratio', 'fold_class', 'best_pocket',
            'pocket_score', 'lining_residues', 'lining_segments', 'lining_span']
    out = os.path.join(HERE, 'dark128_foldclass.tsv')
    with open(out, 'w', newline='', encoding='utf-8') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(cols)
        for r in rows:
            w.writerow([r[c] for c in cols])
    print('wrote dark128_foldclass.tsv  (%d proteins)' % len(rows))

    # ---- fold class and pocket topology onto Table 1 --------------------------
    by_gene = {r['gene_id']: r for r in rows}
    t1 = os.path.join(HERE, 'table1_shortlist.tsv')
    sh = list(csv.DictReader(open(t1, encoding='utf-8'), delimiter='\t'))
    add = ['rg_ratio', 'fold_class', 'lining_segments']
    hdr = [c for c in sh[0].keys() if c not in add] + add
    for r in sh:
        g = by_gene.get(r['gene_id'], {})
        for c in add:
            r[c] = g.get(c, '')
    # rank by the tier a candidate reaches, then by score
    def tier(r):
        if not r['fold_class']:
            return 3
        if r['fold_class'] != 'globular' or int(r['lining_segments']) < MIN_SEGMENTS:
            return 2
        return 0 if float(r['lining_mean_pLDDT']) >= 70 else 1
    sh.sort(key=lambda r: (tier(r), -float(r['pocket_score'])))
    with open(t1, 'w', newline='', encoding='utf-8') as fh:
        w = csv.DictWriter(fh, fieldnames=hdr, delimiter='\t',
                           extrasaction='ignore')
        w.writeheader()
        w.writerows(sh)
    n_pass = sum(1 for r in sh if tier(r) == 0)
    print('updated table1_shortlist.tsv  (%d rows, %d reach the top tier)'
          % (len(sh), n_pass))


if __name__ == '__main__':
    main()
