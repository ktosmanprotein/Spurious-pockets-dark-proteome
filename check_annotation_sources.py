# -*- coding: utf-8 -*-
"""Is the dark set actually dark? UniProt against VEuPathDB, for all 128.

The dark set was defined from UniProt protein names. VEuPathDB curates the
same genes independently and more recently, and it disagrees: several of the
128 carry a specific functional descriptor there.

The trap is the word "putative", which is what inflated the original 1,017.
It attaches to two completely different things:

    "conserved Plasmodium protein, unknown function"   - says nothing
    "cytochrome b-c1 complex subunit 8, putative"      - says a great deal

so the test cannot be whether "putative" appears. It is whether anything
informative remains once "putative" is stripped. A descriptor naming a complex,
a domain, a family or a gene product is informative; "conserved", "hypothetical"
and "unknown function" are not, however they are qualified.

Two counts come out of this, and the paper should report both:

  strict   a specific descriptor in either source disqualifies a protein from
           the dark set. This is what a reviewer checking VEuPathDB will apply.
  lenient  a "putative" VEuPathDB descriptor is an unverified inference, often
           transferred from an ortholog - the very kind of annotation this
           paper argues is unreliable - so it does not disqualify.

Writes dark128_annotation_check.tsv and prints the reconciliation.
"""
import os, re, csv, collections

HERE = os.path.dirname(os.path.abspath(__file__))

# says nothing, in any combination
EMPTY = re.compile(
    r'^(conserved\s+)?(plasmodium\s+)?(exported\s+)?(membrane\s+)?'
    r'protein(\s*\(hyp\d+\))?(,)?\s*(of\s+)?unknown\s+function$'
    r'|^hypothetical\s+protein$'
    r'|^uncharacteri[sz]ed\s+protein$'
    r'|^conserved\s+protein$'
    r'|^protein\s+of\s+unknown\s+function$'
    r'|^$', re.I)


def strip_putative(name):
    """Remove the hedge, not the content."""
    s = (name or '').strip()
    s = re.sub(r'\s*,?\s*putative\s*$', '', s, flags=re.I)
    s = re.sub(r'\s*\|\s*domain-containing protein$', ' domain-containing protein',
               s, flags=re.I)
    return s.strip().rstrip(',').strip()


def is_informative(name):
    core = strip_putative(name)
    return not EMPTY.match(core)


def was_putative(name):
    return bool(re.search(r'putative\s*$', (name or '').strip(), re.I))


def tsv(path):
    return list(csv.DictReader(open(path, encoding='utf-8'), delimiter='\t'))


def main():
    fold = tsv(os.path.join(HERE, 'dark128_foldclass.tsv'))
    dark = [r['gene_id'] for r in fold]
    ann = {r['gene_id']: r for r in
           tsv(os.path.join(HERE, 'paper2_novel211_annotated.tsv'))}

    # UniProt names, from the Paper 1 novelty table via the accession mapping
    nov = tsv(os.path.join(HERE, '..', 'paper1_upload', 'zenodo',
                           '01_novelty_tables', 'pfalciparum_novelty.tsv'))
    uni = {}
    for r in nov:
        m = re.search(r'(PF3D7_\w+)', r['veupathdb'])
        if m:
            uni[m.group(1)] = (r['protein_name'] or '').strip()

    rows, flagged = [], []
    for g in dark:
        vp = (ann.get(g, {}).get('product') or '').strip()
        up = uni.get(g, '')
        vi, ui = is_informative(vp), is_informative(up)
        verdict = ('dark' if not (vi or ui)
                   else 'named in VEuPathDB only' if vi and not ui
                   else 'named in UniProt only' if ui and not vi
                   else 'named in both')
        rows.append(dict(gene_id=g, uniprot_name=up, veupathdb_product=vp,
                         uniprot_informative=int(ui),
                         veupathdb_informative=int(vi),
                         veupathdb_putative=int(was_putative(vp)),
                         verdict=verdict))
        if vi or ui:
            flagged.append(rows[-1])

    with open(os.path.join(HERE, 'dark128_annotation_check.tsv'), 'w',
              newline='', encoding='utf-8') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter='\t')
        w.writeheader()
        w.writerows(rows)

    n = len(rows)
    n_dark = sum(1 for r in rows if r['verdict'] == 'dark')
    n_put = sum(1 for r in flagged if r['veupathdb_putative'])
    print('cross-source annotation check, %d proteins' % n)
    print('  dark in both sources          : %d' % n_dark)
    for k, v in collections.Counter(
            r['verdict'] for r in rows if r['verdict'] != 'dark').most_common():
        print('  %-30s: %d' % (k, v))
    print()
    print('  STRICT  dark set = %d  (any specific descriptor disqualifies)' % n_dark)
    print('  LENIENT dark set = %d  (a "putative" VEuPathDB descriptor does not)'
          % (n_dark + n_put))
    print('          %d of the %d newly named carry "putative"' % (n_put, len(flagged)))
    print()
    print('  the %d with a descriptor VEuPathDB does not hedge:' % (len(flagged) - n_put))
    for r in flagged:
        if not r['veupathdb_putative']:
            print('    %-16s %s' % (r['gene_id'], r['veupathdb_product']))
    print()
    print('  hedged as putative:')
    for r in flagged:
        if r['veupathdb_putative']:
            print('    %-16s %s' % (r['gene_id'], r['veupathdb_product']))


if __name__ == '__main__':
    main()
