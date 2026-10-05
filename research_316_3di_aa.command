#!/bin/bash
# =============================================================================
# Re-search the 316 in Foldseek's default 3Di+AA mode.
#
# WHY THIS EXISTS
# The forward search that produced the 316 was run with --alignment-type 0,
# which is the 3Di alphabet alone. Foldseek's own documentation calls that mode
# "not recommended"; the default and more sensitive local mode is 2 (3Di+AA).
# The manuscript Methods previously described the search as combined 3Di and
# amino-acid while citing mode 0, which was wrong in both papers.
#
# A less sensitive entry search can only ADD proteins to a novel set, never
# remove them, so the 316 is an upper bound. This script measures by how much,
# using the same machinery as the rest of the paper so the number is
# reproducible from this repository rather than quoted from elsewhere.
#
# WHAT IT DOES
#   1. re-searches the 316 against AlphaFold/Swiss-Prot in mode 2 (3Di+AA)
#   2. reports qtmscore, excludes P. falciparum self-matches
#   3. counts how many now reach a known fold at qTM >= 0.50
#   4. splits those into ones the PDB100 step already removed, and ones that
#      would newly leave the final set
#   5. checks every protein in table 1 individually
#
# Output -> research_316_3di_aa.tsv  and  research_316_summary.txt
# =============================================================================
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
AF="$HOME/apico_structural_analysis/data/structures/pfalciparum"
WORK="$HERE/research316_work"
OUT="$HERE/research_316_3di_aa.tsv"
SUM="$HERE/research_316_summary.txt"
mkdir -p "$WORK"

FOLDSEEK=""
for p in /opt/miniconda3/envs/pfstruct/bin/foldseek /opt/miniconda3/bin/foldseek \
         /usr/local/bin/foldseek /opt/homebrew/bin/foldseek \
         "$(command -v foldseek 2>/dev/null)"; do
  [ -n "${p:-}" ] && [ -x "$p" ] && FOLDSEEK="$p" && break
done
[ -z "$FOLDSEEK" ] && { echo "foldseek not found"; read -p "Enter"; exit 1; }
echo "foldseek: $FOLDSEEK"
"$FOLDSEEK" version 2>/dev/null | head -1

# ---- the 316 accessions -----------------------------------------------------
ACC="$WORK/acc_316.txt"
python3 - "$HERE" "$ACC" <<'PYEOF'
import sys, os, csv, re
here, out = sys.argv[1], sys.argv[2]
nov = os.path.join(here, '..', 'paper1_upload', 'zenodo',
                   '01_novelty_tables', 'pfalciparum_novelty.tsv')
keep = []
for r in csv.DictReader(open(nov, encoding='utf-8'), delimiter='\t'):
    if r['hypothetical'] == '1' and r['verdict'] in ('novel', 'no_hit', 'twilight'):
        keep.append(r['accession'])
open(out, 'w').write('\n'.join(keep) + '\n')
print('  %d accessions written' % len(keep))
PYEOF

# ---- stage the structures ---------------------------------------------------
Q="$WORK/queries"
rm -rf "$Q"; mkdir -p "$Q"
n=0
while read -r acc; do
  [ -z "$acc" ] && continue
  src="$AF/AF-$acc-F1-model_v6.pdb"
  [ -f "$src" ] && cp "$src" "$Q/" && n=$((n+1))
done < "$ACC"
echo "  staged $n structures"

# ---- the AlphaFold/Swiss-Prot target database -------------------------------
DB="$WORK/afsp"
if [ ! -f "${DB}.dbtype" ]; then
  echo "downloading Alphafold/Swiss-Prot (one time, several GB)..."
  "$FOLDSEEK" databases Alphafold/Swiss-Prot "$DB" "$WORK/tmp" || {
    echo "database download failed"; read -p "Enter"; exit 1; }
fi

# ---- the search: mode 2, the default, 3Di+AA --------------------------------
echo "searching in 3Di+AA mode (--alignment-type 2)..."
"$FOLDSEEK" easy-search "$Q" "$DB" "$OUT" "$WORK/tmp" \
  --alignment-type 2 -e 0.001 --threads 4 \
  --format-output "query,target,qtmscore,alntmscore,evalue,qcov,taxname" \
  || { echo "search failed"; read -p "Enter"; exit 1; }

# ---- interpret --------------------------------------------------------------
python3 - "$HERE" "$OUT" "$SUM" <<'PYEOF'
import sys, os, csv, re
here, res, summ = sys.argv[1], sys.argv[2], sys.argv[3]

best = {}
for line in open(res):
    f = line.rstrip('\n').split('\t')
    if len(f) < 7:
        continue
    q = f[0].split('-')[1] if f[0].startswith('AF-') else f[0]
    tax = f[6]
    # a Plasmodium falciparum target is a self-match and is excluded, exactly
    # as the PDB100 step excludes own-genus entries
    if 'falciparum' in tax.lower():
        continue
    try:
        qtm = float(f[2])
    except ValueError:
        continue
    if q not in best or qtm > best[q][0]:
        best[q] = (qtm, f[1], tax)

hit = {q: v for q, v in best.items() if v[0] >= 0.50}

ver = os.path.join(here, '..', 'pdb100_verification')
removed = set(open(os.path.join(ver, 'removed_105_accessions.txt')).read().split())
final211 = set(open(os.path.join(ver, 'novel_211_accessions.txt')).read().split())

already = {q for q in hit if q in removed}
newly = {q for q in hit if q in final211}

ids = {r['gene_id']: r['accession'] for r in csv.DictReader(
    open(os.path.join(here, 'dark128_ids.tsv'), encoding='utf-8'), delimiter='\t')}
t1 = [r['gene_id'] for r in csv.DictReader(
    open(os.path.join(here, 'table1_shortlist.tsv'), encoding='utf-8'), delimiter='\t')]

L = []
L.append('Re-search of the 316 in 3Di+AA mode (--alignment-type 2), qtmscore,')
L.append('P. falciparum self-matches excluded.')
L.append('')
L.append('  reach qTM >= 0.50 against AF/Swiss-Prot : %d of 316' % len(hit))
L.append('    already removed by the PDB100 step    : %d' % len(already))
L.append('    still inside the final 211            : %d' % len(newly))
L.append('')
L.append('  so the counts would move: 316 -> %d, 211 -> %d'
         % (316 - len(hit), 211 - len(newly)))
L.append('')
if newly:
    L.append('  the proteins that would newly leave the novel set:')
    for q in sorted(newly):
        qtm, tgt, tax = best[q]
        L.append('    %-12s qTM %.3f  %s  (%s)' % (q, qtm, tgt, tax[:40]))
# the question that decides whether anything downstream has to change
d128 = {r['accession']: r['gene_id'] for r in csv.DictReader(
    open(os.path.join(here, 'dark128_ids.tsv'), encoding='utf-8'), delimiter='\t')}
import json as _json
_N = _json.load(open(os.path.join(here, 'numbers_p2.json')))
strict = set(_N['dark_strict_ids'])
in128 = [q for q in newly if q in d128]
in99 = [q for q in in128 if d128[q] in strict]
L.append('')
L.append('  impact on the dark set:')
L.append('    in the 128 (UniProt-only)  : %d  -> 128 becomes %d'
         % (len(in128), 128 - len(in128)))
L.append('    in the strict 99 (PRIMARY) : %d  -> 99 becomes %d'
         % (len(in99), 99 - len(in99)))
for q in in128:
    L.append('      %s (%s), strict: %s' % (q, d128[q], d128[q] in strict))
L.append('')
L.append('  table 1 candidates:')
bad = []
for g in t1:
    a = ids.get(g)
    if a in hit:
        bad.append(g)
        L.append('    %-16s %-12s NOW MATCHES qTM %.3f  <-- affects the paper'
                 % (g, a, best[a][0]))
    else:
        b = best.get(a)
        L.append('    %-16s %-12s still novel (best qTM %s)'
                 % (g, a, ('%.3f' % b[0]) if b else 'no hit'))
L.append('')
L.append('  VERDICT: %s' % ('TABLE 1 AFFECTED - do not submit until resolved'
                            if bad else 'no table 1 candidate is affected'))
open(summ, 'w').write('\n'.join(L) + '\n')
print('\n'.join(L))
PYEOF

echo ""
echo "wrote $(basename "$SUM")"
read -p "Press Enter to close"
