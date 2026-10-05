#!/bin/bash
# =============================================================================
# Fpocket across the whole dark set — all 128, not a convenience sample.
#
# The existing druggability table covers 36 of 316 proteins and only 6 of the
# 128 that matter, because the earlier driver ran on three hard-coded
# candidates and an older partial run. A shortlist drawn from 5% coverage is
# not a shortlist, so this runs the whole set on the same footing.
#
# Resumable: a protein whose _info.txt already exists is skipped, so an
# interrupted run can simply be restarted.
#
# Output -> paper2_data/fpocket_dark128.tsv
# =============================================================================
set -u
D="$HOME/Downloads"
HERE="$D/paper2_data"
AF="$HOME/apico_structural_analysis/data/structures/pfalciparum"
WORK="$HERE/fpocket_work"
IDS="$HERE/dark128_ids.tsv"
OUT="$HERE/fpocket_dark128.tsv"
mkdir -p "$WORK"
exec > >(tee "$HERE/fpocket_dark128.log") 2>&1

[ -f "$IDS" ] || { echo "missing $IDS"; read -p "Enter"; exit 1; }
[ -d "$AF" ]  || { echo "missing structures at $AF"; read -p "Enter"; exit 1; }

FPOCKET=""
for p in \
  /opt/miniconda3/envs/pfstruct/bin/fpocket \
  /opt/miniconda3/bin/fpocket \
  /usr/local/bin/fpocket \
  /opt/homebrew/bin/fpocket \
  "$(command -v fpocket 2>/dev/null)"; do
  [ -n "${p:-}" ] && [ -x "$p" ] && FPOCKET="$p" && break
done
[ -z "$FPOCKET" ] && { echo "ERROR: fpocket not found. Install it (conda install -c conda-forge fpocket) or edit this script."; read -p "Enter"; exit 1; }
echo "fpocket: $FPOCKET"
echo ""

n=0; ran=0; skipped=0; missing=0
while IFS=$'\t' read -r gene acc; do
  gene="${gene%%$'\r'}"; acc="${acc%%$'\r'}"   # strip CR: a CRLF file
                                              # otherwise appends \r to every path
  [ "$gene" = "gene_id" ] && continue
  [ -z "$acc" ] && continue
  n=$((n+1))
  src="$AF/AF-$acc-F1-model_v6.pdb"
  if [ ! -f "$src" ]; then
    echo "  !! no structure for $gene / $acc"
    missing=$((missing+1)); continue
  fi
  work="$WORK/${acc}.pdb"
  info="$WORK/${acc}_out/${acc}_info.txt"
  if [ -f "$info" ]; then skipped=$((skipped+1)); continue; fi
  cp "$src" "$work"
  "$FPOCKET" -f "$work" >/dev/null 2>&1
  if [ -f "$info" ]; then
    ran=$((ran+1))
    [ $((ran % 20)) -eq 0 ] && echo "  ... $ran run"
  else
    echo "  !! fpocket produced no output for $gene / $acc"
  fi
done < "$IDS"

echo ""
echo "structures: $n   newly run: $ran   already done: $skipped   missing: $missing"

# ---- parse: best pocket per protein -----------------------------------------
python3 - "$WORK" "$IDS" "$OUT" <<'PYEOF'
import os, re, sys, csv
WORK, IDS, OUT = sys.argv[1], sys.argv[2], sys.argv[3]

def best(info):
    """Highest druggability score in an fpocket _info.txt, with its volume."""
    cur, out = {}, []
    for line in open(info, errors='replace'):
        line = line.strip()
        m = re.match(r'Pocket\s+(\d+)', line)
        if m:
            if cur: out.append(cur)
            cur = {}
            continue
        m = re.match(r'Druggability Score\s*:\s*([\d.eE+-]+)', line)
        if m: cur['score'] = float(m.group(1))
        m = re.match(r'Volume\s*:\s*([\d.eE+-]+)', line)
        if m: cur['vol'] = float(m.group(1))
    if cur: out.append(cur)
    out = [p for p in out if 'score' in p]
    if not out: return None, None
    b = max(out, key=lambda p: p['score'])
    return b['score'], b.get('vol')

rows = []
with open(IDS) as fh:
    for r in csv.DictReader(fh, delimiter='\t'):
        info = os.path.join(WORK, f"{r['accession']}_out", f"{r['accession']}_info.txt")
        if os.path.exists(info):
            s, v = best(info)
        else:
            s, v = None, None
        rows.append([r['gene_id'], r['accession'],
                     '' if s is None else round(s, 3),
                     '' if v is None else round(v, 1)])

with open(OUT, 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t')
    w.writerow(['gene_id', 'accession', 'best_druggability_score', 'best_pocket_volume_A3'])
    w.writerows(rows)

scored = [r for r in rows if r[2] != '']
drug = [r for r in scored if float(r[2]) >= 0.5]
print(f"\nparsed {len(scored)} of {len(rows)} with a pocket score")
print(f"  druggable (score >= 0.5): {len(drug)}")
print(f"  wrote {os.path.basename(OUT)}")
PYEOF

echo ""
read -p "Press Enter to close"
