#!/usr/bin/env python3
"""Normalise the raw MEGA exports into analysis-ready files under data/.

Produces:
  data/alignment/metalloproteome_aln.fasta   cleaned protein alignment
  data/tree/mega_tree.nwk                    Newick with bootstrap as node labels
  data/metadata/tips.csv                     one row per tip, ready for annotation
"""
import re, csv, os

RAW_MEG = "data/raw/Divyam_tree.meg"
RAW_NWK = "data/raw/mega_tree_raw.nwk"

# ---------------------------------------------------------------- alignment
def read_meg(path):
    seqs, order, cur = {}, [], None
    for line in open(path):
        line = line.rstrip("\n").strip()
        if not line or line.startswith("!") or line.lower().startswith("#mega"):
            continue
        if line.startswith("#"):
            cur = line[1:].strip()
            if cur not in seqs:
                seqs[cur] = []
                order.append(cur)
        elif cur is not None:
            seqs[cur].append(re.sub(r"\s+", "", line))
    return {k: "".join(v) for k, v in seqs.items()}, order

def clean_name(n):
    """Strip the embedded UniProt description some records carry, and
    normalise MEGA's 'ACC YERPE(2)' spelling back to 'ACC_YERPE(2)'."""
    n = n.split('_"')[0]
    n = re.sub(r"\s+YERPE", "_YERPE", n)
    n = re.sub(r"\s+", "_", n)
    # MEGA marks repeated accessions as "(2)"/"(3)"; parentheses are structural
    # characters in Newick, so recode them as a plain suffix.
    n = re.sub(r"\((\d+)\)$", r"_copy\1", n)
    return n

seqs, order = read_meg(RAW_MEG)
dropped = [n for n in order if set(seqs[n]) <= {"-", "?"}]
kept = [n for n in order if n not in dropped]

L = {len(seqs[n]) for n in kept}
assert len(L) == 1, f"ragged alignment: {L}"

os.makedirs("data/alignment", exist_ok=True)
with open("data/alignment/metalloproteome_aln.fasta", "w") as fh:
    for n in kept:
        fh.write(f">{clean_name(n)}\n")
        s = seqs[n]
        for i in range(0, len(s), 60):
            fh.write(s[i:i+60] + "\n")

# ---------------------------------------------------------------- tree
nwk = open(RAW_NWK).read().strip()
# MEGA writes support as  <prop>[stddev=<sd>]  on internal nodes.
stddev = re.findall(r"\[stddev=([0-9.]+)\]", nwk)
nwk_clean = re.sub(r"\[stddev=[0-9.eE+-]+\]", "", nwk)
# bootstrap proportions 0-1 -> percentages, so the figure legend reads in %
def pct(m):
    return f"{round(float(m.group(1)) * 100):d}"
nwk_pct = re.sub(r"\)([01]\.\d+)", lambda m: ")" + pct(m), nwk_clean)
# normalise quoted tip labels ('ACC YERPE(2)') to the unquoted canonical spelling
def _unquote(m):
    return clean_name(m.group(1))
nwk_pct = re.sub(r"'([^']+)'", _unquote, nwk_pct)
if not nwk_pct.endswith(";"):
    nwk_pct += ";"
os.makedirs("data/tree", exist_ok=True)
open("data/tree/mega_tree.nwk", "w").write(nwk_pct + "\n")

# ---------------------------------------------------------------- tip table
tips = [(a or b).strip() for a, b in
        re.findall(r"[(,]\s*(?:'([^']+)'|([^(),:]+))\s*:", nwk_clean)]
tips = [clean_name(t) for t in tips if t]

rows = []
for t in tips:
    m = re.match(r"(sp|tr)\|([A-Z0-9]+)\|(\S+?)_YERPE(?:_copy(\d+))?$", t)
    if not m:
        rows.append({"tip_label": t, "db": "", "accession": "", "entry_name": "",
                     "copy": 1, "reviewed": ""})
        continue
    db, acc, entry, copy = m.groups()
    rows.append({
        "tip_label": t,
        "db": db,
        "accession": acc,
        "entry_name": entry,
        "copy": int(copy) if copy else 1,
        "reviewed": "yes" if db == "sp" else "no",
    })

os.makedirs("data/metadata", exist_ok=True)
cols = ["tip_label", "db", "accession", "entry_name", "copy", "reviewed",
        "gene_name", "protein_name", "display_label", "clade"]
with open("data/metadata/tips.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols)
    w.writeheader()
    for r in rows:
        r.setdefault("gene_name", "")
        r.setdefault("protein_name", "")
        r.setdefault("display_label", "")
        r.setdefault("clade", "")
        w.writerow(r)

print(f"alignment : {len(kept)} sequences x {L.pop()} columns  (dropped: {dropped})")
print(f"tree      : {len(tips)} tips, {len(stddev)} internal nodes with support")
print(f"metadata  : data/metadata/tips.csv  ({len(rows)} rows, "
      f"{sum(1 for r in rows if r['reviewed']=='yes')} reviewed)")
