# -*- coding: utf-8 -*-
"""Step C: GWAS Catalog ontology-annotated associations.
Phases:
  download : fetch the ontology-annotated full zip, record HTTP metadata, unzip.
  traits   : list all MAPPED_TRAIT containing 'Alzheimer' (+URI+count) -> alzheimer_trait_list.csv
  match    : filter AD rows (strict A-level vs broad), proxy flag, compare vs our 50 SNPs
             by (i) gene and (ii) +/-5kb coordinate, write gwas_ad_matches.csv.
Run:  python fetch_gwas_catalog.py download
      python fetch_gwas_catalog.py traits
      python fetch_gwas_catalog.py match
"""
import os
import sys
import csv
import json
import zipfile
import urllib.request
import datetime

import config as C

STRICT_TRAITS_FILE = os.path.join(C.GWAS_DIR, "strict_traits.txt")

def download():
    os.makedirs(C.GWAS_DIR, exist_ok=True)
    zippath = os.path.join(C.GWAS_DIR, "gwas-catalog-associations_ontology-annotated-full.zip")
    metapath = os.path.join(C.GWAS_DIR, "download_meta.txt")
    if os.path.exists(zippath) and os.path.getsize(zippath) > 0:
        print("zip already present:", zippath)
    else:
        print("downloading", C.GWAS_ZIP_URL)
        req = urllib.request.Request(
            C.GWAS_ZIP_URL,
            headers={"User-Agent": "50snp-annotation/1.0 (contact: local)"},
        )
        with urllib.request.urlopen(req, timeout=300) as r:
            final_url = r.geturl()
            lastmod = r.headers.get("Last-Modified", "")
            n = 0
            with open(zippath, "wb") as out:
                while True:
                    chunk = r.read(1 << 20)
                    if not chunk:
                        break
                    out.write(chunk)
                    n += len(chunk)
            print("downloaded", n, "bytes")
        with open(metapath, "w", encoding="utf-8") as f:
            f.write(f"requested_url={C.GWAS_ZIP_URL}\n")
            f.write(f"final_url={final_url}\n")
            f.write(f"last_modified={lastmod}\n")
            f.write(f"download_time={datetime.datetime.now().isoformat(timespec='seconds')}\n")
        print("meta ->", metapath)
        print("final_url =", final_url)
        print("last_modified =", lastmod)

    tsv_path = None
    with zipfile.ZipFile(zippath) as z:
        names = [n for n in z.namelist() if n.lower().endswith(".tsv")]
        if not names:
            names = z.namelist()
        for n in names:
            if n.lower().endswith(".tsv"):
                tsv_path = os.path.join(C.GWAS_DIR, os.path.basename(n))
                if not os.path.exists(tsv_path):
                    print("extracting", n)
                    z.extract(n, C.GWAS_DIR)
                else:
                    print("tsv already extracted:", tsv_path)
                break
    if tsv_path is None:
        raise SystemExit("no .tsv found in zip")

    return tsv_path

def _tsv_path():
    for n in os.listdir(C.GWAS_DIR):
        if n.lower().endswith(".tsv"):
            return os.path.join(C.GWAS_DIR, n)
    raise SystemExit("no unzipped .tsv in " + C.GWAS_DIR + " — run 'download' first")

def _read_header(path):
    with open(path, encoding="utf-8", newline="") as f:
        header = f.readline().rstrip("\n")
    if header.startswith("﻿"):
        header = header[1:]
    cols = header.split("\t")
    idx = {c.strip(): i for i, c in enumerate(cols)}
    return cols, idx

REQUIRED = ["MAPPED_TRAIT", "MAPPED_TRAIT_URI", "DISEASE/TRAIT", "P-VALUE",
            "PUBMEDID", "CHR_ID", "CHR_POS", "REPORTED GENE(S)", "MAPPED_GENE",
            "STRONGEST SNP-RISK ALLELE", "SNPS"]

def _check_columns(idx):
    missing = [c for c in REQUIRED if c not in idx]
    if missing:
        print("MISSING columns:", missing)
        print("available columns:", sorted(idx))
        raise SystemExit(1)

def traits():
    path = _tsv_path()
    _, idx = _read_header(path)
    _check_columns(idx)
    mi = idx["MAPPED_TRAIT"]
    ui = idx["MAPPED_TRAIT_URI"]
    counts = {}
    with open(path, encoding="utf-8", newline="") as f:
        f.readline()
        for line in f:
            cols = line.rstrip("\n").split("\t")
            if len(cols) <= mi:
                continue
            trait = cols[mi].strip()
            if "alzheimer" in trait.lower():
                uri = cols[ui].strip() if len(cols) > ui else ""
                key = (trait, uri)
                counts[key] = counts.get(key, 0) + 1
    rows = sorted(counts.items(), key=lambda kv: -kv[1])
    out = os.path.join(C.GWAS_DIR, "alzheimer_trait_list.csv")
    with open(out, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["MAPPED_TRAIT", "MAPPED_TRAIT_URI", "count"])
        for (trait, uri), c in rows:
            w.writerow([trait, uri, c])
    print(f"\nAlzheimer MAPPED_TRAIT list ({len(rows)} distinct) -> {out}\n")
    for (trait, uri), c in rows:
        print(f"{c:>6}  {trait}  [{uri}]")
    return rows

def load_strict_traits():
    if not os.path.exists(STRICT_TRAITS_FILE):
        raise SystemExit(f"strict_traits.txt not found at {STRICT_TRAITS_FILE}; "
                         "create it (one MAPPED_TRAIT per line) after reviewing the trait list.")
    with open(STRICT_TRAITS_FILE, encoding="utf-8") as f:
        return {ln.strip() for ln in f if ln.strip() and not ln.startswith("#")}

def _split_chrpos(text):
    """Split CHR_POS on ';' and ' x '; return (parsed_ints, unparseable_list)."""
    nums, bad = [], []
    if not text:
        return nums, bad
    for seg in re_split_multi(text):
        s = seg.strip()
        if not s:
            continue
        try:
            nums.append(int(s))
        except ValueError:
            bad.append(s)
    return nums, bad

def re_split_multi(text):
    import re
    return re.split(r"\s*;\s*|\s+x\s+", text)

def match():
    strict = load_strict_traits()
    path = _tsv_path()
    _, idx = _read_header(path)
    _check_columns(idx)
    mi = idx["MAPPED_TRAIT"]; ui = idx["MAPPED_TRAIT_URI"]; di = idx["DISEASE/TRAIT"]
    pi = idx["P-VALUE"]; pmid = idx["PUBMEDID"]; ci = idx["CHR_ID"]; cpi = idx["CHR_POS"]
    rgi = idx["REPORTED GENE(S)"]; mgi = idx["MAPPED_GENE"]; sai = idx["STRONGEST SNP-RISK ALLELE"]
    sni = idx["SNPS"]

    locus = C.load_locus()
    snp_info = {}
    for rsid, lr in locus.items():
        snp_info[rsid] = {
            "chr": lr["chr"], "pos": lr["pos"], "locus": lr["locus"],
            "genes": C.our_gene_set(lr),
        }

    filtered_out = os.path.join(C.GWAS_DIR, "gwas_ad_filtered.csv")
    matches_out = C.OUT_GWAS_MATCHES

    n_skip_p = 0
    bad_chrpos = []
    filter_fields = ["source_row", "PUBMEDID", "DISEASE/TRAIT", "MAPPED_TRAIT",
                     "MAPPED_TRAIT_URI", "P-VALUE", "tier", "proxy_flag"]

    match_fields = ["our_rsID", "our_locus", "our_chr", "our_pos", "PUBMEDID",
                    "STRONGEST SNP-RISK ALLELE", "P-VALUE", "DISEASE/TRAIT",
                    "MAPPED_TRAIT", "tier", "proxy_flag", "match_type",
                    "gene_source", "source_row"]

    with open(path, encoding="utf-8", newline="") as f, \
         open(filtered_out, "w", encoding="utf-8", newline="") as ff, \
         open(matches_out, "w", encoding="utf-8", newline="") as fm:
        f.readline()
        fw = csv.DictWriter(ff, fieldnames=filter_fields)
        fw.writeheader()
        mw = csv.DictWriter(fm, fieldnames=match_fields)
        mw.writeheader()
        for lineno, line in enumerate(f, start=2):
            cols = line.rstrip("\n").split("\t")
            def g(i):
                return cols[i].strip() if i < len(cols) else ""
            trait = g(mi)
            if "alzheimer" not in trait.lower():
                continue
            pv_text = g(pi)
            try:
                pv = float(pv_text)
            except ValueError:
                n_skip_p += 1
                continue
            if pv >= C.P_THRESH:
                continue
            disease = g(di)
            tier = "A" if trait in strict else "broad"
            proxy = "yes" if ("family history" in (trait + " " + disease).lower()
                              or "proxy" in (trait + " " + disease).lower()) else "no"

            row = {
                "source_row": lineno, "PUBMEDID": g(pmid), "DISEASE/TRAIT": disease,
                "MAPPED_TRAIT": trait, "MAPPED_TRAIT_URI": g(ui),
                "P-VALUE": pv_text, "tier": tier, "proxy_flag": proxy,
            }
            fw.writerow(row)

            gwas_genes = []
            for gt in (g(rgi), g(mgi)):
                for gene, src in C.split_genes_gwas(gt):
                    gwas_genes.append((gene, src))
            gwas_gene_set = {g for g, _ in gwas_genes}

            chr_id = g(ci)
            pos_ints, pos_bad = _split_chrpos(g(cpi))
            for b in pos_bad:
                bad_chrpos.append((lineno, b))

            for rsid, si in snp_info.items():
                hit = []
                gsrc = set()
                if gwas_gene_set and (si["genes"] & gwas_gene_set):
                    hit.append("i")
                    gsrc.update(src for g, src in gwas_genes if g in si["genes"])
                lo = si["pos"] - C.WINDOW_BP
                hi = si["pos"] + C.WINDOW_BP
                if chr_id.strip().upper() == si["chr"].upper() and any(lo <= p <= hi for p in pos_ints):
                    hit.append("ii")
                if hit:
                    mw.writerow({
                        "our_rsID": rsid, "our_locus": si["locus"],
                        "our_chr": si["chr"], "our_pos": si["pos"],
                        "PUBMEDID": g(pmid), "STRONGEST SNP-RISK ALLELE": g(sai),
                        "P-VALUE": pv_text, "DISEASE/TRAIT": disease,
                        "MAPPED_TRAIT": trait, "tier": tier, "proxy_flag": proxy,
                        "match_type": ",".join(hit),
                        "gene_source": ",".join(sorted(gsrc)) if gsrc else "",
                        "source_row": lineno,
                    })

    print(f"\nfiltered AD rows -> {filtered_out}")
    print(f"matches -> {matches_out}")
    print(f"Alzheimer rows skipped for non-numeric P-VALUE: {n_skip_p}")
    print(f"unparseable CHR_POS values (row, value): {len(bad_chrpos)}")
    for r, b in bad_chrpos[:50]:
        print(f"  row {r}: {b!r}")

if __name__ == "__main__":
    phase = sys.argv[1] if len(sys.argv) > 1 else "download"
    if phase == "download":
        download()
    elif phase == "traits":
        traits()
    elif phase == "match":
        match()
    else:
        raise SystemExit("unknown phase: " + phase)
