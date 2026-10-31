# -*- coding: utf-8 -*-
"""Supplement tasks (no re-download). Run from the snp_annotation dir.

1. Enrich gwas_ad_matches.csv with GWAS-side CHR_ID, CHR_POS, SNPS,
   STUDY ACCESSION, INITIAL SAMPLE SIZE, and our-SNP-to-locus distance (bp).
2. nearest_ad_hit.csv: for each of 33 loci, the nearest tier=A (P<5e-8) record
   on the same chromosome, regardless of gene name.
3. pmid_39998322_rows.csv: all AD rows (MAPPED_TRAIT contains 'alzheimer')
   for PUBMEDID 39998322.
"""
import os
import csv
import re

BASE = os.path.dirname(os.path.abspath(__file__))
GWAS_TSV = os.path.join(BASE, "raw", "gwas_catalog",
                        "gwas-catalog-download-associations-alt-full.tsv")
MATCHES_CSV = os.path.join(BASE, "gwas_ad_matches.csv")
FILTERED_CSV = os.path.join(BASE, "raw", "gwas_catalog", "gwas_ad_filtered.csv")
LOCUS_CSV = os.path.join(BASE, "VEP_Ensembl116_summary_by_SNP.csv")
NEAREST_CSV = os.path.join(BASE, "nearest_ad_hit.csv")
PMID_CSV = os.path.join(BASE, "pmid_39998322_rows.csv")

I_PUBMEDID = 1
I_DISEASE = 7
I_SAMPLE = 8
I_CHR_ID = 11
I_CHR_POS = 12
I_MAPPED_GENE = 14
I_SNPS = 21
I_PVALUE = 27
I_MAPPED_TRAIT = 34
I_STUDY_ACC = 36

def g(cols, i):
    return cols[i].strip() if i < len(cols) else ""

def split_chrpos(text):
    """Split CHR_POS on ';' and ' x '; return (int list, unparseable list)."""
    nums, bad = [], []
    if not text:
        return nums, bad
    for seg in re.split(r"\s*;\s*|\s+x\s+", text):
        s = seg.strip()
        if not s:
            continue
        try:
            nums.append(int(s))
        except ValueError:
            bad.append(s)
    return nums, bad

def norm_chr(c):
    c = str(c).strip().lower()
    if c.startswith("chr"):
        c = c[3:]
    return {"23": "x", "24": "y", "m": "mt"}.get(c, c)

def load_locus():
    """rsID -> {'locus': int, 'chr': str, 'pos': int}"""
    out = {}
    with open(LOCUS_CSV, encoding="utf-8-sig", newline="") as f:
        r = csv.reader(f)
        next(r)
        for row in r:
            if not row or len(row) < 7 or not row[1]:
                continue
            rsid = row[1].strip()
            loc = row[2].strip()
            chr_, rest = loc.split(":", 1)
            start = int(rest.split("-", 1)[0])
            out[rsid] = {"locus": int(float(row[6].strip())),
                         "chr": chr_, "pos": start}
    return out

def load_matches():
    fieldnames = None
    rows = []
    with open(MATCHES_CSV, encoding="utf-8-sig", newline="") as f:
        rd = csv.DictReader(f)
        fieldnames = rd.fieldnames
        for r in rd:
            rows.append(r)
    return fieldnames, rows

def load_tier_a_source_rows():
    rows = set()
    with open(FILTERED_CSV, encoding="utf-8-sig", newline="") as f:
        for r in csv.DictReader(f):
            if r.get("tier") == "A":
                rows.add(int(r["source_row"]))
    return rows

def main():
    locus = load_locus()
    match_fields, matches = load_matches()
    match_source_rows = {int(m["source_row"]) for m in matches}
    tier_a_rows = load_tier_a_source_rows()

    enrich = {}
    nearest_rec = {}
    pmid_rows = []
    bad_pos = []

    with open(GWAS_TSV, encoding="utf-8", newline="") as f:
        f.readline()
        for lineno, line in enumerate(f, start=2):
            cols = line.rstrip("\n").split("\t")
            pmid = g(cols, I_PUBMEDID)
            trait = g(cols, I_MAPPED_TRAIT)

            if pmid == "39998322" and "alzheimer" in trait.lower():
                pmid_rows.append([
                    g(cols, I_DISEASE), trait, g(cols, I_STUDY_ACC),
                    g(cols, I_SAMPLE), g(cols, I_SNPS), g(cols, I_PVALUE),
                ])

            if lineno in match_source_rows:
                enrich[lineno] = {
                    "CHR_ID": g(cols, I_CHR_ID),
                    "CHR_POS": g(cols, I_CHR_POS),
                    "SNPS": g(cols, I_SNPS),
                    "STUDY ACCESSION": g(cols, I_STUDY_ACC),
                    "INITIAL SAMPLE SIZE": g(cols, I_SAMPLE),
                }

            if lineno in tier_a_rows:
                pos_ints, pos_bad = split_chrpos(g(cols, I_CHR_POS))
                if pos_bad:
                    bad_pos.append((lineno, pos_bad))
                nearest_rec[lineno] = {
                    "CHR_ID": norm_chr(g(cols, I_CHR_ID)),
                    "positions": pos_ints,
                    "SNPS": g(cols, I_SNPS),
                    "MAPPED_GENE": g(cols, I_MAPPED_GENE),
                    "P-VALUE": g(cols, I_PVALUE),
                    "PUBMEDID": g(cols, I_PUBMEDID),
                    "CHR_POS_raw": g(cols, I_CHR_POS),
                }

    new_fields = ["CHR_ID", "CHR_POS", "SNPS", "STUDY ACCESSION",
                  "INITIAL SAMPLE SIZE", "distance"]
    fieldnames = list(match_fields) + new_fields
    with open(MATCHES_CSV, "w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for m in matches:
            sr = int(m["source_row"])
            e = enrich.get(sr, {})
            row = dict(m)
            for k in new_fields[:5]:
                row[k] = e.get(k, "")
            our_chr = norm_chr(m["our_chr"])
            our_pos = int(m["our_pos"])
            gw_chr = norm_chr(e.get("CHR_ID", ""))
            pos_ints, _ = split_chrpos(e.get("CHR_POS", ""))
            if gw_chr == our_chr and pos_ints:
                row["distance"] = min(abs(our_pos - p) for p in pos_ints)
            else:
                row["distance"] = ""
            w.writerow(row)

    loci = {}
    for rsid, info in locus.items():
        loci.setdefault(info["locus"], []).append((rsid, info["chr"], info["pos"]))

    by_chr = {}
    for sr, rec in nearest_rec.items():
        by_chr.setdefault(rec["CHR_ID"], []).append((sr, rec))

    nearest_out = []
    for lnum in sorted(loci):
        snps = loci[lnum]
        lchr = norm_chr(snps[0][1])
        best = None
        for sr, rec in by_chr.get(lchr, []):
            pv = rec["P-VALUE"]
            try:
                pvf = float(pv)
            except ValueError:
                pvf = float("inf")
            for p in rec["positions"]:
                for (rsid, _ch, pos) in snps:
                    d = abs(pos - p)
                    key = (d, pvf, sr)
                    if best is None or key < (best[0], best[1], best[2]):
                        best = (d, pvf, sr, rsid, pos, rec)
        if best is None:
            nearest_out.append({"locus": lnum, "our_rsID": "", "our_pos": "",
                                "nearest_gwas_snps": "none on chr", "CHR_POS": "",
                                "distance": "", "P-VALUE": "", "PUBMEDID": "",
                                "MAPPED_GENE": "", "source_row": ""})
        else:
            d, _pvf, sr, our_rsid, our_pos, rec = best
            nearest_out.append({"locus": lnum, "our_rsID": our_rsid,
                                "our_pos": our_pos,
                                "nearest_gwas_snps": rec["SNPS"],
                                "CHR_POS": rec["CHR_POS_raw"],
                                "distance": d, "P-VALUE": rec["P-VALUE"],
                                "PUBMEDID": rec["PUBMEDID"],
                                "MAPPED_GENE": rec["MAPPED_GENE"],
                                "source_row": sr})
    nf = ["locus", "our_rsID", "our_pos", "nearest_gwas_snps", "CHR_POS",
          "distance", "P-VALUE", "PUBMEDID", "MAPPED_GENE", "source_row"]
    with open(NEAREST_CSV, "w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=nf)
        w.writeheader()
        for r in nearest_out:
            w.writerow(r)

    with open(PMID_CSV, "w", encoding="utf-8-sig", newline="") as f:
        w = csv.writer(f)
        w.writerow(["DISEASE/TRAIT", "MAPPED_TRAIT", "STUDY ACCESSION",
                    "INITIAL SAMPLE SIZE", "SNPS", "P-VALUE"])
        for r in pmid_rows:
            w.writerow(r)

    n_pmid_sig = 0
    for r in pmid_rows:
        try:
            if float(r[5]) < 5e-8:
                n_pmid_sig += 1
        except ValueError:
            pass
    n_none = sum(1 for r in nearest_out if r["nearest_gwas_snps"] == "none on chr")
    print(f"enriched matches rows      = {len(matches)}")
    print(f"nearest rows               = {len(nearest_out)}  (none-on-chr = {n_none})")
    print(f"tier-A records scanned     = {len(nearest_rec)}  (unparseable CHR_POS rows = {len(bad_pos)})")
    print(f"pmid 39998322 alzheimer rows = {len(pmid_rows)}  (P<5e-8 = {n_pmid_sig})")
    print("wrote:", MATCHES_CSV)
    print("wrote:", NEAREST_CSV)
    print("wrote:", PMID_CSV)

if __name__ == "__main__":
    main()
