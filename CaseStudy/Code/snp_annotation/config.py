# -*- coding: utf-8 -*-
"""Shared paths, constants, and parsing helpers for the 50-SNP annotation task."""
import os
import re
import csv

# PATH: paths below resolve relative to this file (snp_annotation/); place inputs next to config.py or edit as needed.
BASE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(BASE, "raw")
GWAS_DIR = os.path.join(RAW, "gwas_catalog")
FETCH_LOG = os.path.join(RAW, "fetch_log.txt")

RSID_TXT = os.path.join(BASE, "ADNI_top50_rsID.txt")
VEP_TXT = os.path.join(BASE, "results_from_Ensembl_Log_rank.txt")
LOCUS_CSV = os.path.join(BASE, "VEP_Ensembl116_summary_by_SNP.csv")

EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
VAR_API_BASE = "https://api.ncbi.nlm.nih.gov/variation/v0"
GWAS_ZIP_URL = "https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-associations_ontology-annotated-full.zip"

P_THRESH = 5e-8
WINDOW_BP = 5000
LOCUS_WINDOW_BP = 500000

OUT_DBSNP_CSV = os.path.join(BASE, "dbsnp_annotations.csv")
OUT_COMPARE_CSV = os.path.join(BASE, "compare_vep_dbsnp.csv")
OUT_GWAS_MATCHES = os.path.join(BASE, "gwas_ad_matches.csv")
OUT_SUMMARY = os.path.join(BASE, "summary.txt")

def read_rsids():
    with open(RSID_TXT, encoding="utf-8") as f:
        return [ln.strip() for ln in f if ln.strip()]

def numeric_rsid(rsid):
    return rsid[2:] if rsid.lower().startswith("rs") else rsid

def load_locus():
    """Return dict: rsID -> {locus:int, chr:str, pos:int,
                             gene_overlapped:str, gene_within5kb:str}"""
    out = {}
    with open(LOCUS_CSV, encoding="utf-8", newline="") as f:
        r = csv.reader(f)
        header = next(r)
        for row in r:
            if not row or not row[1]:
                continue
            rsid = row[1].strip()
            loc = row[2].strip()
            chr_, pos = _split_loc(loc)
            out[rsid] = {
                "locus": int(float(row[6].strip())) if row[6].strip() else None,
                "chr": chr_,
                "pos": pos,
                "gene_overlapped": row[3],
                "gene_within5kb": row[4],
            }
    return out

def _split_loc(loc):
    """'chr:start-end' -> (chr, int(start))"""
    chr_, rest = loc.split(":", 1)
    start = rest.split("-", 1)[0]
    return chr_, int(start)

def parse_gene_overlapped(raw):
    """Gene_overlapped column -> set of gene symbols.
    Units separated by ';'; each unit 'SYMBOL[type]: ...' -> SYMBOL before '['.
    Drop unnamed '-' (never split on a lone hyphen)."""
    genes = set()
    for seg in raw.split(";"):
        seg = seg.strip()
        if not seg:
            continue
        sym = seg.split("[", 1)[0].strip()
        if sym and sym != "-":
            genes.add(sym)
    return genes

def parse_gene_within5kb(raw):
    """'Gene within 5 kb' column -> set of gene symbols.
    Units 'SYMBOL(upstream,2274bp)' separated by ';' -> SYMBOL before '('.
    Drop unnamed '-'."""
    genes = set()
    for seg in raw.split(";"):
        seg = seg.strip()
        if not seg:
            continue
        sym = seg.split("(", 1)[0].strip()
        if sym and sym != "-":
            genes.add(sym)
    return genes

def our_gene_set(locus_row):
    """Full 'our-side' gene set for one SNP (from the Locus file columns)."""
    genes = set()
    genes |= parse_gene_overlapped(locus_row["gene_overlapped"])
    genes |= parse_gene_within5kb(locus_row["gene_within5kb"])
    return genes

def parse_vep():
    """Return dict: rsID -> {'chr':str, 'pos':int, 'symbols':set}.
    Uses VEP Location (chr:start-end) and SYMBOL columns."""
    vep = {}
    with open(VEP_TXT, encoding="utf-8", newline="") as f:
        f.readline()
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            cols = line.split("\t")
            rsid = cols[0].strip()
            loc = cols[1].strip()
            symbol = cols[5].strip()
            if rsid not in vep:
                chr_, pos = _split_loc(loc)
                vep[rsid] = {"chr": chr_, "pos": pos, "symbols": set()}
            if symbol and symbol != "-":
                vep[rsid]["symbols"].add(symbol)
    return vep

def split_genes_gwas(text):
    """Parse MAPPED_GENE(s)/REPORTED GENE(S).
    Split ONLY on ', ', ';', and ' - ' (space-hyphen-space). Never split a lone '-'
    (gene names may contain hyphens, e.g. ADNP-AS1).
    Returns list of (gene, source) where source is 'direct' or 'x-y'."""
    if text is None:
        return []
    t = str(text).strip()
    if t in ("", "NR", "-", "N/A", "NA", "n/a"):
        return []
    source = "x-y" if re.search(r"\s+-\s+", t) else "direct"
    genes = []
    for seg in re.split(r",\s+|\s*;\s*|\s+-\s+", t):
        g = seg.strip()
        if g and g not in ("", "NR", "-", "N/A", "NA", "n/a"):
            genes.append((g, source))
    return genes

def genes_from_gwas_field(text):
    """Convenience: return set of gene symbols parsed from a GWAS gene field."""
    return {g for g, _ in split_genes_gwas(text)}
