# -*- coding: utf-8 -*-
"""Step B + D.
B: VEP vs dbSNP -> compare_vep_dbsnp.csv (position 3-state, gene 4-state).
D: recompute 33 loci from VEP coords (single-linkage 500kb), verify vs Locus column;
   summarize strict(A)/broad match counts at SNP and locus level -> summary.txt.
"""
import csv
import os
import datetime

import config as C

def load_dbsnp():
    rows = {}
    with open(C.OUT_DBSNP_CSV, encoding="utf-8", newline="") as f:
        for r in csv.DictReader(f):
            rows[r["rsID"]] = r
    return rows

def four_state(vep_genes, dbsnp_genes):
    if vep_genes == dbsnp_genes:
        if not vep_genes:
            return "完全一致", "均无基因名"
        return "完全一致", ""
    if vep_genes & dbsnp_genes:
        return "部分一致", ""
    return "不一致", ""

def step_B():
    locus = C.load_locus()
    vep = C.parse_vep()
    dbsnp = load_dbsnp()
    rsids = C.read_rsids()

    fields = ["rsID", "vep_chr", "vep_pos", "dbsnp_chr", "dbsnp_pos",
              "location_match", "vep_genes", "dbsnp_genes", "gene_match",
              "status", "notes"]
    with open(C.OUT_COMPARE_CSV, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for rsid in rsids:
            lr = locus[rsid]
            vchr = vep[rsid]["chr"]
            vpos = vep[rsid]["pos"]
            vgenes = C.our_gene_set(lr)
            dr = dbsnp.get(rsid)
            notes = []
            if dr is None or dr.get("status") != "ok":
                w.writerow({
                    "rsID": rsid, "vep_chr": vchr, "vep_pos": vpos,
                    "dbsnp_chr": "", "dbsnp_pos": "",
                    "location_match": "unresolved", "vep_genes": ";".join(sorted(vgenes)),
                    "dbsnp_genes": "", "gene_match": "unresolved",
                    "status": "unresolved", "notes": "dbSNP no data",
                })
                continue
            dchr = dr.get("dbsnp_chr", "")
            dpos = dr.get("dbsnp_pos", "")
            dgenes = set(dr["dbsnp_genes"].split(";")) if dr.get("dbsnp_genes") else set()
            loc_match = "一致" if (dchr == vchr and dpos.isdigit() and int(dpos) == vpos) else "不一致"
            gene_match, gnote = four_state(vgenes, dgenes)
            if dr.get("merged_into"):
                notes.append(f"merged into rs{dr['merged_into']} ({dr.get('merge_date','')})")
            if gnote:
                notes.append(gnote)
            w.writerow({
                "rsID": rsid, "vep_chr": vchr, "vep_pos": vpos,
                "dbsnp_chr": dchr, "dbsnp_pos": dpos,
                "location_match": loc_match,
                "vep_genes": ";".join(sorted(vgenes)),
                "dbsnp_genes": ";".join(sorted(dgenes)),
                "gene_match": gene_match, "status": "ok",
                "notes": "; ".join(notes),
            })
    print("wrote", C.OUT_COMPARE_CSV)

def recompute_loci():
    locus = C.load_locus()
    vep = C.parse_vep()
    rsids = C.read_rsids()

    parent = {r: r for r in rsids}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    by_chr = {}
    for r in rsids:
        by_chr.setdefault(vep[r]["chr"], []).append((vep[r]["pos"], r))
    for ch, lst in by_chr.items():
        lst.sort()
        for i in range(1, len(lst)):
            if lst[i][0] - lst[i - 1][0] <= C.LOCUS_WINDOW_BP:
                union(lst[i - 1][1], lst[i][1])

    cluster_id = {}
    next_id = 1
    for r in rsids:
        root = find(r)
        if root not in cluster_id:
            cluster_id[root] = next_id
            next_id += 1
    n_clusters = next_id - 1

    diffs = []
    for r in rsids:
        mine = cluster_id[find(r)]
        theirs = locus[r]["locus"]
        if mine != theirs:
            diffs.append((r, mine, theirs))

    my_clusters = {}
    for r in rsids:
        my_clusters.setdefault(cluster_id[find(r)], []).append(r)
    their_clusters = {}
    for r in rsids:
        their_clusters.setdefault(locus[r]["locus"], []).append(r)
    my_sets = {frozenset(v) for v in my_clusters.values()}
    their_sets = {frozenset(v) for v in their_clusters.values()}
    partition_identical = (my_sets == their_sets)

    return n_clusters, diffs, partition_identical

def load_matches():
    with open(C.OUT_GWAS_MATCHES, encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))

def step_D():
    n_clusters, diffs, partition_identical = recompute_loci()
    mr = load_matches()

    by_tier = {"A": [], "broad": []}
    for r in mr:
        by_tier[r["tier"]].append(r)

    def summarize(tier):
        rows = by_tier[tier]
        rsids = sorted({r["our_rsID"] for r in rows})
        loci = sorted({int(r["our_locus"]) for r in rows})
        types = sorted({r["match_type"] for r in rows})
        return rsids, loci, types

    meta = {}
    metapath = os.path.join(C.GWAS_DIR, "download_meta.txt")
    if os.path.exists(metapath):
        with open(metapath, encoding="utf-8") as f:
            for line in f:
                if "=" in line:
                    k, v = line.rstrip("\n").split("=", 1)
                    meta[k] = v
    dbsnp_dates = []
    if os.path.exists(C.FETCH_LOG):
        with open(C.FETCH_LOG, encoding="utf-8") as f:
            for line in f:
                ts = line.split("\t")[0]
                try:
                    dbsnp_dates.append(datetime.datetime.fromisoformat(ts))
                except ValueError:
                    pass

    lines = []
    lines.append("=== GWAS Catalog version ===")
    lines.append(f"  requested_url = {meta.get('requested_url','')}")
    lines.append(f"  final_url     = {meta.get('final_url','')}")
    lines.append(f"  last_modified = {meta.get('last_modified','')}  (source: HTTP Last-Modified header)")
    lines.append(f"  download_time = {meta.get('download_time','')}")
    lines.append("")
    lines.append("=== dbSNP request window ===")
    if dbsnp_dates:
        lines.append(f"  first_request = {min(dbsnp_dates).isoformat(timespec='seconds')}")
        lines.append(f"  last_request  = {max(dbsnp_dates).isoformat(timespec='seconds')}")
    lines.append("")
    lines.append("=== Locus recompute (single-linkage 500kb on VEP GRCh38) ===")
    lines.append(f"  n_loci = {n_clusters}")
    lines.append(f"  partition_identical_to_Locus_column = {partition_identical}")
    if diffs:
        lines.append(f"  numbering_diffs = {len(diffs)}")
        for d in diffs:
            lines.append(f"    {d[0]}: recomputed={d[1]} Locus_col={d[2]}")
    else:
        lines.append("  numbering_diffs = 0 (all 50 SNPs match Locus column exactly)")
    lines.append("")

    for tier, label in (("A", "strict (A-level)"), ("broad", "broad (reference)")):
        rsids, loci, types = summarize(tier)
        lines.append(f"=== matches: {label} ===")
        lines.append(f"  SNP-level: {len(rsids)} / 50")
        lines.append(f"    rsID = {', '.join(rsids) if rsids else '(none)'}")
        lines.append(f"  locus-level: {len(loci)} / 33")
        lines.append(f"    locus = {', '.join(str(x) for x in loci) if loci else '(none)'}")
        lines.append(f"  match_type breakdown = {types}")
        lines.append("")

    with open(C.OUT_SUMMARY, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))

if __name__ == "__main__":
    step_B()
    step_D()
