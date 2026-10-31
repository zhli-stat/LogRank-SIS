# -*- coding: utf-8 -*-
"""Step A: fetch dbSNP (NCBI E-utilities) annotations for the 50 rsIDs.
- esummary db=snp retmode=json, one request per rsID (numeric id; 'rs' prefix fails).
- raw JSON saved to raw/dbsnp_<rsid>.json; URL+time logged to raw/fetch_log.txt.
- merge info: for merged rsIDs (merged_sort=1 or snp_id != rs number), resolve the
  merge target via the NCBI Variation API (refsnp/<numeric>), saved to raw/.
- outputs dbsnp_annotations.csv.
"""
import json
import os
import csv
import time
import urllib.request
import datetime

import config as C

def log(msg):
    line = msg
    with open(C.FETCH_LOG, "a", encoding="utf-8") as f:
        f.write(line + "\n")
    print(line, flush=True)

def fetch_text(url, retries=2):
    for attempt in range(retries + 1):
        try:
            req = urllib.request.Request(
                url, headers={"User-Agent": "50snp-annotation/1.0 (contact: local)"}
            )
            with urllib.request.urlopen(req, timeout=60) as r:
                body = r.read().decode("utf-8", "replace")
                return r.status, body
        except Exception as e:  # noqa: BLE001
            last = str(e)
            if attempt < retries:
                time.sleep(1.5)
    return None, last

def esummary_snp(rsid):
    numeric = C.numeric_rsid(rsid)
    url = f"{C.EUTILS_BASE}/esummary.fcgi?db=snp&id={numeric}&retmode=json"
    t0 = datetime.datetime.now().isoformat(timespec="seconds")
    status, text = fetch_text(url)
    note = "numeric-id (rs-prefix invalid for db=snp esummary)"
    log(f"{t0}\t{url}\tstatus={status}\t{rsid}\t{note}")
    if status != 200 or text is None:
        return None
    rawpath = os.path.join(C.RAW, f"dbsnp_{rsid}.json")
    with open(rawpath, "w", encoding="utf-8") as f:
        f.write(text)
    return json.loads(text)

def parse_esummary(rsid, data):
    numeric = C.numeric_rsid(rsid)
    result = data.get("result", {})
    rec = result.get(numeric)
    if rec is None or not isinstance(rec, dict):
        return {"status": "unresolved", "error": data.get("error") or "no result record"}
    genes = [g.get("name", "") for g in rec.get("genes", []) if g.get("name")]
    chrpos = rec.get("chrpos", "") or ""
    chr_, pos = "", ""
    if chrpos and ":" in chrpos:
        chr_, pos = chrpos.split(":", 1)
    return {
        "status": "ok",
        "uid": str(rec.get("uid", "")),
        "snp_id": str(rec.get("snp_id", "")),
        "chr": chr_,
        "pos": pos,
        "genes": genes,
        "fxn_class": rec.get("fxn_class", "") or "",
        "merged_sort": str(rec.get("merged_sort", "")),
        "text": rec.get("text", "") or "",
        "acc": rec.get("acc", "") or "",
    }

def resolve_merge(rsid, parsed):
    """Query Variation API for merged rsIDs. Returns (merged_into, merge_date)."""
    numeric = C.numeric_rsid(rsid)
    url = f"{C.VAR_API_BASE}/refsnp/{numeric}"
    t0 = datetime.datetime.now().isoformat(timespec="seconds")
    status, text = fetch_text(url)
    log(f"{t0}\t{url}\tstatus={status}\t{rsid}\tmerge-resolution")
    if status != 200 or text is None:
        return "", ""
    rawpath = os.path.join(C.RAW, f"dbsnp_{rsid}_varapi.json")
    with open(rawpath, "w", encoding="utf-8") as f:
        f.write(text)
    try:
        vd = json.loads(text)
    except Exception:
        return "", ""
    msd = vd.get("merged_snapshot_data") or {}
    merged_into = ",".join(str(x) for x in msd.get("merged_into", []))
    merge_date = msd.get("proxy_time", "")
    return merged_into, merge_date

def main():
    os.makedirs(C.RAW, exist_ok=True)
    rsids = C.read_rsids()
    rows = []
    for i, rsid in enumerate(rsids):
        numeric = C.numeric_rsid(rsid)
        data = esummary_snp(rsid)
        if data is None:
            rows.append({"rsID": rsid, "dbsnp_chr": "", "dbsnp_pos": "",
                         "dbsnp_genes": "", "dbsnp_fxn_class": "", "snp_id": "",
                         "merged_sort": "", "merged_into": "", "merge_date": "",
                         "acc": "", "status": "unresolved", "note": "request failed"})
            time.sleep(0.4)
            continue
        p = parse_esummary(rsid, data)
        merged_into = ""
        merge_date = ""
        note = ""
        if p["status"] == "unresolved":
            note = p["error"] or "unresolved"
        else:
            is_merged = (p["merged_sort"] == "1") or (p["snp_id"] != numeric)
            if is_merged:
                merged_into, merge_date = resolve_merge(rsid, p)
                note = f"merged -> rs{merged_into} (merged_sort={p['merged_sort']}, text={p['text']})" \
                    if merged_into else \
                    f"merged (merged_sort={p['merged_sort']}, text={p['text']})"
        rows.append({
            "rsID": rsid,
            "dbsnp_chr": p.get("chr", ""),
            "dbsnp_pos": p.get("pos", ""),
            "dbsnp_genes": ";".join(p.get("genes", [])),
            "dbsnp_fxn_class": p.get("fxn_class", ""),
            "snp_id": p.get("snp_id", ""),
            "merged_sort": p.get("merged_sort", ""),
            "merged_into": merged_into,
            "merge_date": merge_date,
            "acc": p.get("acc", ""),
            "status": p["status"],
            "note": note,
        })
        time.sleep(0.4)

    fields = ["rsID", "dbsnp_chr", "dbsnp_pos", "dbsnp_genes", "dbsnp_fxn_class",
              "snp_id", "merged_sort", "merged_into", "merge_date", "acc", "status", "note"]
    with open(C.OUT_DBSNP_CSV, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow(r)
    n_ok = sum(1 for r in rows if r["status"] == "ok")
    n_merged = sum(1 for r in rows if r["merged_into"])
    print(f"\nDone: {len(rows)} rsIDs, {n_ok} ok, {n_merged} merged, "
          f"{len(rows)-n_ok} unresolved -> {C.OUT_DBSNP_CSV}")

if __name__ == "__main__":
    main()
