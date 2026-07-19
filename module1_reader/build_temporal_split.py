"""
Build the prospective (temporal) split: train on isolates collected on/before a
cutoff year, test on isolates after it. Collection years come from BV-BRC.

  python -m module1_reader.build_temporal_split \
      --selected data/manifests/selected_genomes.csv \
      --meta     <genome_metadata.csv with cluster_id> \
      --labels   data/manifests/labels.csv \
      --cutoff 2014 \
      --out data/manifests/split_manifest_temporal.csv

The past period is further split into train + a grouped 15% calibration holdout
(whole cgMLST clusters), so calibration stays temporally in-period.
"""
from __future__ import annotations
import argparse, json, re, time, urllib.request
import numpy as np
import pandas as pd


def fetch_years(ids, batch=200):
    rows = []
    for i in range(0, len(ids), batch):
        chunk = ids[i:i + batch]
        q = ("in(genome_id,(%s))&select(genome_id,collection_year,collection_date)"
             "&limit(%d)&http_accept=application/json" % (",".join(chunk), len(chunk)))
        req = urllib.request.Request("https://www.bv-brc.org/api/genome/?" + q,
                                     headers={"Accept": "application/json"})
        for attempt in range(3):
            try:
                rows.extend(json.load(urllib.request.urlopen(req, timeout=90))); break
            except Exception:
                if attempt == 2: raise
                time.sleep(2)
    return pd.DataFrame(rows)


def year_of(r):
    y = r.get("collection_year")
    if pd.notna(y) and str(y).strip() not in ("", "0"):
        try: return int(float(y))
        except Exception: pass
    m = re.match(r"(19|20)\d{2}", str(r.get("collection_date") or ""))
    return int(m.group(0)) if m else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--selected", default="data/manifests/selected_genomes.csv")
    ap.add_argument("--meta", required=True, help="csv with genome_id + cluster_id")
    ap.add_argument("--labels", default="data/manifests/labels.csv")
    ap.add_argument("--cutoff", type=int, default=2014)
    ap.add_argument("--cal-frac", type=float, default=0.15)
    ap.add_argument("--seed", type=int, default=2015)
    ap.add_argument("--out", default="data/manifests/split_manifest_temporal.csv")
    args = ap.parse_args()

    sel = pd.read_csv(args.selected, dtype={"genome_id": str})
    yr = fetch_years(sel.genome_id.tolist())
    yr["year"] = yr.apply(year_of, axis=1)
    yr = yr.dropna(subset=["year"]); yr["year"] = yr.year.astype(int)

    meta = pd.read_csv(args.meta, dtype={"genome_id": str})[["genome_id", "cluster_id"]]
    t = yr.merge(meta, on="genome_id", how="left")
    t["cluster_id"] = t.cluster_id.fillna("solo_" + t.genome_id)
    t["period"] = np.where(t.year <= args.cutoff, "past", "future")

    rng = np.random.default_rng(args.seed)
    past_clusters = t.loc[t.period == "past", "cluster_id"].unique()
    cal = set(rng.choice(past_clusters, size=int(args.cal_frac * len(past_clusters)), replace=False))
    t["split"] = np.where(t.period == "future", "test",
                          np.where(t.cluster_id.isin(cal), "calibration", "train"))

    sm = t[["genome_id", "cluster_id", "split"]]
    sm.to_csv(args.out, index=False)

    # report drift + leak caveat
    counts = sm.split.value_counts().to_dict()
    tc_clusters = set(t[t.split != "test"].cluster_id)
    leak = int((t[t.split == "test"].cluster_id.isin(tc_clusters)).sum())
    n_test = int((t.split == "test").sum())
    lab = pd.read_csv(args.labels, dtype={"genome_id": str}).merge(t[["genome_id", "period"]], on="genome_id")
    lab["y"] = (lab.phenotype == "Resistant").astype(int)
    drift = lab.pivot_table(index="antibiotic", columns="period", values="y", aggfunc="mean").round(3)
    print("split counts:", counts)
    print(f"temporal-only leak: {leak}/{n_test} test genomes share a cluster with a past isolate "
          f"({100*leak//max(n_test,1)}%)")
    print("R-rate drift past->future:\n", drift.to_string())
    print("wrote", args.out)


if __name__ == "__main__":
    main()
