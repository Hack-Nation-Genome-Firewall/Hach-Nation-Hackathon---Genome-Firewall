"""
Track A: assemble the real feature matrix + feature_spec from AMRFinderPlus TSVs,
emitting exactly the Track B contract (module2_predictor/contracts.py):

  features.csv       genome_id + every binary model_feature + target__* + qc_*
  feature_spec.json  schema_version/species/drugs/model_features/marker_evidence/
                     drug_targets/quality_features/quality_policy/expected_label_evidence
  (split_manifest.csv is produced separately by build_split.py)

Marker vocabulary is learned empirically from the cohort's AMRFinderPlus output,
family-collapsed, kept within a prevalence band (drops chromosomal constants and
ultra-rare alleles). Point mutations and porin disruptions are preserved.

Usage:
  python -m module1_reader.assemble_features \
      --tsv-dir <dir of *.tsv> \
      --selected data/manifests/selected_genomes.csv \
      --out-features data/manifests/features.csv \
      --out-spec data/manifests/feature_spec.json \
      --min-prev 0.03 --max-prev 0.97
"""
from __future__ import annotations
import argparse, glob, json, os, re
from pathlib import Path
import pandas as pd

DRUGS = ["meropenem", "ciprofloxacin", "gentamicin", "ceftazidime"]

# AMRFinderPlus 'class' substring -> drug(s) it bears on
CLASS_DRUGS = [
    ("QUINOLONE", ["ciprofloxacin"]),
    ("FLUOROQUINOLONE", ["ciprofloxacin"]),
    ("CARBAPENEM", ["meropenem"]),
    ("CEPHALOSPORIN", ["ceftazidime"]),
    ("BETA-LACTAM", ["meropenem", "ceftazidime"]),  # ESBLs hit ceftazidime; carbapenemases hit meropenem
    ("AMINOGLYCOSIDE", ["gentamicin"]),
]

# essential drug-target genes (present -> gate open); expressed as input feature names
DRUG_TARGETS = {
    "ciprofloxacin": ["target__gyrA", "target__parC"],
    "meropenem": ["target__ftsI", "target__mrdA"],
    "ceftazidime": ["target__ftsI"],
    "gentamicin": ["target__rpsL"],
}
ALL_TARGET_FEATURES = sorted({t for ts in DRUG_TARGETS.values() for t in ts})


def collapse(sym: str):
    """AMRFinderPlus Element symbol -> canonical family/marker, or None to drop."""
    s = str(sym)
    if s.startswith("bla"):
        m = re.match(r"(bla(?:CTX-M|KPC|NDM|OXA|SHV|TEM|IMP|VIM|CMY|DHA))", s)
        if m:
            return "marker__" + m.group(1)
        m = re.match(r"(bla[A-Za-z]+)", s)
        return "marker__" + m.group(1) if m else None
    if s.startswith("aac(6"): return "marker__aac6_Ib"
    if s.startswith("aac(3"): return "marker__aac3"
    if s.startswith("aph"):
        m = re.match(r"aph\((\d+)", s); return f"marker__aph{m.group(1)}" if m else None
    if s.startswith("aad"): return "marker__aadA"
    if s.startswith("ant"):
        m = re.match(r"ant\((\d+)", s); return f"marker__ant{m.group(1)}" if m else None
    m = re.match(r"(rmt[A-Z]|armA|npmA)", s)
    if m: return "marker__" + m.group(1)
    m = re.match(r"(qnr[A-Z])", s)
    if m: return "marker__" + m.group(1)
    if s.startswith("oqx"): return "marker__oqxAB"
    m = re.match(r"(gyrA|parC|parE|gyrB)_([A-Z]\d+)", s)
    if m: return f"marker__{m.group(1)}_{m.group(2)}"
    m = re.match(r"(ompK3[56])", s)
    if m: return f"marker__{m.group(1)}_loss"
    return None


def _cols(df):
    df.columns = [c.strip() for c in df.columns]
    return (("Element symbol" if "Element symbol" in df.columns else "Gene symbol"),
            ("Subtype" if "Subtype" in df.columns else "Element subtype"),
            ("Class" if "Class" in df.columns else "class"),
            ("Method" if "Method" in df.columns else "method"))


def harvest(tsv_dir):
    rows = []
    for tsv in glob.glob(os.path.join(tsv_dir, "*.tsv")):
        gid = os.path.basename(tsv)[:-4]
        try:
            df = pd.read_csv(tsv, sep="\t", dtype=str)
        except Exception:
            continue
        sym_c, sub_c, cls_c, mth_c = _cols(df)
        for _, r in df.iterrows():
            feat = collapse(r.get(sym_c, ""))
            if feat:
                rows.append({"genome_id": gid, "feature": feat,
                             "class": (r.get(cls_c) or ""), "subtype": (r.get(sub_c) or ""),
                             "method": (r.get(mth_c) or "")})
    return pd.DataFrame(rows)


def marker_drugs(harv, marker):
    classes = " ".join(harv.loc[harv.feature == marker, "class"].astype(str)).upper()
    ds = []
    for key, drugs in CLASS_DRUGS:
        if key in classes:
            ds.extend(drugs)
    return sorted(set(ds))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tsv-dir", required=True)
    ap.add_argument("--selected", default="data/manifests/selected_genomes.csv")
    ap.add_argument("--meta", default=None, help="genome metadata csv with checkm_* and contigs")
    ap.add_argument("--out-features", default="data/manifests/features.csv")
    ap.add_argument("--out-spec", default="data/manifests/feature_spec.json")
    ap.add_argument("--min-prev", type=float, default=0.03)
    ap.add_argument("--max-prev", type=float, default=0.97)
    args = ap.parse_args()

    harv = harvest(args.tsv_dir)
    genomes = sorted(harv.genome_id.unique())
    ng = len(genomes)
    prev = harv.groupby("feature").genome_id.nunique() / ng
    markers = sorted(prev[(prev >= args.min_prev) & (prev <= args.max_prev)].index)
    constants = sorted(prev[prev > args.max_prev].index)

    # presence matrix
    pres = (harv[harv.feature.isin(markers)].assign(v=1)
            .pivot_table(index="genome_id", columns="feature", values="v",
                         aggfunc="max", fill_value=0))
    pres = pres.reindex(index=genomes, columns=markers, fill_value=0).astype(int)

    # target gate: essential genes default present; flip on POINT_DISRUPT of that gene
    for tf in ALL_TARGET_FEATURES:
        pres[tf] = 1
    disr = harv[harv.subtype.str.upper().str.contains("POINT_DISRUPT", na=False)]
    for _, r in disr.iterrows():
        base = r.feature.replace("marker__", "").split("_")[0].lower()
        for tf in ALL_TARGET_FEATURES:
            if base and base in tf.lower():
                pres.loc[r.genome_id, tf] = 0

    # QC columns (completeness, contamination, contigs) from metadata/selected file
    meta_path = args.meta or args.selected
    meta = pd.read_csv(meta_path, dtype={"genome_id": str}).set_index("genome_id")
    pres["qc_completeness"] = meta.reindex(pres.index)["checkm_completeness"].astype(float)
    pres["qc_contamination"] = meta.reindex(pres.index)["checkm_contamination"].astype(float)
    pres["qc_contigs"] = meta.reindex(pres.index)["contigs"].astype(float)

    feats = pres.reset_index().rename(columns={"index": "genome_id"})
    Path(args.out_features).parent.mkdir(parents=True, exist_ok=True)
    feats.to_csv(args.out_features, index=False)

    # feature_spec.json in the Track B contract shape
    marker_evidence = {}
    for m in markers:
        ds = marker_drugs(harv, m)
        methods = harv.loc[harv.feature == m, "method"].astype(str)
        mtype = "point_mutation" if "POINT" in " ".join(harv.loc[harv.feature == m, "subtype"].astype(str)).upper() else "acquired_gene"
        marker_evidence[m] = {
            "type": mtype,
            "drugs": ds or DRUGS,  # fall back to all if class unmapped (flag for expert review)
            "source": "AMRFinderPlus --plus, family-collapsed; NCBI Reference Gene Catalog",
            "amrfinder_method": methods.mode().iat[0] if len(methods) else "BLASTX",
        }
    spec = {
        "schema_version": 1,
        "status": "real_pending_domain_review",
        "species": {"name": "Klebsiella pneumoniae", "taxon_id": 573},
        "drugs": DRUGS,
        "model_features": markers,  # target__ and qc_ are inputs, not model_features
        "marker_evidence": marker_evidence,
        "drug_targets": DRUG_TARGETS,
        "quality_features": {"completeness": "qc_completeness",
                             "contamination": "qc_contamination",
                             "contigs": "qc_contigs"},
        "quality_policy": {"minimum_completeness": 90.0, "maximum_contamination": 5.0,
                           "maximum_contigs": 500},
        "expected_label_evidence": "Laboratory Method",
        "annotation": {"tool": "AMRFinderPlus", "tool_version": "4.2.7",
                       "db_version": "2026-05-15.1",
                       "command": "amrfinder -n {fasta} --organism Klebsiella_pneumoniae --plus"},
        "vocabulary_provenance": {"n_genomes": ng, "prevalence_band": [args.min_prev, args.max_prev],
                                  "excluded_constant_genes": constants},
    }
    with open(args.out_spec, "w") as f:
        json.dump(spec, f, indent=2, sort_keys=True)
    print(f"genomes={ng} markers={len(markers)} targets={len(ALL_TARGET_FEATURES)}")
    print("model_features:", markers)
    print("excluded constants:", constants)
    print(f"wrote {args.out_features} {feats.shape} and {args.out_spec}")


if __name__ == "__main__":
    main()
