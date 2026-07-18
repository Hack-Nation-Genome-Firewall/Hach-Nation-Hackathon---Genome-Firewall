"""
GENOME FIREWALL — Module 1 (Track A): homology-grouped split manifest.

Produces split_manifest.csv (genome_id, cluster_id, split) — a Track A deliverable
that Track B refuses to train without. Two steps (per docs/TRACK_A_HANDOFF.md):

  1. cluster genomes by DNA similarity (Mash/sourmash) -> cluster_id   [PLUGGABLE]
  2. assign WHOLE clusters to train / calibration / test              [DETERMINISTIC]

Track B asserts no homology cluster crosses splits (module2_predictor/contracts.py),
so step 2 obeys that BY CONSTRUCTION — clusters are the unit of assignment, never
individual genomes. That is the whole point: near-identical clones can't straddle
train and test, so evaluation can't be inflated by memorised lineages.

Step 1 needs the genome FASTAs + Mash and is left behind a seam: use precomputed
clusters (load_clusters_csv) now; MashClusterer documents how the real step slots in.
"""
from __future__ import annotations

import csv
import random
from collections import defaultdict
from pathlib import Path
from typing import Iterable, Optional

SPLITS = ("train", "calibration", "test")


def load_clusters_csv(path) -> dict:
    """Read a `genome_id,cluster_id` table into {genome_id: cluster_id}."""
    clusters: dict = {}
    with open(path, newline="") as f:
        for rec in csv.DictReader(f):
            clusters[rec["genome_id"]] = rec["cluster_id"]
    return clusters


def assign_clusters_to_splits(cluster_of: dict, *, fractions=(0.70, 0.15, 0.15),
                              seed: int = 0) -> tuple[dict, dict]:
    """
    Assign whole clusters to train/calibration/test to approximate `fractions`.

    Deterministic given `seed`. Returns (genome_id -> split, cluster_id -> split).
    Greedy largest-cluster-first, each cluster going to the split with the biggest
    remaining genome-count deficit — which keeps proportions close while never
    splitting a cluster. Raises if there are fewer clusters than splits.
    """
    if abs(sum(fractions) - 1.0) > 1e-6:
        raise ValueError(f"fractions must sum to 1.0, got {fractions}")

    members: dict[str, list] = defaultdict(list)
    for genome_id, cluster_id in cluster_of.items():
        members[cluster_id].append(genome_id)
    if len(members) < len(SPLITS):
        raise ValueError(
            f"Need at least {len(SPLITS)} homology clusters to fill {list(SPLITS)}; "
            f"got {len(members)}"
        )

    clusters = list(members.items())
    random.Random(seed).shuffle(clusters)                       # reproducible order
    clusters.sort(key=lambda kv: len(kv[1]), reverse=True)      # largest first, for balance

    total = sum(len(g) for _, g in clusters)
    target = {s: frac * total for s, frac in zip(SPLITS, fractions)}
    current = {s: 0 for s in SPLITS}
    cluster_split: dict = {}
    for cluster_id, genomes in clusters:
        pick = max(SPLITS, key=lambda s: (target[s] - current[s], s))
        cluster_split[cluster_id] = pick
        current[pick] += len(genomes)

    if set(cluster_split.values()) != set(SPLITS):
        raise ValueError(f"Split assignment left an empty split (counts={current}); "
                         f"need more/again-sized clusters")

    genome_split = {g: cluster_split[c] for c, gs in clusters for g in gs}
    return genome_split, cluster_split


def build_split_manifest(cluster_of: dict, out_path, *, fractions=(0.70, 0.15, 0.15),
                         seed: int = 0, feature_genome_ids: Optional[Iterable[str]] = None) -> Path:
    """
    Write split_manifest.csv (genome_id, cluster_id, split).

    If `feature_genome_ids` is given, asserts the cluster genome set matches the
    feature genome set exactly (Track B requires one split row per feature genome).
    """
    if feature_genome_ids is not None:
        fset, cset = set(feature_genome_ids), set(cluster_of)
        if fset != cset:
            raise ValueError(
                "cluster genome set != feature genome set; "
                f"missing_from_clusters={sorted(fset - cset)[:5]} "
                f"missing_from_features={sorted(cset - fset)[:5]}"
            )

    genome_split, _ = assign_clusters_to_splits(cluster_of, fractions=fractions, seed=seed)
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["genome_id", "cluster_id", "split"])
        for genome_id in sorted(cluster_of):
            writer.writerow([genome_id, cluster_of[genome_id], genome_split[genome_id]])
    return out_path


class MashClusterer:
    """
    Pluggable step 1: cluster genome FASTAs by Mash distance. STUB — needs Mash
    installed and the genome assemblies, which aren't in the repo. Documented so the
    deterministic assignment above is usable now against precomputed clusters, and
    the real clustering has an obvious home.

    Real implementation plan:
      * `mash sketch` each FASTA -> a small MinHash sketch;
      * `mash dist` pairwise -> a genome-to-genome distance;
      * link any two genomes closer than `threshold`;
      * take connected components as clusters (each ~= one clonal lineage).
    The threshold is an expert-approved knob (per TEAM_OWNERSHIP.md).
    """

    def __init__(self, threshold: float = 0.001, mash_bin: str = "mash"):
        self.threshold = threshold
        self.mash_bin = mash_bin

    def cluster(self, fasta_paths) -> dict:
        raise NotImplementedError(
            "MashClusterer needs Mash + genome FASTAs (not in the repo). Provide "
            "precomputed clusters via load_clusters_csv() for now. See the class "
            "docstring for the sketch -> dist -> connected-components plan."
        )
