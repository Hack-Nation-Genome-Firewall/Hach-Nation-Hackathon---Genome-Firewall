from __future__ import annotations

from copy import deepcopy

import pandas as pd
import pytest

from data.bvbrc import clean_amr_records, evaluate_genome_quality, select_cohort
from data.download_genomes import format_fasta


@pytest.fixture
def config() -> dict:
    return {
        "schema_version": 1,
        "source": {"evidence": "Laboratory Method"},
        "species": {"taxon_id": 573},
        "antibiotics": ["meropenem", "ciprofloxacin"],
        "label_policy": {
            "included": ["Susceptible", "Resistant"],
            "excluded": ["Intermediate"],
        },
        "quality_policy": {
            "minimum_checkm_completeness": 90,
            "maximum_checkm_contamination": 5,
            "maximum_contigs": 500,
            "minimum_genome_length": 4_000_000,
            "maximum_genome_length": 7_000_000,
            "missing_quality_policy": "exclude",
        },
        "selection_policy": {
            "maximum_genomes": 3,
            "minimum_drugs_per_genome": 1,
            "stable_seed": 7,
        },
    }


def amr_record(record_id: str, genome_id: str, drug: str, phenotype: str, **updates) -> dict:
    value = {
        "id": record_id,
        "genome_id": genome_id,
        "antibiotic": drug,
        "resistant_phenotype": phenotype,
        "evidence": "Laboratory Method",
        "taxon_id": 573,
        "measurement": "1",
        "pmid": [123],
    }
    value.update(updates)
    return value


def test_clean_amr_records_collapses_duplicates_and_excludes_conflicts(config):
    records = [
        amr_record("1", "g1", "meropenem", "Resistant"),
        amr_record("2", "g1", "meropenem", "Resistant"),
        amr_record("3", "g2", "ciprofloxacin", "Susceptible"),
        amr_record("4", "g2", "ciprofloxacin", "Resistant"),
        amr_record("5", "g3", "meropenem", "Intermediate"),
        amr_record("6", "g4", "meropenem", "Susceptible", evidence="Computational Method"),
        {key: value for key, value in amr_record("7", "g5", "meropenem", "Resistant").items()
         if key != "resistant_phenotype"},
    ]

    clean, conflicts, stats = clean_amr_records(records, config)

    assert clean[["genome_id", "phenotype", "record_count"]].to_dict("records") == [
        {"genome_id": "g1", "phenotype": "Resistant", "record_count": 2}
    ]
    assert conflicts.loc[0, "genome_id"] == "g2"
    assert stats["duplicate_records_collapsed"] == 1
    assert stats["conflicting_pairs_excluded"] == 1
    assert stats["excluded_phenotype"] == 1
    assert stats["wrong_evidence"] == 1
    assert stats["invalid_phenotype"] == 1


def test_quality_policy_fails_missing_and_out_of_range_values(config):
    metadata = pd.DataFrame(
        [
            {
                "genome_id": "good",
                "checkm_completeness": 99,
                "checkm_contamination": 1,
                "contigs": 10,
                "genome_length": 5_500_000,
            },
            {
                "genome_id": "bad",
                "checkm_completeness": 80,
                "checkm_contamination": 8,
                "contigs": 900,
                "genome_length": 3_000_000,
            },
            {"genome_id": "missing"},
        ]
    )

    result = evaluate_genome_quality(metadata, config["quality_policy"]).set_index("genome_id")

    assert bool(result.loc["good", "quality_pass"])
    assert not bool(result.loc["bad", "quality_pass"])
    assert "completeness_below_minimum" in result.loc["bad", "quality_reasons"]
    assert not bool(result.loc["missing", "quality_pass"])
    assert result.loc["missing", "quality_reasons"].startswith("missing:")


def test_selection_is_stable_and_requires_both_classes(config):
    labels = pd.DataFrame(
        [
            {"genome_id": "g1", "antibiotic": "meropenem", "phenotype": "Resistant"},
            {"genome_id": "g1", "antibiotic": "ciprofloxacin", "phenotype": "Susceptible"},
            {"genome_id": "g2", "antibiotic": "meropenem", "phenotype": "Susceptible"},
            {"genome_id": "g2", "antibiotic": "ciprofloxacin", "phenotype": "Resistant"},
            {"genome_id": "g3", "antibiotic": "meropenem", "phenotype": "Resistant"},
            {"genome_id": "g3", "antibiotic": "ciprofloxacin", "phenotype": "Resistant"},
            {"genome_id": "g4", "antibiotic": "meropenem", "phenotype": "Susceptible"},
            {"genome_id": "g4", "antibiotic": "ciprofloxacin", "phenotype": "Susceptible"},
        ]
    )
    metadata = pd.DataFrame(
        [
            {
                "genome_id": genome_id,
                "checkm_completeness": 99,
                "checkm_contamination": 1,
                "contigs": 10,
                "genome_length": 5_500_000,
            }
            for genome_id in ["g1", "g2", "g3", "g4"]
        ]
    )

    genomes_a, labels_a, _ = select_cohort(labels, metadata, config)
    genomes_b, labels_b, _ = select_cohort(labels.sample(frac=1, random_state=3), metadata, deepcopy(config))

    assert genomes_a["genome_id"].tolist() == genomes_b["genome_id"].tolist()
    assert labels_a[["genome_id", "antibiotic"]].equals(labels_b[["genome_id", "antibiotic"]])


def test_fasta_format_is_stable_and_validates_symbols():
    body = format_fasta(
        [
            {"sequence_id": "contig2", "description": "second", "sequence": "NNAC"},
            {"sequence_id": "contig1", "description": "first", "sequence": "acgt"},
        ]
    )
    assert body.decode("ascii") == ">contig1 first\nACGT\n>contig2 second\nNNAC\n"

    with pytest.raises(ValueError, match="invalid FASTA symbols"):
        format_fasta([{"sequence_id": "bad", "sequence": "ACGTZ"}])
