"""
Track A tests — lock down the frozen feature contract and the parsing logic.

These are the fragile, high-value bits: a silent change to feature order or a hole
in the validation gate would corrupt Tracks B and C without crashing. Everything
here is fast and deterministic — it uses the bundled fixtures, never AMRFinderPlus.

Run:  pytest module1_reader/tests
"""
import pytest

import feature_annotator as fa
from feature_annotator import (
    load_spec, validate_feature_row, parse_amrfinder_tsv, marker_columns,
    target_columns, get_annotator, PrecomputedAnnotator, ContractError,
    load_project_config, spec_project_discrepancies,
)
from build_features import run_genome_reader

FIXTURE_TSV = fa.MODULE_DIR / "fixtures" / "sample_amrfinder.tsv"


@pytest.fixture
def spec():
    return load_spec()


# --------------------------------------------------------------------------- #
# 1. Contract shape — the tripwire for sacred rule #1 (feature order)
# --------------------------------------------------------------------------- #
def test_run_genome_reader_returns_contract_in_order(spec):
    row = run_genome_reader(genome_id="G1", backend="amrfinderplus",
                            spec=spec, tsv_override=FIXTURE_TSV)
    model_cols = [c for c in row if c not in spec["qc_columns"]]
    assert model_cols == spec["feature_order"], "feature columns/order drifted from the contract"
    for qc in spec["qc_columns"]:
        assert qc in row


def test_flags_are_binary_ints(spec):
    row = run_genome_reader(genome_id="G1", backend="amrfinderplus",
                            spec=spec, tsv_override=FIXTURE_TSV)
    for col in spec["feature_order"]:
        assert row[col] in (0, 1)
        assert isinstance(row[col], int)


# --------------------------------------------------------------------------- #
# 2. AMRFinderPlus TSV parsing + alias mapping
# --------------------------------------------------------------------------- #
def test_parse_detects_markers_and_maps_alias(spec):
    flags = parse_amrfinder_tsv(FIXTURE_TSV, spec)
    # direct hits in the fixture
    assert flags["blaKPC-2"] == 1
    assert flags["gyrA_S83L"] == 1
    assert flags["aac(3)-IIa"] == 1
    # OqxB in the TSV must map to our combined marker oqxAB via marker_aliases
    assert flags["oqxAB"] == 1
    # a marker not in the fixture stays absent
    assert flags["blaNDM-1"] == 0


def test_unknown_symbol_is_dropped_not_added(tmp_path, spec):
    tsv = tmp_path / "unknown.tsv"
    tsv.write_text(
        "Element symbol\tType\n"
        "blaFOO-999\tAMR\n"        # not in our vocabulary
        "blaKPC-2\tAMR\n"
    )
    flags = parse_amrfinder_tsv(tsv, spec)
    assert set(flags.keys()) == set(marker_columns(spec)), "row shape changed on unknown marker"
    assert flags["blaKPC-2"] == 1
    assert "blaFOO-999" not in flags


def test_targets_default_present(spec):
    row = run_genome_reader(genome_id="G1", backend="amrfinderplus",
                            spec=spec, tsv_override=FIXTURE_TSV)
    for col in target_columns(spec):
        assert row[col] == 1  # present-until-proven-absent default (documented TODO)


# --------------------------------------------------------------------------- #
# 3. The validation gate — reject off-contract rows LOUDLY
# --------------------------------------------------------------------------- #
def test_validation_rejects_missing_column(spec):
    row = {c: 0 for c in spec["feature_order"]}
    del row[spec["feature_order"][0]]           # drop a required feature
    row["qc_complete"], row["qc_contigs"] = 1.0, 0
    with pytest.raises(ContractError):
        validate_feature_row(row, spec)


def test_validation_rejects_non_binary_flag(spec):
    row = {c: 0 for c in spec["feature_order"]}
    row[spec["feature_order"][0]] = 7           # not 0/1
    row["qc_complete"], row["qc_contigs"] = 1.0, 0
    with pytest.raises(ContractError):
        validate_feature_row(row, spec)


def test_validation_coerces_string_and_numeric_types(spec):
    row = {c: "0" for c in spec["feature_order"]}
    row[spec["feature_order"][0]] = "1"
    row["qc_complete"], row["qc_contigs"] = "0.95", "12"
    out = validate_feature_row(row, spec)
    assert out[spec["feature_order"][0]] == 1 and isinstance(out[spec["feature_order"][0]], int)
    assert out["qc_complete"] == pytest.approx(0.95) and isinstance(out["qc_complete"], float)
    assert out["qc_contigs"] == 12 and isinstance(out["qc_contigs"], int)


# --------------------------------------------------------------------------- #
# 4. Precomputed backend (bring-your-own-dataset)
# --------------------------------------------------------------------------- #
def test_precomputed_loads_row(tmp_path, spec):
    table = tmp_path / "byo.csv"
    table.write_text("genome_id,blaNDM-1,armA\nBYO1,1,1\n")
    ann = PrecomputedAnnotator(spec, table)
    row = ann.annotate_validated("BYO1", None)
    assert row["blaNDM-1"] == 1 and row["armA"] == 1
    assert row["blaKPC-2"] == 0


def test_precomputed_unknown_genome_id_raises(tmp_path, spec):
    table = tmp_path / "byo.csv"
    table.write_text("genome_id,blaNDM-1\nBYO1,1\n")
    ann = PrecomputedAnnotator(spec, table)
    with pytest.raises(ContractError):
        ann.annotate_validated("NOPE", None)


# --------------------------------------------------------------------------- #
# 5. Backend registry + spec resolution
# --------------------------------------------------------------------------- #
def test_get_annotator_unknown_backend_raises(spec):
    with pytest.raises(ValueError):
        get_annotator("does-not-exist", spec)


def test_spec_falls_back_to_bundled_sample():
    # No shared Phase-0 spec exists in this repo yet, so load_spec() must resolve
    # to Track A's bundled sample.
    assert not fa.SHARED_SPEC_PATH.exists()
    s = load_spec()
    assert s["species"] == "Klebsiella pneumoniae"
    assert set(s["drugs"]) == {"meropenem", "ciprofloxacin", "gentamicin", "ceftazidime"}


# --------------------------------------------------------------------------- #
# 6. Connection to Phase 0's shared project config
# --------------------------------------------------------------------------- #
def test_project_config_is_readable():
    project = load_project_config()
    if project is None:
        pytest.skip("data/config/project.json not present in this checkout")
    assert "antibiotics" in project and len(project["antibiotics"]) >= 1


def test_sample_spec_stays_aligned_with_project_config(spec):
    # Tripwire: fails if our sample drifts from the team's declared species/drugs.
    project = load_project_config()
    if project is None:
        pytest.skip("data/config/project.json not present in this checkout")
    assert spec_project_discrepancies(spec, project) == []
