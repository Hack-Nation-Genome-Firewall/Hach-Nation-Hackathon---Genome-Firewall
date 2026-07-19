"""
Track A tests — lock down the feature contract, the parsing logic, and (crucially)
that our output satisfies Track B's contract in module2_predictor/contracts.py.

These are the fragile, high-value bits: a silent change to the feature schema or a
hole in the validation gate would corrupt Tracks B and C without crashing. Everything
here is fast and deterministic — it uses the bundled fixtures, never AMRFinderPlus.

Run:  pytest module1_reader/tests
"""
import pytest

import feature_annotator as fa
from feature_annotator import (
    load_spec, validate_feature_row, parse_amrfinder_tsv, parse_amrfinder_markers,
    marker_columns, target_columns, quality_columns, get_annotator, load_qc_map,
    PrecomputedAnnotator, ContractError, load_project_config, spec_project_discrepancies,
)
from build_features import run_genome_reader, build_features_table

FIXTURE_TSV = fa.MODULE_DIR / "fixtures" / "sample_amrfinder.tsv"


@pytest.fixture
def spec():
    return load_spec()


# --------------------------------------------------------------------------- #
# 1. Contract shape — the tripwire for a silent schema/order change
# --------------------------------------------------------------------------- #
def test_run_genome_reader_returns_contract_in_order(spec):
    row = run_genome_reader(genome_id="G1", backend="amrfinderplus",
                            spec=spec, tsv_override=FIXTURE_TSV)
    flag_cols = [c for c in row if c not in quality_columns(spec)]
    assert flag_cols == marker_columns(spec) + target_columns(spec), \
        "feature columns/order drifted from the contract"
    for qc in quality_columns(spec):
        assert qc in row


def test_flags_are_binary_ints(spec):
    row = run_genome_reader(genome_id="G1", backend="amrfinderplus",
                            spec=spec, tsv_override=FIXTURE_TSV)
    for col in marker_columns(spec) + target_columns(spec):
        assert row[col] in (0, 1)
        assert isinstance(row[col], int)


# --------------------------------------------------------------------------- #
# 2. AMRFinderPlus TSV parsing + alias mapping
# --------------------------------------------------------------------------- #
def test_parse_detects_markers_and_maps_alias(spec):
    flags = parse_amrfinder_tsv(FIXTURE_TSV, spec)
    assert flags["blaKPC-2"] == 1
    assert flags["gyrA_S83L"] == 1
    assert flags["aac(3)-IIa"] == 1
    # OqxB in the TSV must map to our combined marker oqxAB via marker_aliases
    assert flags["oqxAB"] == 1
    # a marker not in the fixture stays absent
    assert flags["blaNDM-1"] == 0


def test_lowercase_oqx_maps_to_oqxab(tmp_path, spec):
    # AMRFinderPlus 4.2.7 emits lowercase oqxA/oqxB; they must map to oqxAB and
    # NOT land in unknowns (regression for the alias-case bug found in the pilot).
    tsv = tmp_path / "oqx.tsv"
    tsv.write_text("Element symbol\tType\noqxA\tAMR\noqxB\tAMR\n")
    flags, unknown = parse_amrfinder_markers(tsv, spec)
    assert flags["oqxAB"] == 1
    assert "oqxA" not in unknown and "oqxB" not in unknown


def test_unknown_symbol_is_preserved_not_dropped(tmp_path, spec):
    tsv = tmp_path / "unknown.tsv"
    tsv.write_text(
        "Element symbol\tType\n"
        "blaFOO-999\tAMR\n"        # not in our vocabulary
        "blaKPC-2\tAMR\n"
    )
    flags, unknown = parse_amrfinder_markers(tsv, spec)
    # the fixed feature vector must not grow a column for the unknown...
    assert set(flags.keys()) == set(marker_columns(spec)), "row shape changed on unknown marker"
    assert flags["blaKPC-2"] == 1
    assert "blaFOO-999" not in flags
    # ...but the unknown must be PRESERVED for review, not silently dropped
    assert unknown == ["blaFOO-999"]


def test_run_genome_reader_surfaces_unknown_markers(tmp_path, spec):
    tsv = tmp_path / "u.tsv"
    tsv.write_text("Element symbol\tType\nblaFOO-999\tAMR\nblaKPC-2\tAMR\n")
    unknown: list = []
    run_genome_reader(genome_id="G", backend="amrfinderplus", spec=spec,
                      tsv_override=tsv, unknown_markers_out=unknown)
    assert unknown == ["blaFOO-999"]


def test_batch_writes_unknown_markers_sidecar(tmp_path, spec):
    tsv = tmp_path / "u.tsv"
    tsv.write_text("Element symbol\tType\nblaFOO-999\tAMR\nblaKPC-2\tAMR\n")
    out = tmp_path / "features.csv"
    build_features_table(
        [{"genome_id": "G1", "source": str(tsv)}],
        backend="amrfinderplus", out_path=out, spec=spec, precomputed_tsv=True,
    )
    sidecar = out.parent / "unknown_markers.csv"
    assert sidecar.exists()
    text = sidecar.read_text()
    assert "genome_id,unknown_marker" in text
    assert "G1,blaFOO-999" in text


def test_targets_default_present(spec):
    row = run_genome_reader(genome_id="G1", backend="amrfinderplus",
                            spec=spec, tsv_override=FIXTURE_TSV)
    for col in target_columns(spec):
        assert row[col] == 1  # present-until-proven-absent default (documented TODO)


# --------------------------------------------------------------------------- #
# 3. The validation gate — reject off-contract rows LOUDLY
# --------------------------------------------------------------------------- #
def _good_row(spec):
    row = {c: 0 for c in marker_columns(spec) + target_columns(spec)}
    qf = spec["quality_features"]
    row[qf["completeness"]], row[qf["contamination"]], row[qf["contigs"]] = 100.0, 0.0, 0
    return row


def test_validation_rejects_missing_column(spec):
    row = _good_row(spec)
    del row[marker_columns(spec)[0]]            # drop a required feature
    with pytest.raises(ContractError):
        validate_feature_row(row, spec)


def test_validation_rejects_non_binary_flag(spec):
    row = _good_row(spec)
    row[marker_columns(spec)[0]] = 7            # not 0/1
    with pytest.raises(ContractError):
        validate_feature_row(row, spec)


def test_validation_coerces_string_and_numeric_types(spec):
    row = {c: "0" for c in marker_columns(spec) + target_columns(spec)}
    first = marker_columns(spec)[0]
    row[first] = "1"
    qf = spec["quality_features"]
    row[qf["completeness"]], row[qf["contamination"]], row[qf["contigs"]] = "95.0", "1.2", "12"
    out = validate_feature_row(row, spec)
    assert out[first] == 1 and isinstance(out[first], int)
    assert out[qf["completeness"]] == pytest.approx(95.0) and isinstance(out[qf["completeness"]], float)
    assert out[qf["contigs"]] == 12 and isinstance(out[qf["contigs"]], int)


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
    assert not fa.SHARED_SPEC_PATH.exists()
    s = load_spec()
    assert s["species"]["name"] == "Klebsiella pneumoniae"
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
    project = load_project_config()
    if project is None:
        pytest.skip("data/config/project.json not present in this checkout")
    assert spec_project_discrepancies(spec, project) == []


# --------------------------------------------------------------------------- #
# 7. Compatibility with Track B's contract (module2_predictor/contracts.py)
#    The real proof of alignment: our spec + output must satisfy Track B's own
#    validators. Skips cleanly if Track B / pandas aren't importable.
# --------------------------------------------------------------------------- #
def test_sample_spec_passes_track_b_feature_spec_validation():
    contracts = pytest.importorskip("module2_predictor.contracts")
    contracts.validate_feature_spec(load_spec())  # raises if our schema is off-contract


def test_feature_row_passes_track_b_inference_validation(spec):
    contracts = pytest.importorskip("module2_predictor.contracts")
    row = run_genome_reader(genome_id="G1", backend="amrfinderplus",
                            spec=spec, tsv_override=FIXTURE_TSV)
    contracts.validate_inference_row(row, spec)   # raises if Track B would reject the row


# --------------------------------------------------------------------------- #
# 8. Homology-grouped split manifest
# --------------------------------------------------------------------------- #
def _clusters(n_clusters=5, per=4):
    return {f"g{c}_{i}": f"cluster_{c}" for c in range(n_clusters) for i in range(per)}


def test_split_assigns_whole_clusters_and_all_splits_used():
    from split_manifest import assign_clusters_to_splits
    cluster_of = _clusters()
    genome_split, cluster_split = assign_clusters_to_splits(cluster_of, seed=20260718)
    # every genome in a cluster shares that cluster's split (no cluster crosses splits)
    for gid, cid in cluster_of.items():
        assert genome_split[gid] == cluster_split[cid]
    assert set(cluster_split.values()) == {"train", "calibration", "test"}


def test_split_is_deterministic():
    from split_manifest import assign_clusters_to_splits
    cluster_of = _clusters(8, 3)
    a, _ = assign_clusters_to_splits(cluster_of, seed=20260718)
    b, _ = assign_clusters_to_splits(cluster_of, seed=20260718)
    assert a == b


def test_split_too_few_clusters_raises():
    from split_manifest import assign_clusters_to_splits
    with pytest.raises(ValueError):
        assign_clusters_to_splits({"g0": "c0", "g1": "c1"})   # only 2 clusters


def test_split_manifest_satisfies_track_b_contract(tmp_path, spec):
    pd = pytest.importorskip("pandas")
    contracts = pytest.importorskip("module2_predictor.contracts")
    from split_manifest import build_split_manifest

    cluster_of = _clusters(6, 3)                       # 18 genomes, 6 clusters
    manifest = build_split_manifest(cluster_of, tmp_path / "split.csv",
                                    seed=20260718, feature_genome_ids=list(cluster_of))

    features = pd.DataFrame(
        [{"genome_id": g, **{m: 0 for m in spec["model_features"]}} for g in cluster_of]
    )
    labels = pd.DataFrame(
        [{"genome_id": g, "antibiotic": d, "phenotype": "Susceptible",
          "evidence": spec["expected_label_evidence"]}
         for g in cluster_of for d in spec["drugs"]]
    )
    splits = pd.read_csv(manifest, dtype={"genome_id": str, "cluster_id": str})

    # The real proof: Track B's own validator accepts our manifest + matching frames,
    # including its "no homology cluster crosses splits" assertion.
    contracts.validate_training_frames(features, labels, splits, spec)


# --------------------------------------------------------------------------- #
# 9. QC wiring — real values from the cohort manifest; unknown -> no-call
# --------------------------------------------------------------------------- #
def test_qc_unknown_defaults_to_none(spec):
    # No qc_source -> QC columns are present but None (unknown), never faked clean.
    row = run_genome_reader(genome_id="G1", backend="amrfinderplus",
                            spec=spec, tsv_override=FIXTURE_TSV)
    for col in quality_columns(spec):
        assert col in row and row[col] is None


def test_qc_filled_from_source(spec):
    qf = spec["quality_features"]
    qc_source = {"G1": {"completeness": 98.5, "contamination": 1.2, "contigs": 87}}
    row = run_genome_reader(genome_id="G1", backend="amrfinderplus", spec=spec,
                            tsv_override=FIXTURE_TSV, qc_source=qc_source)
    assert row[qf["completeness"]] == 98.5
    assert row[qf["contamination"]] == 1.2
    assert row[qf["contigs"]] == 87


def test_load_qc_map_reads_checkm_columns(tmp_path):
    csvf = tmp_path / "selected_genomes.csv"
    csvf.write_text(
        "genome_id,checkm_completeness,checkm_contamination,contigs\n"
        "573.1,99.4,0.3,84\n"
        "573.2,88.0,6.1,700\n"
    )
    qc = load_qc_map(csvf)
    assert qc["573.1"] == {"completeness": 99.4, "contamination": 0.3, "contigs": 84}
    assert qc["573.2"]["contigs"] == 700


def test_qc_aligns_with_track_b_quality_gate(spec):
    predict = pytest.importorskip("module2_predictor.predict")
    qf = spec["quality_features"]
    base = {qf["completeness"]: 99.0, qf["contamination"]: 0.5, qf["contigs"]: 90}
    assert predict._quality_gate(base, spec)["status"] == "pass"
    bad = {qf["completeness"]: 80.0, qf["contamination"]: 0.5, qf["contigs"]: 90}
    assert predict._quality_gate(bad, spec)["status"] == "fail"
    unknown = {qf["completeness"]: None, qf["contamination"]: None, qf["contigs"]: None}
    assert predict._quality_gate(unknown, spec)["status"] == "unknown"
