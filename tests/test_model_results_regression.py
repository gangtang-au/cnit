"""Compare calibrated model results with the implementation on ``main``."""

import csv
import json
import os
from dataclasses import fields
from pathlib import Path
import subprocess
import sys
import tarfile

import numpy as np
import pytest

from cnit.configs.model_config import CarbonNitrogenCycleModelConfig


TEST_DIR = Path(__file__).parent
REPO_ROOT = TEST_DIR.parent
PARAMETER_FILE = REPO_ROOT / "data" / "cali_cmip6_best_fitting.csv"
FORCING_FILE = REPO_ROOT / "data" / "CMIP6_cali_forcing.nc"
RUNNER = TEST_DIR / "_cmip6_version_runner.py"
EXPERIMENTS = ("1pctCO2", "hist_ssp585")


def _serialize_config(config):
    """Convert configured quantities to JSON values for the baseline process."""
    values = {}
    for config_field in fields(config):
        value = getattr(config, config_field.name)
        if hasattr(value, "units") and hasattr(value, "magnitude"):
            values[config_field.name] = {
                "magnitude": float(value.magnitude),
                "unit": str(value.units),
            }
        else:
            values[config_field.name] = value
    return values


def _extract_baseline_source(destination):
    baseline_ref = os.environ.get("CNIT_BASELINE_REF", "main")
    destination.mkdir(parents=True, exist_ok=True)
    archive_path = destination / "baseline.tar"
    with archive_path.open("wb") as archive_file:
        subprocess.run(
            ["git", "archive", baseline_ref, "src/cnit"],
            cwd=REPO_ROOT,
            stdout=archive_file,
            check=True,
        )

    with tarfile.open(archive_path) as archive:
        archive.extractall(destination)
    return destination / "src"


def _run_revision(source, revision, config_dir, output_path):
    subprocess.run(
        [
            sys.executable,
            str(RUNNER),
            str(source),
            str(PARAMETER_FILE),
            str(FORCING_FILE),
            str(config_dir),
            str(output_path),
            revision,
        ],
        cwd=REPO_ROOT,
        check=True,
    )


def test_cmip6_results_match_main_for_calibrated_models(tmp_path):
    """Current and main code produce identical outputs under one runtime."""
    try:
        baseline_source = _extract_baseline_source(tmp_path / "baseline")
    except subprocess.CalledProcessError as exc:
        pytest.fail(
            "Could not extract the baseline source. Ensure the baseline ref "
            "exists locally, or set CNIT_BASELINE_REF to an available commit.\n"
            f"{exc}"
        )

    with PARAMETER_FILE.open(newline="", encoding="utf-8") as parameter_file:
        model_names = next(csv.reader(parameter_file))[2:]

    config_dir = tmp_path / "configs"
    config_dir.mkdir()
    for model_name in model_names:
        config = CarbonNitrogenCycleModelConfig.from_parameter_csv(
            PARAMETER_FILE,
            model_name=model_name,
        )
        (config_dir / f"{model_name}.json").write_text(
            json.dumps(_serialize_config(config)),
            encoding="utf-8",
        )

    current_output = tmp_path / "current.npz"
    baseline_output = tmp_path / "baseline.npz"
    _run_revision(REPO_ROOT / "src", "current", config_dir, current_output)
    _run_revision(baseline_source, "baseline", config_dir, baseline_output)

    with np.load(current_output) as current, np.load(baseline_output) as baseline:
        current_keys = set(current.files)
        baseline_keys = set(baseline.files)
        assert current_keys == baseline_keys, (
            f"Output keys differ: current-only={sorted(current_keys - baseline_keys)}, "
            f"main-only={sorted(baseline_keys - current_keys)}"
        )

        assert any("::1pctCO2::" in key for key in current_keys)
        assert any("::hist_ssp585::" in key for key in current_keys)
        for key in sorted(current_keys):
            np.testing.assert_allclose(
                current[key],
                baseline[key],
                rtol=1e-8,
                atol=1e-8,
                err_msg=f"Current and main differ for {key}",
            )
