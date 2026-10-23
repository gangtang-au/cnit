"""Compare default and calibrated model results with the implementation on ``main``."""

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
ACCESS_MODEL = "ACCESS-ESM1-5"
EXPERIMENTS = ("1pctCO2", "hist_ssp585")
ACCESS_DIAGNOSTICS = (
    "CplsP",
    "CplsL",
    "CplsS",
    "NplsP",
    "NplsL",
    "NplsS",
    "NplsM",
    "CflxNPP",
    "NflxBNF",
    "NflxPU",
)


def _serialize_config(config):
    """Convert configured quantities to JSON values for the runner process."""
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


def _run_revision(source, config_path, output_path):
    subprocess.run(
        [
            sys.executable,
            str(RUNNER),
            str(source),
            str(FORCING_FILE),
            str(config_path),
            str(output_path),
        ],
        cwd=REPO_ROOT,
        check=True,
    )


def _max_errors(actual, expected):
    """Return maximum absolute and element-wise relative errors.

    Relative error is calculated where the expected value is nonzero. If every
    expected value is zero, it is undefined and returned as NaN.
    """
    absolute_error = np.abs(actual - expected)
    nonzero_expected = expected != 0
    if np.any(nonzero_expected):
        relative_error = absolute_error[nonzero_expected] / np.abs(
            expected[nonzero_expected]
        )
        max_relative = float(np.max(relative_error))
    else:
        max_relative = float("nan")
    return float(np.max(absolute_error)), max_relative


def test_default_model_matches_main_and_report_access_errors(tmp_path):
    """Check default results exactly and report calibrated ACCESS differences."""
    try:
        baseline_source = _extract_baseline_source(tmp_path / "baseline")
    except subprocess.CalledProcessError as exc:
        pytest.fail(
            "Could not extract the baseline source. Ensure the baseline ref "
            "exists locally, or set CNIT_BASELINE_REF to an available commit.\n"
            f"{exc}"
        )

    default_config = CarbonNitrogenCycleModelConfig()
    access_config = CarbonNitrogenCycleModelConfig.from_parameter_csv(
        PARAMETER_FILE,
        model_name=ACCESS_MODEL,
    )
    config_path = tmp_path / "configs.json"
    config_path.write_text(
        json.dumps(
            {
                "default": _serialize_config(default_config),
                "ACCESS calibrated": _serialize_config(access_config),
            }
        ),
        encoding="utf-8",
    )

    current_output = tmp_path / "current.npz"
    baseline_output = tmp_path / "baseline.npz"
    _run_revision(REPO_ROOT / "src", config_path, current_output)
    _run_revision(baseline_source, config_path, baseline_output)

    with np.load(current_output) as current, np.load(baseline_output) as baseline:
        current_keys = set(current.files)
        baseline_keys = set(baseline.files)
        assert current_keys == baseline_keys, (
            f"Output keys differ: current-only={sorted(current_keys - baseline_keys)}, "
            f"main-only={sorted(baseline_keys - current_keys)}"
        )

        # The default configuration has an integer-year regrowth timescale, so
        # this comparison checks that the refactor preserves default model results.
        for experiment in EXPERIMENTS:
            prefix = f"default::{experiment}::"
            for key in sorted(key for key in current_keys if key.startswith(prefix)):
                np.testing.assert_allclose(
                    current[key],
                    baseline[key],
                    rtol=1e-8,
                    atol=1e-8,
                    err_msg=f"Default model differs from main for {key}",
                )

        print("\nACCESS calibrated emulation: maximum differences from main")
        print("experiment       variable    max_abs_error    max_rel_error")
        for experiment in EXPERIMENTS:
            for variable in ACCESS_DIAGNOSTICS:
                key = f"ACCESS calibrated::{experiment}::{variable}"
                max_absolute, max_relative = _max_errors(current[key], baseline[key])
                print(
                    f"{experiment:<16} {variable:<10} "
                    f"{max_absolute:>14.6g} {max_relative:>15.6g}"
                )
