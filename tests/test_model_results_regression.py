"""Compare default and calibrated model results with the implementation on ``main``."""

import csv
import json
import os
from pathlib import Path
import subprocess
import sys
import tarfile

from attrs import fields
import numpy as np
import pytest

from cnit.configs.model_config import CarbonNitrogenCycleModelConfig


TEST_DIR = Path(__file__).parent
REPO_ROOT = TEST_DIR.parent
PARAMETER_FILE = REPO_ROOT / "data" / "cali_cmip6_best_fitting.csv"
FORCING_FILE = REPO_ROOT / "data" / "CMIP6_cali_forcing.nc"
RUNNER = TEST_DIR / "_cmip6_version_runner.py"
EXPERIMENTS = ("1pctCO2", "hist_ssp585")
KEY_OUTPUTS = (
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


def _cmip6_model_names():
    with PARAMETER_FILE.open("r", encoding="utf-8", newline="") as csv_file:
        return csv.DictReader(csv_file).fieldnames[2:]


def _serialize_config(config):
    """Convert configured quantities to JSON values for the runner process."""
    values = {}
    for config_field in fields(type(config)):
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


def _run_revision(source, config_path, output_path, *, truncate_tau_rgr=False):
    command = [
        sys.executable,
        str(RUNNER),
        str(source),
        str(FORCING_FILE),
        str(config_path),
        str(output_path),
    ]
    if truncate_tau_rgr:
        command.append("truncate_tau_rgr")
    subprocess.run(
        command,
        cwd=REPO_ROOT,
        check=True,
    )


def test_default_and_all_cmip6_calibrations_report_differences_from_main(
    tmp_path,
):
    """Report key output differences from ``main`` for all model configurations."""
    try:
        baseline_source = _extract_baseline_source(tmp_path / "baseline")
    except subprocess.CalledProcessError as exc:
        pytest.fail(
            "Could not extract the baseline source. Ensure the baseline ref "
            "exists locally, or set CNIT_BASELINE_REF to an available commit.\n"
            f"{exc}"
        )

    default_config = CarbonNitrogenCycleModelConfig()
    calibrated_configs = {
        model_name: _serialize_config(
            CarbonNitrogenCycleModelConfig.from_parameter_csv(
                PARAMETER_FILE,
                model_name=model_name,
            )
        )
        for model_name in _cmip6_model_names()
    }
    config_path = tmp_path / "configs.json"
    config_path.write_text(
        json.dumps(
            {
                "default": _serialize_config(default_config),
                **calibrated_configs,
            }
        ),
        encoding="utf-8",
    )

    current_legacy_output = tmp_path / "current_legacy_tau.npz"
    current_rounded_output = tmp_path / "current_rounded_tau.npz"
    baseline_output = tmp_path / "baseline.npz"
    _run_revision(
        REPO_ROOT / "src",
        config_path,
        current_legacy_output,
        truncate_tau_rgr=True,
    )
    _run_revision(REPO_ROOT / "src", config_path, current_rounded_output)
    _run_revision(baseline_source, config_path, baseline_output)

    with (
        np.load(current_legacy_output) as current_legacy,
        np.load(current_rounded_output) as current_rounded,
        np.load(baseline_output) as baseline,
    ):
        current_keys = set(current_legacy.files)
        rounded_keys = set(current_rounded.files)
        baseline_keys = set(baseline.files)
        assert current_keys == baseline_keys == rounded_keys, (
            f"Output keys differ: current-only={sorted(current_keys - baseline_keys)}, "
            f"rounded-only={sorted(rounded_keys - baseline_keys)}, "
            f"main-only={sorted(baseline_keys - current_keys)}"
        )

        # With integer-year timescales, the current implementation should retain
        # the legacy duration behavior used by main.
        for key in sorted(current_keys):
            np.testing.assert_allclose(
                current_legacy[key],
                baseline[key],
                rtol=1e-8,
                atol=1e-8,
                err_msg=f"Integer tau_rgr differs from main for {key}",
            )

        # Collect maximum differences for key pools and fluxes across all
        # configurations and experiments for the rounded-duration implementation.
        # Relative differences are calculated only where the baseline value is
        # nonzero; zero-baseline discrepancies are represented by the absolute
        # difference.
        differences = {
            name: {
                "max_abs": 0.0,
                "max_rel": 0.0,
                "abs_case": "",
                "rel_case": "",
            }
            for name in KEY_OUTPUTS
        }
        for config_name in ("default", *calibrated_configs):
            for experiment in EXPERIMENTS:
                prefix = f"{config_name}::{experiment}::"
                for name in KEY_OUTPUTS:
                    key = f"{prefix}{name}"
                    if key not in rounded_keys:
                        continue
                    current_values = np.asarray(current_rounded[key], dtype=float)
                    baseline_values = np.asarray(baseline[key], dtype=float)
                    assert current_values.shape == baseline_values.shape, (
                        f"Output shape differs for {key}: "
                        f"current={current_values.shape}, "
                        f"main={baseline_values.shape}"
                    )

                    absolute_difference = np.abs(current_values - baseline_values)
                    max_abs = float(np.nanmax(absolute_difference))
                    metric = differences[name]
                    if max_abs > metric["max_abs"]:
                        metric["max_abs"] = max_abs
                        metric["abs_case"] = f"{config_name}::{experiment}"

                    nonzero_baseline = baseline_values != 0
                    if np.any(nonzero_baseline):
                        relative_difference = (
                            absolute_difference[nonzero_baseline]
                            / np.abs(baseline_values[nonzero_baseline])
                        )
                        max_rel = float(np.nanmax(relative_difference))
                        if max_rel > metric["max_rel"]:
                            metric["max_rel"] = max_rel
                            metric["rel_case"] = f"{config_name}::{experiment}"

        print(
            "\nRounded-duration maximum differences from main across default "
            "and CMIP6 "
            "configurations (relative differences use nonzero main values):"
        )
        print(
            f"{'Output':<14} {'Max absolute':>14} {'Case':<32} "
            f"{'Max relative':>14} {'Case':<32}"
        )
        for name, metric in differences.items():
            print(
                f"{name:<14} {metric['max_abs']:>14.6g} "
                f"{metric['abs_case']:<32} {metric['max_rel']:>13.6%} "
                f"{metric['rel_case']:<32}"
            )
