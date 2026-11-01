import csv
from pathlib import Path

from attrs import fields
import numpy as np
import pint
import pytest

from cnit.configs.model_config import CarbonNitrogenCycleModelConfig
from cnit.physics.carbon_nitrogen_cycle import CarbonNitrogenCycleModel
from cnit.physics.carbon_cycle import (
    CarbonCycleCalculator,
    CarbonCycleModel,
    CarbonTurnoverCalculator,
    CarbonNPPLPRCalculator,
)
from cnit.physics.effects import (
    EffectCO2Calculator,
    EffectTemperatureCalculator,
    EffectCarbonNitrogenCouplingCalculator,
    EffectLanduseCalculator,
)
from cnit.physics.landuse import LanduseCalculator
from cnit.physics.nitrogen_cycle import (
    NitrogenCycleCalculator,
    NitrogenCycleModel,
    NitrogenPUBNFCalculator,
    NitrogenTurnoverCalculator,
)
from cnit.utils.units import Q


def _calculator_defaults():
    carbon_cycle = CarbonCycleCalculator()
    carbon_turnover = CarbonTurnoverCalculator()
    nitrogen_cycle = NitrogenCycleCalculator()
    nitrogen_turnover = NitrogenTurnoverCalculator()
    return [
        LanduseCalculator(),
        EffectCO2Calculator(),
        EffectTemperatureCalculator(),
        EffectCarbonNitrogenCouplingCalculator(),
        EffectLanduseCalculator(),
        CarbonNPPLPRCalculator(),
        carbon_turnover,
        carbon_cycle,
        CarbonCycleModel(carbon_turnover, carbon_cycle),
        NitrogenPUBNFCalculator(),
        nitrogen_turnover,
        nitrogen_cycle,
        NitrogenCycleModel(nitrogen_turnover, nitrogen_cycle),
    ]


def test_explicit_model_defaults_match_calculator_defaults():
    config = CarbonNitrogenCycleModelConfig()
    defaults = _calculator_defaults()

    for item in fields(type(config)):
        value = getattr(config, item.name)
        if item.name in {"switch_Cpls", "switch_Npls", "switch_N"}:
            continue
        matching_values = [
            getattr(calculator, item.name)
            for calculator in defaults
            if hasattr(calculator, item.name)
        ]
        assert matching_values, f"No calculator default found for {item.name}"
        assert value == matching_values[0], item.name


def test_parameter_csv_uses_explicit_file_and_model_name():
    csv_path = Path(__file__).parents[1] / "data" / "cali_cmip6_best_fitting.csv"
    config = CarbonNitrogenCycleModelConfig.from_parameter_csv(
        csv_path, model_name="ACCESS-ESM1-5"
    )
    with csv_path.open(encoding="utf-8", newline="") as csv_file:
        rows = {row["parameter"]: row for row in csv.DictReader(csv_file)}

    assert config.frc_rgr == Q(float(rows["frc_rgr"]["ACCESS-ESM1-5"]), "1")
    assert config.switch_N == Q(float(rows["switch_N"]["ACCESS-ESM1-5"]), "1")
    assert config.switch_Npls == (1, 1, 1, 1)


def test_parameter_csv_requires_explicit_file_and_model_name():
    with pytest.raises(TypeError):
        CarbonNitrogenCycleModelConfig.from_parameter_csv()
    with pytest.raises(TypeError):
        CarbonNitrogenCycleModelConfig.from_parameter_csv("calibration.csv")


@pytest.mark.parametrize(
    "header, missing_column",
    [
        ("unit,example\n", "parameter"),
        ("parameter,example\n", "unit"),
    ],
)
def test_parameter_csv_requires_parameter_and_unit_columns(
    tmp_path, header, missing_column
):
    csv_path = tmp_path / "missing_column.csv"
    csv_path.write_text(header, encoding="utf-8")

    with pytest.raises(ValueError, match=f"missing required column.*{missing_column}"):
        CarbonNitrogenCycleModelConfig.from_parameter_csv(csv_path, "example")


def test_parameter_csv_rejects_duplicate_parameter_rows(tmp_path):
    csv_path = tmp_path / "duplicate_parameter.csv"
    csv_path.write_text(
        "parameter,unit,example\nNPP0,GtC/yr,50\nNPP0,GtC/yr,60\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="NPP0.*more than once"):
        CarbonNitrogenCycleModelConfig.from_parameter_csv(csv_path, "example")


def test_missing_nitrogen_pool_parameter_freezes_nitrogen_without_changing_switch(
    tmp_path,
):
    csv_path = tmp_path / "carbon_only.csv"
    csv_path.write_text(
        "parameter,unit,example\n" "switch_N,dimensionless,1\n" "NPP0,GtC/yr,50\n",
        encoding="utf-8",
    )

    config = CarbonNitrogenCycleModelConfig.from_parameter_csv(csv_path, "example")

    assert config.switch_N == Q(1, "1")
    assert config.switch_Npls == (0, 0, 0, 0)


def test_config_can_be_created_from_parameter_dict():
    config = CarbonNitrogenCycleModelConfig.from_dict({"NPP0": Q(42, "GtC/yr")})

    assert config.NPP0 == Q(42, "GtC/yr")
    assert config.CplsP0 == Q(400, "GtC")


def test_model_config_validates_parameter_units_in_field_validator():
    with pytest.raises(pint.errors.DimensionalityError):
        CarbonNitrogenCycleModelConfig.from_dict({"tau_CplsP": Q(10, "GtC")})


def test_from_config_does_not_delegate_to_from_dict(monkeypatch):
    def fail_if_called(cls, params):
        raise AssertionError("from_config should construct the model directly")

    monkeypatch.setattr(
        CarbonNitrogenCycleModel, "from_dict", classmethod(fail_if_called)
    )
    config = CarbonNitrogenCycleModelConfig()

    model = CarbonNitrogenCycleModel.from_config(config)

    assert model.calc_lu.frc_rgr == config.frc_rgr
    assert model.calc_eff_CO2.CO2ref == config.CO2ref
    assert model.calc_eff_dT.sns_NPP2dT == config.sns_NPP2dT
    assert model.calc_eff_CN.sns_NPP2PUdef == config.sns_NPP2PUdef
    assert model.calc_eff_CN.NPP0 == config.NPP0
    assert model.calc_npp_lpr.NPP0 == config.NPP0
    assert model.calc_pu_bnf.PU0 == config.PU0
    assert model.mdl_c_cycle.CplsP0 == config.CplsP0
    assert model.mdl_c_cycle.switch_Cpls == config.switch_Cpls
    assert model.mdl_n_cycle.NplsP0 == config.NplsP0
    assert model.mdl_n_cycle.switch_Npls == config.switch_Npls
    assert model.mdl_n_cycle.calc_n_turnover.frc_nLSgas == config.frc_nLSgas
    assert model.switch_N == config.switch_N


@pytest.mark.parametrize(
    "switch_n, pool_n, expected", [(0, 10, 0), (1, 0, 0), (1, 10, 1)]
)
def test_nitrogen_pool_switch_uses_switch_and_initial_pool(
    tmp_path, switch_n, pool_n, expected
):
    csv_path = tmp_path / "nitrogen.csv"
    csv_path.write_text(
        "parameter,unit,example\n"
        f"switch_N,dimensionless,{switch_n}\n"
        f"NplsP0,GtN,{pool_n}\n",
        encoding="utf-8",
    )

    config = CarbonNitrogenCycleModelConfig.from_parameter_csv(csv_path, "example")

    assert config.switch_N == Q(switch_n, "1")
    assert config.switch_Npls == (expected,) * 4


def test_carbon_pool_switch_uses_initial_plant_pool(tmp_path):
    csv_path = tmp_path / "carbon_pool.csv"
    csv_path.write_text(
        "parameter,unit,example\n" "CplsP0,GtC,0\n",
        encoding="utf-8",
    )

    config = CarbonNitrogenCycleModelConfig.from_parameter_csv(csv_path, "example")

    assert config.switch_Cpls == (0, 0, 0)


def test_switches_are_immutable_but_can_be_reassigned():
    config = CarbonNitrogenCycleModelConfig.from_dict(
        {"switch_Cpls": [0, 1, 0], "switch_Npls": [1, 0, 1, 0]}
    )
    model = CarbonCycleModel(
        CarbonTurnoverCalculator(),
        CarbonCycleCalculator(),
        switch_Cpls=[0, 1, 0],
    )
    nitrogen_model = NitrogenCycleModel(
        NitrogenTurnoverCalculator(),
        NitrogenCycleCalculator(),
        switch_Npls=[1, 0, 1, 0],
    )

    assert config.switch_Cpls == (0, 1, 0)
    assert config.switch_Npls == (1, 0, 1, 0)
    assert model.switch_Cpls == (0, 1, 0)
    assert nitrogen_model.switch_Npls == (1, 0, 1, 0)

    with pytest.raises(TypeError):
        config.switch_Cpls[0] = 1
    with pytest.raises(TypeError):
        model.switch_Cpls[0] = 1
    with pytest.raises(TypeError):
        nitrogen_model.switch_Npls[0] = 0

    config.switch_Cpls = [1, 0, 1]
    assert config.switch_Cpls == (1, 0, 1)
    with pytest.raises(ValueError, match="0 or 1"):
        config.switch_Cpls = [1, 2, 1]

    model.switch_Cpls = [1, 0, 1]
    assert model.switch_Cpls == (1, 0, 1)

    with pytest.raises(ValueError, match="0 or 1"):
        model.switch_Cpls = [1, 2, 1]


@pytest.mark.parametrize(
    "name, value",
    [
        ("frc_rgr", Q(1.01, "1")),
        ("frc_nLSgas", Q(-0.01, "1")),
        ("frc_NPP2P", Q(float("nan"), "1")),
    ],
)
def test_config_rejects_invalid_fraction_values(name, value):
    with pytest.raises(ValueError, match=name):
        CarbonNitrogenCycleModelConfig.from_dict({name: value})


@pytest.mark.parametrize(
    "name, values",
    [
        ("switch_Cpls", [1, 2, 1]),
        ("switch_Cpls", [1, 0.5, 1]),
        ("switch_Npls", [1, 1, -1, 1]),
        ("switch_Npls", [1, 1, 1, 2]),
    ],
)
def test_config_rejects_non_binary_pool_switches(name, values):
    with pytest.raises(ValueError, match="0 or 1"):
        CarbonNitrogenCycleModelConfig.from_dict({name: values})


@pytest.mark.parametrize(
    "model_factory, switch_name, values",
    [
        (
            lambda **kwargs: CarbonCycleModel(
                CarbonTurnoverCalculator(), CarbonCycleCalculator(), **kwargs
            ),
            "switch_Cpls",
            [1, 2, 1],
        ),
        (
            lambda **kwargs: NitrogenCycleModel(
                NitrogenTurnoverCalculator(), NitrogenCycleCalculator(), **kwargs
            ),
            "switch_Npls",
            [1, 1, 2, 1],
        ),
    ],
)
def test_cycle_models_reject_non_binary_pool_switches(
    model_factory, switch_name, values
):
    with pytest.raises(ValueError, match="0 or 1"):
        model_factory(**{switch_name: values})


@pytest.mark.parametrize(
    "parameters, message",
    [
        ({"frc_NPP2P": Q(0.7, "1"), "frc_NPP2L": Q(0.4, "1")}, "frc_NPP2P"),
        ({"frc_PU2P": Q(0.8, "1"), "frc_PU2L": Q(0.3, "1")}, "frc_PU2P"),
    ],
)
def test_config_rejects_fraction_partitions_above_one(parameters, message):
    with pytest.raises(ValueError, match=message):
        CarbonNitrogenCycleModelConfig.from_dict(parameters)


@pytest.mark.parametrize(
    "name, value",
    [
        ("tau_rgr", Q(0, "yr")),
        ("tau_CplsP", Q(-1, "yr")),
        ("tau_NplsM", Q(float("inf"), "yr")),
    ],
)
def test_config_rejects_nonpositive_or_nonfinite_turnover_times(name, value):
    with pytest.raises(ValueError, match=name):
        CarbonNitrogenCycleModelConfig.from_dict({name: value})


@pytest.mark.parametrize("name", ["CplsP0", "NplsM0"])
def test_config_rejects_negative_initial_pool_sizes(name):
    unit = "GtC" if name.startswith("C") else "GtN"
    with pytest.raises(ValueError, match=name):
        CarbonNitrogenCycleModelConfig.from_dict({name: Q(-1, unit)})


@pytest.mark.parametrize(
    "name", [item.name for item in fields(CarbonNitrogenCycleModelConfig)]
)
@pytest.mark.parametrize("invalid", [np.nan, np.inf, -np.inf])
def test_all_numeric_parameters_reject_nonfinite_values(name, invalid):
    default = getattr(CarbonNitrogenCycleModelConfig(), name)
    if isinstance(default, pint.Quantity):
        value = Q(invalid, default.units)
    else:
        value = list(default)
        value[0] = invalid
    with pytest.raises(ValueError, match=name):
        CarbonNitrogenCycleModelConfig.from_dict({name: value})


def test_finite_signed_sensitivities_and_zero_fluxes_remain_valid():
    config = CarbonNitrogenCycleModelConfig.from_dict(
        {
            "sns_NPP2dT": Q(-0.1, "1/K"),
            "NPP0": Q(0, "GtC/yr"),
        }
    )
    assert config.sns_NPP2dT == Q(-0.1, "1/K")
    assert config.NPP0 == Q(0, "GtC/yr")
