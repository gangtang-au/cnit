import csv
import math
from pathlib import Path

from attrs import define, field, fields
import pint

from ..utils.units import Q, check_units


@define
class CarbonNitrogenCycleModelConfig:
    """Model parameters and initial pool settings.

    Field validators check parameter units. ``__attrs_post_init__`` checks fraction
    ranges and sums, positive timescales, and non-negative initial pool sizes.
    """

    frc_rgr: pint.Quantity = field(default=Q(0.9, "1"), validator=check_units("1"))
    tau_rgr: pint.Quantity = field(default=Q(100, "yr"), validator=check_units("yr"))
    method_cLU: pint.Quantity = field(default=Q(0, "1"), validator=check_units("1"))
    method_nLU: pint.Quantity = field(default=Q(1, "1"), validator=check_units("1"))

    CO2ref: pint.Quantity = field(default=Q(284.875, "ppm"), validator=check_units("ppm"))
    CO2b: pint.Quantity = field(default=Q(31, "ppm"), validator=check_units("ppm"))
    sns_CO2_log: pint.Quantity = field(default=Q(1, "1"), validator=check_units("1"))
    sns_CO2_sig: pint.Quantity = field(default=Q(1, "1"), validator=check_units("1"))
    eff_CO2_sig_max: pint.Quantity = field(default=Q(2, "1"), validator=check_units("1"))
    method_CO2_NPP: pint.Quantity = field(default=Q(0, "1"), validator=check_units("1"))

    sns_NPP2dT: pint.Quantity = field(default=Q(0, "1/K"), validator=check_units("1/K"))
    sns_NPP2dT_sig: pint.Quantity = field(default=Q(0, "1/K"), validator=check_units("1/K"))
    method_dT_NPP: pint.Quantity = field(default=Q(0, "1"), validator=check_units("1"))
    sns_LPR2dT: pint.Quantity = field(default=Q(0, "1/K"), validator=check_units("1/K"))
    sns_cLP2dT: pint.Quantity = field(default=Q(0, "1/K"), validator=check_units("1/K"))
    sns_cLD2dT: pint.Quantity = field(default=Q(0, "1/K"), validator=check_units("1/K"))
    sns_cSR2dT: pint.Quantity = field(default=Q(0, "1/K"), validator=check_units("1/K"))
    sns_PU2dT: pint.Quantity = field(default=Q(0, "1/K"), validator=check_units("1/K"))
    sns_BNF2dT: pint.Quantity = field(default=Q(0, "1/K"), validator=check_units("1/K"))
    sns_nLP2dT: pint.Quantity = field(default=Q(0, "1/K"), validator=check_units("1/K"))
    sns_nLD2dT: pint.Quantity = field(default=Q(0, "1/K"), validator=check_units("1/K"))
    sns_nSR2dT: pint.Quantity = field(default=Q(0, "1/K"), validator=check_units("1/K"))
    sns_nLSgas2dT: pint.Quantity = field(default=Q(0, "1/K"), validator=check_units("1/K"))

    sns_NPP2PUdef: pint.Quantity = field(default=Q(0, "yr/GtN"), validator=check_units("yr/GtN"))
    sns_cLP2PUdef: pint.Quantity = field(default=Q(0, "yr/GtN"), validator=check_units("yr/GtN"))
    sns_cLD2PUdef: pint.Quantity = field(default=Q(0, "yr/GtN"), validator=check_units("yr/GtN"))
    sns_cSR2PUdef: pint.Quantity = field(default=Q(0, "yr/GtN"), validator=check_units("yr/GtN"))
    sns_PU2PUdef: pint.Quantity = field(default=Q(0, "yr/GtN"), validator=check_units("yr/GtN"))
    sns_BNF2PUdef: pint.Quantity = field(default=Q(0, "yr/GtN"), validator=check_units("yr/GtN"))
    sns_nLP2PUdef: pint.Quantity = field(default=Q(0, "yr/GtN"), validator=check_units("yr/GtN"))
    sns_nLD2PUdef: pint.Quantity = field(default=Q(0, "yr/GtN"), validator=check_units("yr/GtN"))
    sns_nSR2PUdef: pint.Quantity = field(default=Q(0, "yr/GtN"), validator=check_units("yr/GtN"))
    eff_C_PU_max: pint.Quantity = field(default=Q(2, "1"), validator=check_units("1"))
    sns_PU2NPPrd: pint.Quantity = field(default=Q(2.5, "1"), validator=check_units("1"))
    eff_C_BNF_max: pint.Quantity = field(default=Q(2, "1"), validator=check_units("1"))
    sns_BNF2NPPrd: pint.Quantity = field(default=Q(2.5, "1"), validator=check_units("1"))

    sns_cLP2LUrgr: pint.Quantity = field(default=Q(0, "yr/GtC"), validator=check_units("yr/GtC"))
    sns_cLD2LUrgr: pint.Quantity = field(default=Q(0, "yr/GtC"), validator=check_units("yr/GtC"))
    sns_cSR2LUrgr: pint.Quantity = field(default=Q(0, "yr/GtC"), validator=check_units("yr/GtC"))
    sns_nLP2LUrgr: pint.Quantity = field(default=Q(0, "yr/GtC"), validator=check_units("yr/GtC"))
    sns_nLD2LUrgr: pint.Quantity = field(default=Q(0, "yr/GtC"), validator=check_units("yr/GtC"))
    sns_nSR2LUrgr: pint.Quantity = field(default=Q(0, "yr/GtC"), validator=check_units("yr/GtC"))

    NPP0: pint.Quantity = field(default=Q(50, "GtC/yr"), validator=check_units("GtC/yr"))
    LPR0: pint.Quantity = field(default=Q(0, "GtC/yr"), validator=check_units("GtC/yr"))

    frc_NPP2P: pint.Quantity = field(default=Q(0.55, "1"), validator=check_units("1"))
    frc_NPP2L: pint.Quantity = field(default=Q(0.35, "1"), validator=check_units("1"))
    frc_cLP2L: pint.Quantity = field(default=Q(0.7, "1"), validator=check_units("1"))
    frc_cLD2S: pint.Quantity = field(default=Q(0.3, "1"), validator=check_units("1"))
    frc_cLUgrs2P: pint.Quantity = field(default=Q(0.6, "1"), validator=check_units("1"))
    frc_cLUgrs2L: pint.Quantity = field(default=Q(0.3, "1"), validator=check_units("1"))

    tau_CplsP: pint.Quantity = field(default=Q(10, "yr"), validator=check_units("yr"))
    tau_CplsL: pint.Quantity = field(default=Q(1, "yr"), validator=check_units("yr"))
    tau_CplsS: pint.Quantity = field(default=Q(100, "yr"), validator=check_units("yr"))

    CplsP0: pint.Quantity = field(default=Q(400, "GtC"), validator=check_units("GtC"))
    CplsL0: pint.Quantity = field(default=Q(50, "GtC"), validator=check_units("GtC"))
    CplsS0: pint.Quantity = field(default=Q(1550, "GtC"), validator=check_units("GtC"))

    PU0: pint.Quantity = field(default=Q(0.5, "GtN/yr"), validator=check_units("GtN/yr"))
    BNF0: pint.Quantity = field(default=Q(0.1, "GtN/yr"), validator=check_units("GtN/yr"))
    NCratio_PUreq: pint.Quantity = field(default=Q(0.02, "GtN/GtC"), validator=check_units("GtN/GtC"))
    NetMIN0: pint.Quantity = field(default=Q(0.5, "GtN/yr"), validator=check_units("GtN/yr"))
    sns_NetMIN2dT: pint.Quantity = field(default=Q(0, "1/K"), validator=check_units("1/K"))
    NetMINlt0: pint.Quantity = field(default=Q(0.5, "GtN/yr"), validator=check_units("GtN/yr"))
    sns_NetMINlt2AD: pint.Quantity = field(default=Q(0, "yr/GtN"), validator=check_units("yr/GtN"))

    tau_NplsP: pint.Quantity = field(default=Q(10, "yr"), validator=check_units("yr"))
    tau_NplsL: pint.Quantity = field(default=Q(1, "yr"), validator=check_units("yr"))
    tau_NplsS: pint.Quantity = field(default=Q(100, "yr"), validator=check_units("yr"))
    tau_NplsM: pint.Quantity = field(default=Q(1, "yr"), validator=check_units("yr"))

    frc_BNF2P: pint.Quantity = field(default=Q(0.55, "1"), validator=check_units("1"))
    frc_BNF2L: pint.Quantity = field(default=Q(0.35, "1"), validator=check_units("1"))
    frc_PU2P: pint.Quantity = field(default=Q(0.55, "1"), validator=check_units("1"))
    frc_PU2L: pint.Quantity = field(default=Q(0.35, "1"), validator=check_units("1"))
    frc_nLP2L: pint.Quantity = field(default=Q(0.7, "1"), validator=check_units("1"))
    frc_nLD2S: pint.Quantity = field(default=Q(0.3, "1"), validator=check_units("1"))
    frc_nLSgas: pint.Quantity = field(default=Q(1, "1"), validator=check_units("1"))
    frc_nLUgrs2P: pint.Quantity = field(default=Q(0.6, "1"), validator=check_units("1"))
    frc_nLUgrs2L: pint.Quantity = field(default=Q(0.3, "1"), validator=check_units("1"))

    NplsP0: pint.Quantity = field(default=Q(5, "GtN"), validator=check_units("GtN"))
    NplsL0: pint.Quantity = field(default=Q(1, "GtN"), validator=check_units("GtN"))
    NplsS0: pint.Quantity = field(default=Q(150, "GtN"), validator=check_units("GtN"))
    NplsM0: pint.Quantity = field(default=Q(1, "GtN"), validator=check_units("GtN"))

    time0: pint.Quantity = field(default=Q(1850, "yr"), validator=check_units("yr"))
    switch_N: pint.Quantity = field(default=Q(1, "1"), validator=check_units("1"))
    switch_Cpls: list = field(factory=lambda: [1, 1, 1])
    switch_Npls: list = field(factory=lambda: [1, 1, 1, 1])

    def __attrs_post_init__(self) -> None:
        """Validate physical parameter ranges when the config is created."""
        fraction_names = (
            "frc_rgr",
            "frc_NPP2P",
            "frc_NPP2L",
            "frc_cLP2L",
            "frc_cLD2S",
            "frc_cLUgrs2P",
            "frc_cLUgrs2L",
            "frc_BNF2P",
            "frc_BNF2L",
            "frc_PU2P",
            "frc_PU2L",
            "frc_nLP2L",
            "frc_nLD2S",
            "frc_nLSgas",
            "frc_nLUgrs2P",
            "frc_nLUgrs2L",
        )
        fractions = {}
        for name in fraction_names:
            value = float(getattr(self, name).to("1").m)
            if not math.isfinite(value) or not 0 <= value <= 1:
                raise ValueError(f"{name} must be finite and between 0 and 1.")
            fractions[name] = value

        partition_pairs = (
            ("frc_NPP2P", "frc_NPP2L"),
            ("frc_cLUgrs2P", "frc_cLUgrs2L"),
            ("frc_BNF2P", "frc_BNF2L"),
            ("frc_PU2P", "frc_PU2L"),
            ("frc_nLUgrs2P", "frc_nLUgrs2L"),
        )
        for first, second in partition_pairs:
            if fractions[first] + fractions[second] > 1 + 1e-12:
                raise ValueError(f"{first} + {second} must not exceed 1.")

        time_names = (
            "tau_rgr",
            "tau_CplsP",
            "tau_CplsL",
            "tau_CplsS",
            "tau_NplsP",
            "tau_NplsL",
            "tau_NplsS",
            "tau_NplsM",
        )
        for name in time_names:
            value = float(getattr(self, name).to("yr").m)
            if not math.isfinite(value) or value <= 0:
                raise ValueError(f"{name} must be a positive finite duration in years.")

        pool_names = (
            "CplsP0",
            "CplsL0",
            "CplsS0",
            "NplsP0",
            "NplsL0",
            "NplsS0",
            "NplsM0",
        )
        for name in pool_names:
            unit = "GtC" if name.startswith("C") else "GtN"
            value = float(getattr(self, name).to(unit).m)
            if not math.isfinite(value) or value < 0:
                raise ValueError(f"{name} must be finite and non-negative.")

    @classmethod
    def from_dict(cls, params: dict) -> "CarbonNitrogenCycleModelConfig":
        """Create a config from parameter overrides, using defaults for omissions."""
        valid_fields = {item.name for item in fields(cls)}
        unknown = set(params) - valid_fields
        if unknown:
            raise ValueError(f"Unknown parameter(s) in config: {sorted(unknown)}")
        return cls(**params)

    @classmethod
    def from_parameter_csv(
        cls,
        csv_file,
        model_name: str,
    ) -> "CarbonNitrogenCycleModelConfig":
        """Build a model config from a user-selected calibration table.

        Parameters present in the selected CSV column override the standard
        calculator defaults. Parameters absent from the table retain those defaults.

        Both the CSV path and source model column name must be provided.
        """
        parameter_path = Path(csv_file)
        csv_context = parameter_path.open("r", encoding="utf-8", newline="")
        config_fields = {item.name for item in fields(cls)}

        with csv_context:
            reader = csv.DictReader(csv_context)
            if not reader.fieldnames or model_name not in reader.fieldnames:
                available = reader.fieldnames[2:] if reader.fieldnames else []
                raise ValueError(
                    f"Unknown model {model_name!r}. Available models: "
                    f"{', '.join(available)}"
                )

            values = {}
            for row in reader:
                name = row["parameter"].strip()
                if name not in config_fields:
                    raise ValueError(
                        f"Parameter {name!r} in {parameter_path.name} is not a "
                        "field of CarbonNitrogenCycleModelConfig"
                    )
                raw_value = row[model_name].strip()
                if not raw_value:
                    continue
                unit = row["unit"].strip()
                if unit.lower() == "dimensionless":
                    unit = "1"
                try:
                    values[name] = Q(float(raw_value), unit)
                except (TypeError, ValueError, pint.errors.PintError) as exc:
                    raise ValueError(
                        f"Invalid value for parameter {name!r} and model "
                        f"{model_name!r}: {raw_value!r} {unit}"
                    ) from exc

        # Resolve table values against config defaults before deriving pool switches.
        config = cls.from_dict(values)
        carbon_active = config.CplsP0.m != 0
        config.switch_Cpls = [1, 1, 1] if carbon_active else [0, 0, 0]

        # Freeze nitrogen pools when nitrogen is disabled or NplsP0 is absent/zero.
        nitrogen_active = (
            "NplsP0" in values and config.switch_N.m != 0 and config.NplsP0.m != 0
        )
        config.switch_Npls = [1, 1, 1, 1] if nitrogen_active else [0, 0, 0, 0]

        return config
