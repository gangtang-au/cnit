from dataclasses import dataclass, field, fields
import csv
from pathlib import Path
import pint

from ..utils.units import Q


@dataclass
class CarbonNitrogenCycleModelConfig:
    frc_rgr: pint.Quantity = Q(0.9, "1")
    tau_rgr: pint.Quantity = Q(100, "yr")
    method_cLU: pint.Quantity = Q(0, "1")
    method_nLU: pint.Quantity = Q(1, "1")

    CO2ref: pint.Quantity = Q(284.875, "ppm")
    CO2b: pint.Quantity = Q(31, "ppm")
    sns_CO2_log: pint.Quantity = Q(1, "1")
    sns_CO2_sig: pint.Quantity = Q(1, "1")
    eff_CO2_sig_max: pint.Quantity = Q(2, "1")
    method_CO2_NPP: pint.Quantity = Q(0, "1")

    sns_NPP2dT: pint.Quantity = Q(0, "1/K")
    sns_NPP2dT_sig: pint.Quantity = Q(0, "1/K")
    method_dT_NPP: pint.Quantity = Q(0, "1")
    sns_LPR2dT: pint.Quantity = Q(0, "1/K")
    sns_cLP2dT: pint.Quantity = Q(0, "1/K")
    sns_cLD2dT: pint.Quantity = Q(0, "1/K")
    sns_cSR2dT: pint.Quantity = Q(0, "1/K")
    sns_PU2dT: pint.Quantity = Q(0, "1/K")
    sns_BNF2dT: pint.Quantity = Q(0, "1/K")
    sns_nLP2dT: pint.Quantity = Q(0, "1/K")
    sns_nLD2dT: pint.Quantity = Q(0, "1/K")
    sns_nSR2dT: pint.Quantity = Q(0, "1/K")
    sns_nLSgas2dT: pint.Quantity = Q(0, "1/K")

    sns_NPP2PUdef: pint.Quantity = Q(0, "yr/GtN")
    sns_cLP2PUdef: pint.Quantity = Q(0, "yr/GtN")
    sns_cLD2PUdef: pint.Quantity = Q(0, "yr/GtN")
    sns_cSR2PUdef: pint.Quantity = Q(0, "yr/GtN")
    sns_PU2PUdef: pint.Quantity = Q(0, "yr/GtN")
    sns_BNF2PUdef: pint.Quantity = Q(0, "yr/GtN")
    sns_nLP2PUdef: pint.Quantity = Q(0, "yr/GtN")
    sns_nLD2PUdef: pint.Quantity = Q(0, "yr/GtN")
    sns_nSR2PUdef: pint.Quantity = Q(0, "yr/GtN")
    eff_C_PU_max: pint.Quantity = Q(2, "1")
    sns_PU2NPPrd: pint.Quantity = Q(2.5, "1")
    eff_C_BNF_max: pint.Quantity = Q(2, "1")
    sns_BNF2NPPrd: pint.Quantity = Q(2.5, "1")

    sns_cLP2LUrgr: pint.Quantity = Q(0, "yr/GtC")
    sns_cLD2LUrgr: pint.Quantity = Q(0, "yr/GtC")
    sns_cSR2LUrgr: pint.Quantity = Q(0, "yr/GtC")
    sns_nLP2LUrgr: pint.Quantity = Q(0, "yr/GtC")
    sns_nLD2LUrgr: pint.Quantity = Q(0, "yr/GtC")
    sns_nSR2LUrgr: pint.Quantity = Q(0, "yr/GtC")

    NPP0: pint.Quantity = Q(50, "GtC/yr")
    LPR0: pint.Quantity = Q(0, "GtC/yr")

    frc_NPP2P: pint.Quantity = Q(0.55, "1")
    frc_NPP2L: pint.Quantity = Q(0.35, "1")
    frc_cLP2L: pint.Quantity = Q(0.7, "1")
    frc_cLD2S: pint.Quantity = Q(0.3, "1")
    frc_cLUgrs2P: pint.Quantity = Q(0.6, "1")
    frc_cLUgrs2L: pint.Quantity = Q(0.3, "1")

    tau_CplsP: pint.Quantity = Q(10, "yr")
    tau_CplsL: pint.Quantity = Q(1, "yr")
    tau_CplsS: pint.Quantity = Q(100, "yr")

    CplsP0: pint.Quantity = Q(400, "GtC")
    CplsL0: pint.Quantity = Q(50, "GtC")
    CplsS0: pint.Quantity = Q(1550, "GtC")

    PU0: pint.Quantity = Q(0.5, "GtN/yr")
    BNF0: pint.Quantity = Q(0.1, "GtN/yr")
    NCratio_PUreq: pint.Quantity = Q(0.02, "GtN/GtC")
    NetMIN0: pint.Quantity = Q(0.5, "GtN/yr")
    sns_NetMIN2dT: pint.Quantity = Q(0, "1/K")
    NetMINlt0: pint.Quantity = Q(0.5, "GtN/yr")
    sns_NetMINlt2AD: pint.Quantity = Q(0, "yr/GtN")

    tau_NplsP: pint.Quantity = Q(10, "yr")
    tau_NplsL: pint.Quantity = Q(1, "yr")
    tau_NplsS: pint.Quantity = Q(100, "yr")
    tau_NplsM: pint.Quantity = Q(1, "yr")

    frc_BNF2P: pint.Quantity = Q(0.55, "1")
    frc_BNF2L: pint.Quantity = Q(0.35, "1")
    frc_PU2P: pint.Quantity = Q(0.55, "1")
    frc_PU2L: pint.Quantity = Q(0.35, "1")
    frc_nLP2L: pint.Quantity = Q(0.7, "1")
    frc_nLD2S: pint.Quantity = Q(0.3, "1")
    frc_nLSgas: pint.Quantity = Q(1, "1")
    frc_nLUgrs2P: pint.Quantity = Q(0.6, "1")
    frc_nLUgrs2L: pint.Quantity = Q(0.3, "1")

    NplsP0: pint.Quantity = Q(5, "GtN")
    NplsL0: pint.Quantity = Q(1, "GtN")
    NplsS0: pint.Quantity = Q(150, "GtN")
    NplsM0: pint.Quantity = Q(1, "GtN")

    time0: pint.Quantity = Q(1850, "yr")
    switch_N: pint.Quantity = Q(1, "1")
    switch_Cpls: list = field(default_factory=lambda: [1, 1, 1])
    switch_Npls: list = field(default_factory=lambda: [1, 1, 1, 1])

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
