"""Experiment forcing validation at the public configuration boundary."""

from attrs import evolve
import numpy as np
import pytest

from cnit import CNitExpConfig, Q


FORCING_UNITS = {
    "dT_s": "K",
    "CO2_s": "ppm",
    "CemsLUnet_s": "GtC/yr",
    "CemsLUgrs_s": "GtC/yr",
    "NflxAD_s": "GtN/yr",
    "NflxFT_s": "GtN/yr",
    "NemsLUnet_s": "GtN/yr",
    "NemsLUgrs_s": "GtN/yr",
    "NemsLUmin_s": "GtN/yr",
}


def experiment():
    return CNitExpConfig(
        name="baseline",
        switch_N=Q(0, "1"),
        time_axis=Q([1850, 1851, 1852], "yr"),
        **{name: Q(np.zeros(3), unit) for name, unit in FORCING_UNITS.items()},
    )


@pytest.mark.parametrize("name,unit", FORCING_UNITS.items())
@pytest.mark.parametrize("value", [np.nan, np.inf, -np.inf])
def test_rejects_nonfinite_forcing_including_unused_inputs(name, unit, value):
    with pytest.raises(ValueError, match=f"{name} must contain only finite values"):
        evolve(experiment(), **{name: Q([0, value, 0], unit)})


def test_accepts_finite_forcing_and_unused_nonzero_placeholders():
    config = experiment()
    config = evolve(
        config,
        dT_s=Q([-1, 0, 1], "K"),
        CO2_s=Q([280, 281, 282], "ppm"),
        CemsLUgrs_s=Q([1, 2, 3], "GtC/yr"),
        NflxAD_s=Q([1, 2, 3], "GtN/yr"),
    )
    assert config.CemsLUgrs_s.m.tolist() == [1, 2, 3]
    assert config.NflxAD_s.m.tolist() == [1, 2, 3]
