import xarray as xr


def make_dataset_from_var_dict(var_dict, time_axis):
    time_values_yr = time_axis.to("yr").m
    return xr.Dataset(
        data_vars={
            var: (["time"], value.m, {"units": str(value.units)})
            for var, value in var_dict.items()
        },
        coords={"time": (["time"], time_values_yr, {"units": "yr"})},
    )
