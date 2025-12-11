import xarray as xr


def make_dataset_from_var_dict(var_dict, time_axis):
    return xr.Dataset(
        data_vars={
            var: (
                ["time"],
                value.m,
                {"units": str(value.units)}
            )
            for var, value in var_dict.items()
        },
        coords={"time": time_axis},
    )