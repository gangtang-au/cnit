"""Switch validation helpers."""

from typing import Any, Callable

import numpy as np


def check_binary_switches(
    expected_length: int,
) -> Callable[[Any, Any, Any], None]:
    """Create an attrs validator for a fixed-length sequence of 0/1 switches."""

    def validate(instance, attribute, value):
        name = attribute.name
        if not isinstance(value, (list, tuple, np.ndarray)):
            raise TypeError(f"{name} must be a sequence of binary switches.")

        message = f"{name} must contain exactly {expected_length} binary values."

        try:
            values = np.asarray(value)
        except (TypeError, ValueError) as exc:
            raise ValueError(message) from exc

        if values.shape != (expected_length,) or not np.isin(values, (0, 1)).all():
            raise ValueError(message)

    return validate
