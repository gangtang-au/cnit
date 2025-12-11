from typing import Any, Callable
import openscm_units
import attr
import pint


UNIT_REGISTRY: pint.UnitRegistry = openscm_units.unit_registry
"""
Unit registry used internally
"""

Q = UNIT_REGISTRY.Quantity


def check_units(
        required_dimension: str,
) -> Callable[[Any, attr.Attribute, pint.Quantity], None]:
    """
    Check units of class attribute

    Intended to be used as a validator with :func:`attrs.field`

    Parameters
    ----------
    required_dimension
        Dimension that the input is required to have

    Returns
    -------
        Function that will validate that the intended dimension is passed
    """

    def check_unit_internal(
            self: Any, attribute: attr.Attribute, value: pint.Quantity
    ) -> None:
        """
        Check units of attribute

        Parameters
        ----------
        self
            Object instance

        attribute
            Attribute to check

        value
            Value to check

        Raises
        ------
        :obj:`pint.errors.DimensionalityError`
            Units are not the correct dimensionality
        """
        if not value.check(required_dimension):
            raise pint.errors.DimensionalityError(
                value,
                "a quantity with dimensionality",
                value.dimensionality,
                UNIT_REGISTRY.get_dimensionality(required_dimension),
                extra_msg=f" to set attribute `{attribute.name}`",
            )

    return check_unit_internal
