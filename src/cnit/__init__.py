"""
CNit - Carbon-Nitrogen Interactions in Terrestrial ecosystems.

CNit is a process-based terrestrial biogeochemistry model that simulates
coupled carbon and nitrogen cycles in terrestrial ecosystems.
"""

# Version
try:
    from importlib.metadata import version, PackageNotFoundError
except ImportError:
    from importlib_metadata import version, PackageNotFoundError

try:
    __version__ = version("cnit")
except PackageNotFoundError:
    __version__ = "unknown"

# Suppress expected warnings
import warnings

try:
    from pint import UnitStrippedWarning
    warnings.filterwarnings('ignore', category=UnitStrippedWarning)
except ImportError:
    pass

# Core model classes
from cnit.physics.carbon_nitrogen_cycle import (
    CarbonNitrogenCycleModel,
    CarbonNitrogenCycleModelConfig,
    CarbonNitrogenCycleExperimentConfig,
)

# Convenient aliases
CNitModel = CarbonNitrogenCycleModel
CNitModelConfig = CarbonNitrogenCycleModelConfig
CNitExpConfig = CarbonNitrogenCycleExperimentConfig

# Utility imports
from cnit.utils.units import Q

# Define public API
__all__ = [
    # Version
    "__version__",
    # Main classes (full names)
    "CarbonNitrogenCycleModel",
    "CarbonNitrogenCycleModelConfig",
    "CarbonNitrogenCycleExperimentConfig",
    # Convenient aliases
    "CNitModel",
    "CNitModelConfig",
    "CNitExpConfig",
    # Utilities
    "Q",
]