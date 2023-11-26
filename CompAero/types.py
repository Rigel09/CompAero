"""
This module defines the basic types (Enumerations) for the module
"""
from enum import Enum


class FlowState(Enum):
    """Enum that defines wether a flow is sub / super sonic"""

    SUB_SONIC = "SUB_SONIC"
    SUPER_SONIC = "SUPER_SONIC"


class ShockType(Enum):
    """Enum that defines wether a shock is strong or weak. (For oblique shocks)"""

    WEAK = "WEAK"
    STRONG = "STRONG"
