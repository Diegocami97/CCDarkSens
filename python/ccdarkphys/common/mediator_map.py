# ============================================================================
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  File: mediator_map.py
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  mediator_map.py -- Mediator labels shared by ``qedark.entry`` and
#  ``qcdark.entry``.
# ============================================================================

"""
Mediator labels shared by ``qedark.entry`` and ``qcdark.entry``.

FDM exponent n in (α m_e / q)^n for heavy (n=0) vs light (n=2) mediators.
"""

from __future__ import annotations

# ----------------------------------------------------------------------------
# MEDIATOR_TO_FDM_INDEX
#   Mediator label (string or number) -> the FDM exponent n: heavy/massive -> 0, light/massless -> 2.
# ----------------------------------------------------------------------------
MEDIATOR_TO_FDM_INDEX: dict[str | int, int] = {
    "heavy": 0,
    "massive": 0,
    "0": 0,
    0: 0,
    "light": 2,
    "massless": 2,
    "2": 2,
    2: 2,
}
