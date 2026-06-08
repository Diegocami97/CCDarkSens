"""
Mediator labels shared by ``qedark.entry`` and ``qcdark.entry``.

FDM exponent n in (α m_e / q)^n for heavy (n=0) vs light (n=2) mediators.
"""

from __future__ import annotations

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
