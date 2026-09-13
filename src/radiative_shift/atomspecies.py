import numpy as np
from dataclasses import dataclass, field
from .constants import C


@dataclass
class AtomSpecies:
    F0: float = 0.0
    F: float = 0.0
    J0: float = 0.0
    J: float = 0.0
    I: float = 0.0
    lambda_nm: float = 1.0
    gamma: float = 0.0

    wavenumber: float = field(init=False)
    omega: float = field(init=False)
    m0: list = field(init=False)
    m: list = field(init=False)

    def __post_init__(self):
        self.wavenumber = 2 * np.pi / (self.lambda_nm * 1e-7)
        self.omega = C * self.wavenumber
        self.lbar_nm = self.lambda_nm / 2 / np.pi
        self.lbar = 1 / self.wavenumber  # cm, like the model coordinates
        self.m0 = _spin_sublevels(self.F0)
        self.m = _spin_sublevels(self.F)


def _spin_sublevels(F: float) -> list[float]:
    n2 = 2 * F
    if F < 0 or abs(n2 - round(n2)) > 1e-12:
        raise ValueError("F must be a non-negative integer or half-integer")
    n2 = int(round(n2))
    return [m / 2 for m in range(-n2, n2 + 1, 2)]
