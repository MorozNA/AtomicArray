import numpy as np
from .sigma_matrix import SigmaMatrix
from src.radiative_shift.constants import DDI, HBAR
from src.radiative_shift.model import GeneralModel
from src.radiative_shift.tools import d_up, d_down, reshape_to_matrix


class SigmaEin(SigmaMatrix):
    """
    Self-energy matrix for a (F0=0, F=1) - atomic medium.

    Matrix elements have units of energy, while the corresponding
    resolvent elements have units of inverse energy.

    The medium basis is ordered by atom. For each atom, the magnetic
    sublevels appear in the order m = -1, 0, +1; for example:
    (atom 0, -1), (atom 0, 0), (atom 0, +1),
    (atom 1, -1), (atom 1, 0), (atom 1, +1), ...
    """

    def __init__(self, model: GeneralModel):
        self.medium_atom = model.medium_atom
        if (self.medium_atom.F0, self.medium_atom.F) != (0, 1):
            raise ValueError("The medium must be a V atom with F0=0 and F=1")
        self.wavenumber = model.reference_atom.wavenumber

        xm, x0, rr = model.calculate_distances()
        nat = len(model.x)

        # Covariant spherical components of the unit separation vector.
        x = np.zeros((3, nat, nat), dtype="complex128")
        x[0] = xm
        x[1] = x0
        x[2] = -np.conj(xm)

        g = np.array([[0, 0, -1], [0, 1, 0], [-1, 0, 0]])

        # D^(E) / HBAR from Eq. (2.23), with DDI = 1 (dipole-dipole interaction is ON)
        # Vacuum self-decay is included separately in the resolvent.
        k_medium = self.wavenumber

        d1 = ((DDI * 1 - 1j * k_medium * rr - (k_medium * rr) ** 2) / ((rr + np.identity(nat)) ** 3)
                * np.exp(1j * k_medium * rr)) * (np.ones(nat) - np.identity(nat))
        d2 = -1 * ((DDI * 3 - 3 * 1j * k_medium * rr - (k_medium * rr) ** 2) / ((rr + np.identity(nat)) ** 3)
                * np.exp(1j * k_medium * rr)) * (np.ones(nat) - np.identity(nat))

        outer1 = np.einsum("ij,ab->ijab", g, d1)
        outer2 = np.einsum("iab,jab->ijab", x, x) * d2
        D = outer1 + outer2

        m = self.medium_atom.m


        # d_down = <e|d|g>, d_up = <g|d|e>
        d_a = np.array([d_down(self.medium_atom, 0, mi) for mi in m])
        d_b = np.array([d_up(self.medium_atom, 0, mi) for mi in m])

        # (atom a, atom b, state m, state n) -> ((a, m), (b, n)).
        di = np.einsum("mi,nj,ijab->abmn", d_a, d_b, D)
        self.sigma = reshape_to_matrix(di)

    def get_resolvent_for_v(self, omega):
        """Return the medium resolvent in inverse-energy units."""
        atom = self.medium_atom
        return np.linalg.inv(HBAR * (omega - atom.omega + 1j * atom.gamma / 2) * np.eye(len(self.sigma)) - self.sigma)

    def get_sigma_outside(self, model: GeneralModel, resolvent):
        """
        Return the reference atom self-energy in energy units;.
        Omega dependence enters through resolvent, evaluated at omega.
        """
        nat = len(model.x)
        atom = model.reference_atom
        m0 = atom.m0
        m = atom.m
        k_medium = self.wavenumber

        xm, x0, rr = model.calculate_distances_to_reference_atom()

        x = np.zeros((3, nat), dtype="complex128")
        x[0] = xm
        x[1] = x0
        x[2] = -np.conj(xm)

        g = np.array([[0, 0, -1], [0, 1, 0], [-1, 0, 0]])

        d1 = ((DDI * 1 - 1j * k_medium * rr - (k_medium * rr) ** 2)
              / (rr ** 3) * np.exp(1j * k_medium * rr))
        d2 = -1 * ((DDI * 3 - 3 * 1j * k_medium * rr - (k_medium * rr) ** 2)
                   / (rr ** 3) * np.exp(1j * k_medium * rr))

        outer1 = np.einsum("ij,a->ija", g, d1)
        outer2 = np.einsum("ia,ja->ija", x, x) * d2
        D = outer1 + outer2

        mV = [-1, 0, 1]
        d_b1 = np.array([d_up(self.medium_atom, 0, mi) for mi in mV])
        d_a2 = np.array([d_down(self.medium_atom, 0, mi) for mi in mV])

        sigma_out = -1j * HBAR * atom.gamma / 2 * np.identity(len(m))

        # Sum independent intermediate signal-ground channels.
        for i in range(len(m0)):
            d_a1 = np.array([d_down(atom, m0[i], mi) for mi in m])
            d_b2 = np.array([d_up(atom, m0[i], mi) for mi in m])

            # di_1: (atom, signal state, medium state).
            # di_2: (atom, medium state, signal state).
            di_1 = np.einsum("mi,nj,ija->amn", d_a1, d_b1, D)
            di_2 = np.einsum("mi,nj,ija->amn", d_a2, d_b2, D)

            # c: medium -> signal; b: signal -> medium.
            c = di_1.transpose(1, 0, 2).reshape(len(m), nat * len(mV))
            b = di_2.reshape(nat * len(mV), len(m))

            sigma_out += c @ resolvent @ b

        return sigma_out
