import numpy as np
from .sigma_matrix import MediumSelfEnergyMatrix
from src.radiative_shift.constants import HBAR
from src.radiative_shift.model import GeneralModel
from src.radiative_shift.tools import dipole_mn, dipole_nm
from src.radiative_shift.tools import blocks_to_matrix


class EinsumVMediumSelfEnergyMatrix(MediumSelfEnergyMatrix):
    """
    Self-energy matrix for a (F0=0, F=1) - atomic medium.

    Matrix elements have units of energy, while the corresponding
    resolvent elements have units of inverse energy.

    The medium basis is ordered by atom. For each atom, the magnetic
    sublevels appear in the order M_e = -1, 0, +1; for example:
    (atom 0, -1), (atom 0, 0), (atom 0, +1),
    (atom 1, -1), (atom 1, 0), (atom 1, +1), ...
    """

    def __init__(self, model: GeneralModel):
        self.medium_atom = model.medium_atom
        if (self.medium_atom.F0, self.medium_atom.F) != (0, 1):
            raise ValueError("The medium must be a V atom with F0=0 and F=1")
        self.wavenumber = model.reference_atom.wavenumber

        X_minus_norm, X_0_norm, R_ab = model.calculate_distances()
        N = len(model.x)

        # Covariant spherical components X_mu/R, normalized by the model.
        X_plus_norm = -np.conj(X_minus_norm)
        X_ab_norm = np.array([X_minus_norm, X_0_norm, X_plus_norm], dtype="complex128")

        spherical_metric = np.array([[0, 0, -1], [0, 1, 0], [-1, 0, 0]])

        # D^(E) / HBAR from Eq. (2.23)
        # Vacuum self-decay is included separately in the resolvent.
        k0 = self.wavenumber

        # D^(E)/HBAR = D_coeff1 * spherical_metric + D_coeff2 * (X_norm outer X_norm).
        D_coeff1 = (
            (1 - 1j * k0 * R_ab - (k0 * R_ab) ** 2)
            / ((R_ab + np.identity(N)) ** 3) * np.exp(1j * k0 * R_ab)
        ) * (np.ones(N) - np.identity(N))
        D_coeff2 = -1 * (
            (3 - 3 * 1j * k0 * R_ab - (k0 * R_ab) ** 2)
            / ((R_ab + np.identity(N)) ** 3) * np.exp(1j * k0 * R_ab)
        ) * (np.ones(N) - np.identity(N))

        metric_term = np.einsum("uv,ab->uvab", spherical_metric, D_coeff1)
        radial_term = np.einsum(
            "uab,vab->uvab", X_ab_norm, X_ab_norm
        ) * D_coeff2
        D_E_over_hbar = metric_term + radial_term

        medium_excited_sublevels = self.medium_atom.m


        # Medium dipoles f_eg = <e|d|g>, f_ge = <g|d|e>, Eq. (2.30).
        f_eg = np.array([dipole_nm(self.medium_atom, 0, M_e) for M_e in medium_excited_sublevels])
        f_ge = np.array([dipole_mn(self.medium_atom, 0, M_e) for M_e in medium_excited_sublevels])

        # (atom a, atom b, state e, state e') -> ((a, e), (b, e')).
        medium_self_energy_blocks = np.einsum("eu,pv,uvab->abep", f_eg, f_ge, D_E_over_hbar)
        self.sigma = blocks_to_matrix(medium_self_energy_blocks)

    def get_medium_resolvent(self, omega):
        """Return the medium resolvent at E = HBAR*omega, in inverse energy."""
        medium_atom = self.medium_atom
        return np.linalg.inv(
            HBAR * (omega - medium_atom.omega + 1j * medium_atom.gamma / 2)
            * np.eye(len(self.sigma)) - self.sigma
        )

    def get_reference_self_energy(self, model: GeneralModel, resolvent):
        """
        Return the dressed reference self-energy of Eqs. (2.13)-(2.14).

        The result is in energy units and includes vacuum decay, Eq. (2.18).
        resolvent is the medium resolvent evaluated at the chosen omega.
        """
        N = len(model.x)
        atom = model.reference_atom
        m0 = atom.m0  # Reference ground sublevels.
        m = atom.m  # Reference excited sublevels.
        k0 = self.wavenumber

        X_minus_norm, X_0_norm, R_0a = model.calculate_distances_to_reference_atom()

        X_plus_norm = -np.conj(X_minus_norm)
        X_0a_norm = np.array([X_minus_norm, X_0_norm, X_plus_norm], dtype="complex128")

        spherical_metric = np.array([[0, 0, -1], [0, 1, 0], [-1, 0, 0]])

        # D^(E)/HBAR = D_coeff1 * spherical_metric + D_coeff2 * (X_norm outer X_norm).
        D_coeff1 = (
            (1 - 1j * k0 * R_0a - (k0 * R_0a) ** 2)
            / (R_0a ** 3) * np.exp(1j * k0 * R_0a)
        )
        D_coeff2 = -1 * (
            (3 - 3 * 1j * k0 * R_0a - (k0 * R_0a) ** 2)
            / (R_0a ** 3) * np.exp(1j * k0 * R_0a)
        )

        metric_term = np.einsum("uv,a->uva", spherical_metric, D_coeff1)
        radial_term = np.einsum("ua,va->uva", X_0a_norm, X_0a_norm) * D_coeff2
        D_E_over_hbar = metric_term + radial_term

        mV = [-1, 0, 1]  # Medium excited sublevels.
        f_ge = np.array([dipole_mn(self.medium_atom, 0, M_e) for M_e in mV])
        f_eg = np.array([dipole_nm(self.medium_atom, 0, M_e) for M_e in mV])

        sigma_ref = (
            -1j * HBAR * atom.gamma / 2 * np.identity(len(m))
        )

        # Sum independent intermediate reference-ground channels.
        for g in range(len(m0)):
            d_nm = np.array([
                dipole_nm(atom, m0[g], M_e)
                for M_e in m
            ])
            d_mn = np.array([
                dipole_mn(atom, m0[g], M_e)
                for M_e in m
            ])

            # Axes: (atom a, reference state e, medium state e'), and its reverse.
            dc = np.einsum("eu,pv,uva->aep", d_nm, f_ge, D_E_over_hbar)
            db = np.einsum("pu,ev,uva->ape", f_eg, d_mn, D_E_over_hbar)

            # dc: medium -> reference (2.29); db: reference -> medium (2.27).
            dc_matrix = dc.transpose(1, 0, 2).reshape(
                len(m), N * len(mV)
            )
            db_matrix = db.reshape(
                N * len(mV), len(m)
            )

            sigma_ref += dc_matrix @ resolvent @ db_matrix

        return sigma_ref
