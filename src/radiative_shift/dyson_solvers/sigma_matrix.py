import numpy as np
from src.radiative_shift.constants import HBAR
from src.radiative_shift.model import GeneralModel
from src.radiative_shift.tools import dipole_mn, dipole_nm
from src.radiative_shift.tools import matrix_to_blocks, blocks_to_matrix
from abc import ABC


class MediumSelfEnergyMatrix(ABC):
    """Medium pair self-energy Sigma^(ab), Eq. (2.30), in energy units."""

    sigma: np.array


# TODO: get rid of the class, only leave functions
class VMediumSelfEnergyMatrix(MediumSelfEnergyMatrix):
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
        X_ab_norm = np.array([X_minus_norm, X_0_norm, X_plus_norm], dtype='complex128')

        spherical_metric = np.array([[0, 0, -1], [0, 1, 0], [-1, 0, 0]])

        # D^(E) / HBAR from Eq. (2.23)
        # Vacuum self-decay is included separately in the resolvent.
        k0 = self.wavenumber

        # D^(E)/HBAR = D_coeff1 * spherical_metric + D_coeff2 * (X_norm outer X_norm).
        D_coeff1 = ((1 - 1j * k0 * R_ab - (k0 * R_ab) ** 2)
                    / ((R_ab + np.identity(N)) ** 3) * np.exp(1j * k0 * R_ab)) * (np.ones(N) - np.identity(N))
        D_coeff2 = -1 * ((3 - 3 * 1j * k0 * R_ab - (k0 * R_ab) ** 2)
                         / ((R_ab + np.identity(N)) ** 3) * np.exp(1j * k0 * R_ab)) * (np.ones(N) - np.identity(N))

        medium_self_energy_blocks = np.zeros((N, N, 3, 3), dtype=complex)

        medium_excited_sublevels = self.medium_atom.m

        # Medium dipoles f_eg = <e|d|g>, f_ge = <g|d|e>, Eq. (2.30).
        f_eg = np.array([dipole_nm(self.medium_atom, 0, M_e) for M_e in medium_excited_sublevels])
        f_ge = np.array([dipole_mn(self.medium_atom, 0, M_e) for M_e in medium_excited_sublevels])

        for e in range(len(medium_excited_sublevels)):
            for e_prime in range(len(medium_excited_sublevels)):
                medium_self_energy_blocks[:, :, e, e_prime] = (
                    np.dot(spherical_metric @ f_eg[e], f_ge[e_prime])
                    * D_coeff1
                )
                medium_dipole_dyad = np.outer(f_eg[e], f_ge[e_prime])
                for mu in range(3):
                    for nu in range(3):
                        medium_self_energy_blocks[:, :, e, e_prime] += (
                            medium_dipole_dyad[mu, nu] * X_ab_norm[mu] * X_ab_norm[nu]
                            * D_coeff2
                        )

        # (atom a, atom b, state e, state e') -> ((a, e), (b, e')).
        self.sigma = blocks_to_matrix(medium_self_energy_blocks)


    def get_medium_resolvent(self, omega):
        """Return the medium resolvent at E = HBAR*omega, in inverse energy."""
        medium_atom = self.medium_atom
        return np.linalg.inv(HBAR * (omega - medium_atom.omega + 1j * medium_atom.gamma / 2) * np.eye(len(self.sigma)) - self.sigma)


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
        mV = [-1, 0, 1]  # Medium excited sublevels.
        k0 = self.wavenumber

        X_minus_norm, X_0_norm, R_0a = model.calculate_distances_to_reference_atom()
        X_plus_norm = -np.conj(X_minus_norm)
        X_0a_norm = np.array(
            [X_minus_norm, X_0_norm, X_plus_norm], dtype='complex128'
        )

        spherical_metric = np.array([[0, 0, -1], [0, 1, 0], [-1, 0, 0]])

        # D^(E)/HBAR = D_coeff1 * spherical_metric
        #             + D_coeff2 * (X_norm outer X_norm).
        D_coeff1 = (
                (1 - 1j * k0 * R_0a - (k0 * R_0a) ** 2)
                / (R_0a ** 3) * np.exp(1j * k0 * R_0a)
        )
        D_coeff2 = -1 * (
                (3 - 3 * 1j * k0 * R_0a - (k0 * R_0a) ** 2)
                / (R_0a ** 3) * np.exp(1j * k0 * R_0a)
        )

        # db: reference -> medium, Eq. (2.27).
        # dc: medium -> reference, Eq. (2.29).
        db = matrix_to_blocks(
            np.zeros((N * len(mV), len(m)), dtype=complex),
            len(mV), len(m),
        )
        dc = matrix_to_blocks(
            np.zeros((len(m), N * len(mV)), dtype=complex),
            len(m), len(mV),
        )

        d_nm = np.zeros((len(m0), len(m), 3), dtype=complex)
        d_mn = np.zeros((len(m0), len(m), 3), dtype=complex)

        f_ge = np.array([dipole_mn(self.medium_atom, 0, M_e) for M_e in mV])
        f_eg = np.array([dipole_nm(self.medium_atom, 0, M_e) for M_e in mV])

        sigma_ref = -1j * HBAR * atom.gamma / 2 * np.identity(len(m))

        for g in range(len(m0)):
            for e in range(len(m)):
                d_nm[g, e] = dipole_nm(atom, m0[g], m[e])
                d_mn[g, e] = dipole_mn(atom, m0[g], m[e])

                for e_prime in range(len(mV)):
                    dc[:, :, e, e_prime] = (
                            np.dot(spherical_metric @ d_nm[g, e], f_ge[e_prime])
                            * D_coeff1
                    )
                    db[:, 0, e_prime, e] = (
                            np.dot(spherical_metric @ f_eg[e_prime], d_mn[g, e])
                            * D_coeff1
                    )

                    outer_c = np.outer(d_nm[g, e], f_ge[e_prime])
                    outer_b = np.outer(f_eg[e_prime], d_mn[g, e])

                    for mu in range(3):
                        for nu in range(3):
                            dc[:, :, e, e_prime] += (
                                    outer_c[mu, nu] * X_0a_norm[mu]
                                    * X_0a_norm[nu] * D_coeff2
                            )
                            db[:, 0, e_prime, e] += (
                                    outer_b[mu, nu] * X_0a_norm[mu]
                                    * X_0a_norm[nu] * D_coeff2
                            )

            sigma_ref += blocks_to_matrix(dc) @ resolvent @ blocks_to_matrix(db)

        return sigma_ref
