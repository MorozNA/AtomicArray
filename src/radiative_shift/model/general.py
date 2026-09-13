import numpy as np
from abc import ABC, abstractmethod
import logging
from src.radiative_shift.atomspecies import AtomSpecies


class Properties:
    width: float
    length: float
    density: float
    n_atoms: int


class GeneralModel(ABC):
    def __init__(self):
        # Coordinates are in cm; angles are in radians.
        self.x = np.array([], dtype=float)
        self.y = np.array([], dtype=float)
        self.z = np.array([], dtype=float)
        self.medium_atom = None
        self.properties = Properties()
        self.reference_atom = None
        self.x_ref = np.array([])
        self.y_ref = np.array([])
        self.z_ref = np.array([])


    def set_reference_atom_cartesian(self, atom_type:AtomSpecies, x, y, z):
        x_ref = np.array([x], dtype=float)
        y_ref = np.array([y], dtype=float)
        z_ref = np.array([z], dtype=float)
        self.reference_atom = atom_type
        self.x_ref, self.y_ref, self.z_ref = x_ref, y_ref, z_ref


    def set_reference_position_cylindrical(self, atom_type:AtomSpecies, r=0, phi=0, z=None):
        x = r * np.cos(phi)
        y = r * np.sin(phi)
        if z is None:
            z = (np.amin(self.z) + np.amax(self.z)) / 2 if len(self.z) else 0
        self.set_reference_atom_cartesian(atom_type, x, y, z)


    def add_medium_atom_cylindrical(self, r, phi=0, z=None):
        x = r * np.cos(phi)
        y = r * np.sin(phi)
        if z is None:
            z = (np.amin(self.z) + np.amax(self.z)) / 2
        self.x = np.concatenate([self.x, np.array([x])])
        self.y = np.concatenate([self.y, np.array([y])])
        self.z = np.concatenate([self.z, np.array([z])])
        self._refresh_properties()

    def add_atom_xyz(self, x, y, z):
        self.x = np.concatenate([self.x, np.array([x])])
        self.y = np.concatenate([self.y, np.array([y])])
        self.z = np.concatenate([self.z, np.array([z])])
        self._refresh_properties()


    def set_medium_atom_i_position(self, r, phi=0, pos=-1):
        self.x[pos], self.y[pos] = r * np.cos(phi), r * np.sin(phi)
        self._refresh_properties()

    def replace_atom_z(self, z, pos=-1):
        self.z[pos] = z
        self._refresh_properties()


    def rotate_about_z(self, phi):
        c, s = np.cos(phi), np.sin(phi)
        self.x, self.y = self.x * c - self.y * s, self.x * s + self.y * c
        self.x_ref, self.y_ref = self.x_ref * c - self.y_ref * s, self.x_ref * s + self.y_ref * c
        self._refresh_properties()


    def rotate_about_y(self, angle):
        c, s = np.cos(angle), np.sin(angle)
        self.x, self.z = c * self.x + s * self.z, -s * self.x + c * self.z
        self.x_ref, self.z_ref = c * self.x_ref + s * self.z_ref, -s * self.x_ref + c * self.z_ref
        self._refresh_properties()


    def calculate_distances(self):
        # Calculates covariant (x_{-1}, x_{0}) distances between atoms i and j
        distance_xm = (self.x[:, None] - self.x - 1j * self.y[:, None] + 1j * self.y) / np.sqrt(2)  # covariant
        distance_z = self.z[:, None] - self.z
        arr_r = np.sqrt(2 * np.square(np.abs(distance_xm)) + np.square(distance_z))
        if np.any((arr_r == 0) & ~np.eye(len(self.x), dtype=bool)):
            raise ValueError("Distinct medium atoms must not coincide")
        safe_r = arr_r + np.identity(len(self.x))  # avoid division by zero only on the diagonal
        return distance_xm / safe_r, distance_z / safe_r, arr_r

    def calculate_distances_to_reference_atom(self):
        if len(self.x_ref) != 1:
            raise ValueError("Add a reference atom before calculating its self-energy")
        distance_xm = (self.x_ref - self.x - 1j * self.y_ref + 1j * self.y) / np.sqrt(2)
        distance_z = self.z_ref - self.z
        arr_r = np.sqrt(2 * np.square(np.abs(distance_xm)) + np.square(distance_z))
        if np.any(arr_r == 0):
            raise ValueError("The reference atom must not coincide with a medium atom")
        return distance_xm / arr_r, distance_z / arr_r, arr_r

    def remove_near_duplicates(self, r_lbar=0.1):
        r = r_lbar * self.medium_atom.lbar
        excess = []
        for i in range(len(self.x)):
            for j in range(i):
                if j not in excess and (self.x[i] - self.x[j]) ** 2 \
                        + (self.y[i] - self.y[j]) ** 2 + (self.z[i] - self.z[j]) ** 2 < r ** 2:
                    excess.append(i)
                    break
        self.x = np.delete(self.x, excess)
        self.y = np.delete(self.y, excess)
        self.z = np.delete(self.z, excess)
        self._refresh_properties()

    @abstractmethod
    def _refresh_properties(self):
        # Use medium_atom.lbar (cm) to normalize lengths and density.
        # TODO: add parameters: length_parameter, width_parametr, density_parameter
        pass

    def write_log(self):
        logging.info("=========================================================")
        logging.info("Model parameters")
        logging.info(f"Length = {self.properties.length:.2f} λ / 2 π")
        logging.info(f"Width = {self.properties.width:.2f} λ / 2 π")
        logging.info(f"Number of atoms {self.properties.n_atoms}")
        logging.info(f"Density {self.properties.density:.2f} nλbar^3")
