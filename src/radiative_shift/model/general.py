import numpy as np
from .param import LBAR
from abc import ABC
import logging


class Properties:
    radius: float
    length: float
    density: float
    noa: float


class GeneralModel(ABC):
    def __init__(self):
        self.x = np.array([])
        self.y = np.array([])
        self.z = np.array([])
        self.properties = Properties()

    def add_atom(self, r, phi=0):
        x = r * np.cos(phi)
        y = r * np.sin(phi)
        z = np.amax(self.z) / 2
        self.x = np.concatenate([self.x, np.array([x])])
        self.y = np.concatenate([self.y, np.array([y])])
        self.z = np.concatenate([self.z, np.array([z])])

    def replace_atom(self, r, phi=0):
        self.x[-1], self.y[-1] = r * np.cos(phi), r * np.sin(phi)

    def rotate_model(self, phi):
        self.x, self.y = self.x * np.cos(phi) - self.y * np.sin(phi), self.x * np.sin(phi) + self.y * np.cos(phi)

    def calculate_distances(self):
        distance_xm = (self.x[:, None] - self.x - 1j * self.y[:, None] + 1j * self.y) / np.sqrt(2)  # covariant
        distance_z = self.z[:, None] - self.z
        arr_r = np.sqrt(2 * np.square(np.abs(distance_xm)) + np.square(distance_z)) + np.identity(
            len(self.x))  # ignoring division by zero
        return distance_xm / arr_r, distance_z / arr_r, arr_r

    def calculate_distances_to_atom(self, atom_index=-1):
        distance_xm = (self.x[atom_index] - self.x - 1j * self.y[atom_index] + 1j * self.y) / np.sqrt(2)
        distance_z = self.z[atom_index] - self.z
        arr_r = np.sqrt(2 * np.square(np.abs(distance_xm)) + np.square(distance_z))
        arr_r[atom_index] = 1
        return distance_xm / arr_r, distance_z / arr_r, arr_r

    def remove_duplicates(self, r=0.1 * LBAR):
        excess = []
        for i in range(len(self.x)):
            for j in range(len(self.x)):
                if i != j and (self.x[i] - self.x[j]) ** 2 \
                        + (self.y[i] - self.y[j]) ** 2 + (self.z[i] - self.z[j]) ** 2 < r ** 2:
                    excess.append(i)
                    continue
        self.x = np.delete(self.x, excess)
        self.y = np.delete(self.y, excess)
        self.z = np.delete(self.z, excess)

    def measure_properties(self):
        self.properties = Properties()
        self.properties.radius = np.amax(np.sqrt([x ** 2 + y ** 2 for x, y in zip(self.x, self.y)])) / LBAR
        self.properties.length = np.amax(self.z) / LBAR
        self.properties.noa = len(self.x)
        self.properties.density = len(self.x) / self.properties.length / self.properties.radius ** 2 / np.pi

    def write_log(self):
        logging.info("=========================================================")
        logging.info("Cylinder parameters")
        logging.info(f"Length = {self.properties.length:.2f} λ / 2 π")
        logging.info(f"Radius = {self.properties.radius:.2f} λ / 2 π")
        logging.info(f"Number of atoms {self.properties.noa}")
        logging.info(f"Density {self.properties.density:.2f} nλbar^3")
