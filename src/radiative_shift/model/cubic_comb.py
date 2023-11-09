import numpy as np
from .general import GeneralModel

class CubicComb(GeneralModel):

    # TODO: parameters should be changed to -> number of etches, period, density
    def __init__(self, length, period, density):
        super().__init__()
        a = density ** (- 1 / 3)
        b = 0.5 * period
        r = b % a / int(b / a)
        a = a + r
        unitsize = a

        width = 0.5 * period
        height = 0.5 * period
        length_etched = 1.5 * period

        # Create the cubic lattice
        nz = round(length / unitsize)
        ny = round(height / unitsize)
        zc = np.arange(0, nz + 1) * unitsize  # because nz * unitsize = length
        yc = np.ones(nz + 1) * 0
        zc_temp = zc
        for i in range(1, ny + 1):
            zc = np.concatenate((zc, zc_temp))
            yc = np.concatenate((yc, np.ones(nz + 1) * i * unitsize))
        self.z = zc
        self.y = yc
        nx = round(width / unitsize)
        self.x = np.zeros(len(self.z))
        z_temp = self.z
        y_temp = self.y
        for i in range(1, nx + 1):
            self.x = np.concatenate((self.x, np.ones(len(z_temp)) * i * unitsize))
            self.z = np.concatenate((self.z, z_temp))
            self.y = np.concatenate((self.y, y_temp))

        '''
        added to test borders
        '''
        # minus_shift = np.amax(self.z) / 2
        # self.z = self.z + minus_shift

        # This part of the code if for etches
        num_etched = round(length / period)

        nx_etched = round(length_etched / unitsize)
        ny_etched = round(height / unitsize)
        xc_echted = np.arange(1, nx_etched + 1) * unitsize + width  # from 1 to ignore duplicates
        yc_ethced = np.zeros(len(xc_echted))
        xc_temp = xc_echted
        for i in range(1, ny_etched + 1):
            xc_echted = np.concatenate((xc_echted, xc_temp))
            yc_ethced = np.concatenate((yc_ethced, np.ones(len(xc_temp)) * i * unitsize))

        '''
        shift changed to change borders
        '''
        add_shift = (length / period - round(length / period)) * period / 2
        shift = abs(period - width) / 2 + add_shift

        nz_etched = int(width / unitsize)
        for i in range(num_etched):
            for j in range(0, nz_etched + 1):
                shift_i = np.ones(len(xc_echted)) * (i * period + shift)
                self.z = np.concatenate((self.z, np.ones(len(xc_echted)) * j * unitsize + shift_i))
                self.x = np.concatenate((self.x, xc_echted))
                self.y = np.concatenate((self.y, yc_ethced))

        self.y = self.y - height / 2
        self.x = -self.x

        # self.remove_duplicates(unitsize / 5)
        self.measure_properties()
        self.write_log()

    def measure_properties(self):
        # TODO: maybe add another property: height
        self.properties.radius = np.amax(abs(self.x))
        self.properties.length = np.amax(self.z)  # / LBAR
        self.properties.noa = len(self.x)
        self.properties.density = len(self.x) / self.properties.length / self.properties.radius ** 2 / np.pi
