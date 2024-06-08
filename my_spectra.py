""" Analyse a quasar spectra """

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.transforms as transforms
from astropy.convolution import Box1DKernel, convolve
import linelist


class Quasar:
    """ This class reads the spectra """

    def __init__(self, lamb, flux, eflux):

        self.lamb = lamb
        self.flux = flux
        self.eflux = eflux

    @classmethod
    def from_file(cls, filename, columns=[0, 1, 2]):
        lamb, flux, eflux = np.loadtxt(
            filename, usecols=columns, unpack=True)
        return cls(lamb, flux, eflux)

    def tick_plot(self, test_z):
        """ Given a redshift, plots tick for some ion transitions """

        my_lines = ['HI_1215', 'HI_1025', 'HI_972', 'HI_949', 'HI_937',
                    'HI_930', 'HI_926', 'HI_923', 'HI_920', 'HI_919',
                    'HI_918', 'HI_917', 'HI_916', 'HI_915a', 'HI_912a',
                    'CIV_1548', 'CIV_1550', 'SiII_1526', 'MgII_2803',
                    'MgII_2796', 'CoII_1480', 'CaII_K', 'CaII_H',
                    'AlII_1670', 'AlIII_1854', 'CoII_1552', 'FeII_2586',
                    'FeII_1143',
                    'FeII_1260', 'FeII_1608', 'FeII_2260', 'FeII_2600',
                    'FeII_2382']

        # Check if plot is open, else do it
        ax = plt.gca()
        ax.clear()

        # Tick plots
        ax.step(self.lamb, self.flux)

        for i in range(len(my_lines)):
            x = (1 + test_z)*linelist.wl[my_lines[i]]

            if (min(self.lamb) <= x) & (max(self.lamb) >= x):
                ax.axvline(x, ymin=0, ymax=1, ls="dashed", lw=.5, color="r")
                trans = transforms.blended_transform_factory(
                    ax.transData, ax.transAxes)
                ax.text(x, 0.5,  my_lines[i].split('_')[0] +
                        ' ' + my_lines[i].split('_')[1], transform=trans,
                        rotation=90, fontsize=10, color="r")

        ax.set_xlabel(r"$\rm Wavelength$")
        ax.set_ylabel(r"$\rm Flux$")
        plt.show()

    def smooth_spec(self, width):
        ''' smoothes the spectral flux by certain width'''

        kernel1D = Box1DKernel(width=width)
        smooth_flux = convolve(self.flux, kernel1D)
        return smooth_flux
