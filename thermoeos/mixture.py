"""
Mixture
=======

Provides
--------
    * Classes that represents studied mixture or its properties.
"""

from .equations import Phase
import numpy as np
from scipy.constants import R
from copy import deepcopy


class Mixture(object):
    """
    Contains data on the studied mixture.
    Attributes
    ----------
    name : str
        The name of the mixture.
    x : list of float
        A list of mole fractions of the inidividual components in liquid phase.
        The list mast exclude the last compoenent.
    y : list of float
        A list of mole fractions of the inidividual components in vapor phase.
        The list mast exclude the last compoenent.
    Tb : list of float [K]
        Boiling temperatures of individual components.
    Tc : list of float [K]
        Critical temperatures of individual components.
    Pc : list of float [Pa]
        Critical pressures of individual components.
    omega : list of float
        Pitzer acentric factors of individual components.
    phase : Phase, optional
        A phase of the mixture which has constanat mole fraction of
        the components.
    k : 2d array of float, optional
        Binary interaction parameters of components pairs.
        If mixture contains three components(A, B and C) k will look like
        following
          | A | B | C
        --+---+---+---
        A | 0 |kAB|kAC
        --+---+---+---
        B |kAB| 0 |kBC
        --+---+---+---
        C |kAC|kBC| 0
    """

    def __init__(self, name, x, y, Tb, Tc, Pc, omega,
                 phase=Phase.LIQUID, k=None):
        x = np.append(np.array(x), 1 - np.sum(x))
        y = np.append(np.array(y), 1 - np.sum(y))
        self.name = name
        self.x = np.array(x, ndmin=2)
        self.y = np.array(y, ndmin=2)
        self.Tb = np.array(Tb, ndmin=2)
        self.Tc = np.array(Tc, ndmin=2)
        self.Pc = np.array(Pc, ndmin=2)
        self.omega = np.array(omega, ndmin=2)
        self.phase = phase
        if k is None:
            nd = np.max(self.x.shape)
            self.k = np.zeros((nd, nd))
        else:
            self.k = np.array(self.k, ndmin=2)
        # Precalculated for speed
        self.RTc = R * self.Tc
        self.RTc2 = self.RTc**2
        self.RTc25 = R**2 * self.Tc**2.5

    @staticmethod
    def from_thermo_mix(mix):
        name = "-".join(mix.components)
        mmix = Mixture(
            name, mix.zs[:-1],
            mix.zs[:-1],
            mix.Tbs,
            mix.Tcs,
            mix.Pcs,
            mix.omegas
        )
        return mmix

    phi_vap = None
    vol_vap = None
    phi_liq = None
    vol_liq = None
    T = None
    P = None
    s = None

    @property
    def k(self):
        return self._k

    @k.setter
    def k(self, k):
        self._k = np.array(k)

    @property
    def x(self):
        return self._x

    @x.setter
    def x(self, x):
        if isinstance(x, np.ndarray):
            if len(x.shape) < 2:
                x.shape = (1, len(x))
            self._x = x
        else:
            self._x = np.array(x, ndmin=2)
        return None

    @property
    def y(self):
        return self._y

    @y.setter
    def y(self, y):
        if isinstance(y, np.ndarray):
            if len(y.shape) < 2:
                y.shape = (1, len(y))
            self._y = y
        else:
            self._y = np.array(y, ndmin=2)
        return None

    def set_x(self, x, inplace=False):
        if inplace:
            self.x = x
            return None
        else:
            mix = deepcopy(self)
            mix.x = x
            return mix

    def set_y(self, y, inplace=False):
        if inplace:
            self.y = y
            return None
        else:
            mix = deepcopy(self)
            mix.y = y
            return mix

    def set_k(self, k, inplace=False):
        if inplace:
            self.k = k
        else:
            mix = deepcopy(self)
            mix.k = k
            return mix

    def __str__(self):
        s = f"Mix: {self.name}\n"
        s += f"x: {self.x[0]}\n"
        s += f"y: {self.y[0]}\n"
        s += f"s: {self.s}\n"
        return s
