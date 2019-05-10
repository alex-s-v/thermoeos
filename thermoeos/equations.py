"""
Equations
=========

Provides
--------
    * Various equations of state.
    * Classes that represents different states of the studied mixtures.
"""
import numpy as np
from scipy.constants import R
from enum import Enum
from copy import deepcopy
from abc import ABC

# TODO: re-investigate difference between constant mole fraction and
# changing with iterations


class Predict(Enum):
    PRESSURE = 1
    TEMPERATURE = 2


class Phase(Enum):
    LIQUID = 1
    VAPOR = 2


class Equation(ABC):
    """
    Abstract class that represents general equation of state behavior.
    """

    def __init__(self, mix):
        self.mix = deepcopy(mix)

    def _calculate_parameters(self, T, P, phase):
        raise NotImplementedError

    def _calculate_coefficients(self, RT, P, A, B, A_, B_):
        raise NotImplementedError

    def _calculate_hs(self, Z, V, RT, A, B, A_, B_):
        raise NotImplementedError

    def _calculate_compressibility_factor(self, coefs):
        """
        Calculate compressibility factor from the cubic equation with
        given coefficients.
        Parameters
        ----------
        coefs : list of float
            Coefficients of the cubic equation a + b*x + c*x^2 + d*x^3 = 0,
            where coefs = [a, b, c, d]
        Returns
        -------
        z : np.array of float
            Compressibility factors for different states of the mixture.
        : bool
            True if all factors are real numbers, False otherwise.
        """
        z = np.roots(coefs)
        if all(np.isreal(z)):
            return z, True
        return z, False

    def _calculate_fugacity_coefficient(self, T, P, phase):
        """
        Calculate fugacity coefficient and molar volume for the given phase
        of the mixture.
        Parameters
        ----------
        T : float
            A temperature of the mixture.
        P : float
            A pressure of the mixture.
        phase : Phase
            A phase for which calculations will be done.
        Returns
        -------
        phi : np.array of float
            Fugacity coefficients for each of the components in
            the given phase.
        V : float
            A molar volume of the mixture in the given phase.
        calculated : bool
            True if calculations succeeded, False otherwise.
        """
        A, B, A_, B_ = self._calculate_parameters(T, P, phase)
        RT = R * T
        coefs = self._calculate_coefficients(RT, P, A, B, A_, B_)
        zs, calculated = self._calculate_compressibility_factor(coefs)
        if not calculated:
            return [], 0, False
        Z = min(zs) if phase == Phase.LIQUID else max(zs)
        V = Z * RT / P
        h1, h2, h3, h4 = self._calculate_hs(Z, V, RT, A, B, A_, B_)
        phi = np.exp(h1 - np.log(h2) + h3 * np.log(h4))
        return phi, V, calculated

    def _recalculate_mole_fraction_liquid(self, ep):
        """
        Recalculate mole fraction of the components in the vapor phase.
        Parameters
        ----------
        ep : Mixture
            A data about equilibrium point for which recalculation
            will be done.
        """
        ep = deepcopy(ep)
        kmf = np.multiply(np.divide(ep.phi_liq, ep.phi_vap), ep.x)
        ep.s = np.sum(kmf)
        ep.y = kmf / ep.s
        return ep

    def _recalculate_mole_fraction_vapor(self, ep):
        """
        Recalculate mole fraction of the components in the liquid phase.
        Parameters
        ----------
        ep : Mixture
            A data about equilibrium point for which recalculation
            will be done.
        """
        ep = deepcopy(ep)
        kmf = np.multiply(np.divide(ep.phi_vap, ep.phi_liq), ep.y)
        ep.s = np.sum(kmf)
        ep.x = kmf / ep.s
        return ep

    def _rmf(self, ep):
        """
        Recalculate mole fraction of the components.
        Parameters
        ----------
        ep : Mixture
            A data about equilibrium point for which recalculation
            will be done.
        """
        if self.mix.phase.value == Phase.LIQUID.value:
            ep = self._recalculate_mole_fraction_liquid(ep)
        else:
            ep = self._recalculate_mole_fraction_vapor(ep)
        return ep

    def solve(self, T, P):
        """
        Solves EOS and reutrns data about specified point.
        Parameters
        ----------
        T : float [K]
            The temperature of the mixture.
        P : float [Pa]
            The pressure of the mixture.
        Returns
        -------
        Mixture
            Contains data about specified point.
        bool
            True if calculations succeeded, False otherwise.
        """
        pv, vv, c1 = self._calculate_fugacity_coefficient(T, P, Phase.VAPOR)
        pl, vl, c2 = self._calculate_fugacity_coefficient(T, P, Phase.LIQUID)
        ep = deepcopy(self.mix)
        ep.T = T
        ep.P = P
        ep.phi_vap = pv
        ep.phi_liq = pl
        ep.vol_vap = vv
        ep.vol_liq = vl
        if c1 and c2:
            ep = self._rmf(ep)
        return ep, (c1 and c2)


class SRK(Equation):
    """
    Represents logic of Soave-Redlich-Kwang EOS.
    """

    def __init__(self, mix):
        self.m = 0.48508 + 1.55171 * mix.omega - 0.15613 * mix.omega**2
        self.ca = 0.42747
        self.cb = 0.08664
        self.b_ = self.cb * mix.RTc / mix.Pc
        self.a_ = self.ca * mix.RTc2 / mix.Pc
        super(SRK, self).__init__(mix)

    def _calculate_parameters(self, T, P, phase):
        """
        Calculate four main parameters.
        Parameres
        ---------
        T : float [K]
            The temperature of the mixture.
        P : float [Pa]
            The pressure of the mixture.
        phase : Phase
            A phase of the mixture.
        Returns
        -------
        Returns parameters in the following order: A, B, A_, B_,
        where A and B are single `float` values, A_ and B_
        are lists of `float`.
        """
        frac = self.mix.x if phase == Phase.LIQUID else self.mix.y

        alpha = (1 + self.m * (1 - np.sqrt(T / self.mix.Tc)))**2
        a_ = self.a_ * alpha

        aa = (1 - self.mix.k) * np.sqrt(a_.T @ a_)
        top = 2 * (aa @ frac.T).T
        a = np.sum(aa * (frac.T @ frac))
        b = np.sum(frac * self.b_)

        RT = R * T
        A_ = top / a
        B_ = self.b_ * P / RT

        A = a * P / RT**2
        B = b * P / RT

        return A, B, A_, B_

    def _calculate_coefficients(self, RT, P, A, B, A_, B_):
        """
        Calculates coefficients for cubic form of the EOS.
        """
        return [1, -1, A - B - B**2, -A * B]

    def _calculate_hs(self, Z, V, RT, A, B, A_, B_):
        """
        Calculates parts of the fugacity equation.
        """
        h1 = (Z - 1) * B_ / B
        h2 = Z - B
        h3 = (A / B) * (B_ / B - A_)
        h4 = 1 + B / Z

        return h1, h2, h3, h4


class PR(Equation):
    def __init__(self, mix):
        self.m = 0.37464 + 1.54226 * mix.omega - 0.26992 * mix.omega**2
        self.ca = 0.457235
        self.cb = 0.077796
        self.b_ = self.cb * mix.RTc / mix.Pc
        self.a_ = self.ca * mix.RTc2 / mix.Pc
        super(PR, self).__init__(mix)

    def _calculate_parameters(self, T, P, phase):
        frac = self.mix.x if phase == Phase.LIQUID else self.mix.y

        alpha = (1 + self.m * (1 - np.sqrt(T / self.mix.Tc)))**2
        a_ = self.a_ * alpha

        aa = (1 - self.mix.k) * np.sqrt(a_.T @ a_)
        top = 2 * (aa @ frac.T).T
        a = np.sum(aa * (frac.T @ frac))
        b = np.sum(frac * self.b_)

        RT = R * T
        A_ = top / a
        B_ = self.b_ * P / RT

        A = a * P / RT**2
        B = b * P / RT

        return A, B, A_, B_

    def _calculate_coefficients(self, RT, P, A, B, A_, B_):
        return [1, -1 + B, A - 2 * B - 3 * B**2, -A * B + B**2 + B**3]

    def _calculate_hs(self, Z, V, RT, A, B, A_, B_):
        h1 = B_ * (Z - 1) / B
        h2 = Z - B
        h3 = (A / (2.8284 * B)) * (B_ / B - A_)
        h4 = (Z + 2.4142 * B) / (Z - 0.4142 * B)
        return h1, h2, h3, h4


class RK(Equation):
    def __init__(self, mix):
        self.ca = 0.42747
        self.cb = 0.08664
        self.a_ = self.ca * mix.RTc25 / mix.Pc
        self.b_ = self.cb * mix.RTc / mix.Pc
        super(RK, self).__init__(mix)

    def _calculate_parameters(self, T, P, phase):
        frac = self.mix.x if phase == Phase.LIQUID else self.mix.y

        A_ = self.a_ * P / (R**2 * T**2.5)
        B_ = self.b_ * P / (R * T)

        A = np.sum(frac * np.sqrt(A_))**2
        B = np.sum(frac * B_)

        return (A, B, A_, B_)

    def _calculate_coefficients(self, RT, P, A, B, A_, B_):
        return [1, -1, A - B - B**2, -A * B]

    def _calculate_hs(self, Z, V, RT, A, B, A_, B_):
        h1 = (Z - 1) * B_ / B
        h2 = Z - B
        h3 = (A / B) * (B_ / B - 2 * np.sqrt(A_ / A))
        h4 = 1 + B / Z
        return h1, h2, h3, h4


class VdW(Equation):
    def __init__(self, mix):
        self.ca = 0.421875
        self.cb = 0.125000
        self.A_ = self.ca * mix.RTc2 / mix.Pc
        self.B_ = self.cb * mix.RTc / mix.Pc
        super(VdW, self).__init__(mix)

    def _calculate_parameters(self, T, P, phase):
        frac = self.mix.x if phase == Phase.LIQUID else self.mix.y

        A = np.sum(frac * np.sqrt(self.A_))**2
        B = np.sum(frac * self.B_)

        return (A, B, self.A_, self.B_)

    def _calculate_coefficients(self, RT, P, A, B, A_, B_):
        return [1, -(B * P / RT + 1), A * P / RT**2, -A * B * P**2 / RT**3]

    def _calculate_hs(self, Z, V, RT, A, B, A_, B_):
        h1 = B / (V - B)
        h2 = Z * (1 - B / V)
        h3 = -2 * np.sqrt(A_ * A) / (RT * V)
        h4 = np.exp(1)
        return h1, h2, h3, h4
