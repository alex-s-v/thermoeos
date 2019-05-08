"""
Equilibrium
===========

Provides
--------
    * Algorithm for calculating vapor-liquid phase equilibrium.
    * Algorithm for calculating binary interaction parameters.
"""
import numpy as np
from scipy import optimize
from .equations import Phase, Predict, SRK, PR
from .mixture import Mixture


# TODO: improve algorithm that finds borders
def initialise(f, init, step_size, max_steps):
    """
    Find initial bracketing interval for optimization algorithm.
    Parameters
    ----------
    f : Function
        A function for which bracketing interval needs to be found.
    init : float
        A starting point for the iterative algorithm.
    step_size : float
        A size of the step of the iterative algorithm.
    max_size : float
        A maximum number of steps of the iterative algorithm.
    Returns
    -------
    float
        A left value of a bracketing interval.
    float
        A right value of a bracketing interval.
    """
    a = 1e-10
    b = 1e-10
    av = 0
    bv = 0
    for _ in range(max_steps):
        ep, calculated = f(init)
        if calculated:
            s = ep.s - 1
            if a > s:
                a = s
                av = init
            elif b < s:
                b = s
                bv = init
        init += step_size
        if np.sign(a) != np.sign(b):
            return av, bv
    return None


def find_equilibrium(eq, T, P, predict=Predict.PRESSURE):
    """
    Calculate an equilibrium point using the specified equation of state.
    Parameters
    ----------
    eq : Equation
        An equation of state to calculate an equilibrium point.
    T : float
        A temperature of the mixture.
    P : float
        A pressure of the mixture.
    predict : Predict, optional
        Which parameter to predict, a pressure or a temperature.
    Returns
    -------
    Mixture
        An equilibrium point with predicted properties.
    bool
        True if equilibrium point is found, False otherwise.
    """
    if predict == Predict.PRESSURE:
        step = 1e4
        max_steps = 100
        init = P

        def solve(val): return eq.solve(T, val)
    elif predict == Predict.TEMPERATURE:
        step = 5
        max_steps = 15
        init = T

        def solve(val): return eq.solve(val, P)
    else:
        raise Exception(f"Unknown value to predict: {predict}")

    bounds = initialise(solve, init, step, max_steps)

    def to_opt(val):
        ep, calculated = solve(val)
        if calculated:
            return ep.s - 1
        else:
            return None
    try:
        x, r = optimize.brentq(to_opt, min(bounds), max(
            bounds), xtol=1e-3, rtol=1e-3, disp=False, full_output=True)
    except:
        return None, False

    if r.converged:
        return solve(x)
    else:
        return None, False


def optim_function(k, eqs, predict, exp, ws):
    """
    Optimization function for the bip search algotithm.
    Parameters
    ----------
    k : float
        Binary interaction parameter.
    eqs : list of SRK or PR
        EoS for which binary interaction parameter will be optimized.
    predict : Predict
        Which property of the mixture to predict.
    exp : dict {"x", "y", "t", "p"}
        Experimental property data for the error calculation.
    ws : dict {"x", "y", "t", "p"}
        Weights for calculation of the error.
    Returns
    -------
    mse : float
        Weighted sum of all mean square errors.
    """
    for eq in eqs:
        eq.mix.k = [[0, k], [k, 0]]
    eps = [find_equilibrium(eq, eq.mix.T, eq.mix.P, predict) for eq in eqs]
    mse = 0
    if "t" in exp:
        ts = [(ep[0].T, t) for ep, t in zip(eps, exp["t"]) if ep[1]]
        res = sum([(p - e)**2 for p, e in ts]) / len(ts) * ws["t"]
        mse += sum(exp["t"]) / len(exp["t"]) * \
            ws["t"] if np.isnan(res) else res
    if "p" in exp:
        ps = [(ep[0].P, p) for ep, p in zip(eps, exp["p"]) if ep[1]]
        res = sum([(p - e)**2 for p, e in ps]) / len(ps) * ws["p"]
        mse += sum(exp["p"]) / len(exp["p"]) * \
            ws["p"] if np.isnan(res) else res
    if "x" in exp:
        xs = [(ep[0].x[0][0], x) for ep, x in zip(eps, exp["x"]) if ep[1]]
        res = sum([(p - e)**2 for p, e in xs]) / len(xs) * ws["x"]
        mse += sum(exp["x"]) / len(exp["x"]) * \
            ws["x"] if np.isnan(res) else res
    if "y" in exp:
        ys = [(ep[0].y[0][0], y) for ep, y in zip(eps, exp["y"]) if ep[1]]
        res = sum([(p - e)**2 for p, e in ys]) / len(ys) * ws["y"]
        mse += sum(exp["y"]) / len(exp["y"]) * \
            ws["y"] if np.isnan(res) else res

    return mse


def fit(eq, mix, const, exp, ws=None):
    """
    Obtain binary interaction parameter for the given EOS
    and experimental data.
    Parameters
    ----------
    eq : SRK, PR
        Equation of state to optimize.
    mix : Mixture
        Investigating mixture.
    const : float
        Value of the constant parameter.
    exp : dict
        Experimental values (without constant parameter) of the:
        p, t, x, y
    ws : dict, optional
        Weigths of the coefficients for the optimization algorithm.
        p, t, x, y
        Default is {"x": 1, "y": 1, "t": 1e-6, "p": 1e-12}.
    Returns
    -------
    k : float
        Optimized binary interaction parameter.
    mse : float
        Mean square error of the optimized prediction.
    r2 : float
        R squared of the optimized prediction.
    """
    if ws is None:
        ws = {
            "x": 1,
            "y": 1,
            "t": 1e-6,
            "p": 1e-12
        }

    if "x" in exp:
        mixs = [mix.set_x([x, 1-x]).set_y([x, 1-x]) for x in exp["x"]]
        for mix in mixs:
            mix.phase = Phase.LIQUID
    else:
        mixs = [mix.set_x([y, 1-y]).set_y([y, 1-y]) for y in exp["y"]]
        for mix in mixs:
            mix.phase = Phase.VAPOR

    if "t" in exp:
        for mix in mixs:
            mix.T = exp["t"][0]*0.8
            mix.P = const
        predict = Predict.TEMPERATURE
    else:
        for mix in mixs:
            mix.P = exp["p"][0]*0.2
            mix.T = const
        predict = Predict.PRESSURE

    eqs = [eq(mix) for mix in mixs]

    k, r, it, nf = optimize.fminbound(
        optim_function, -0.2, 0.2,
        args=(eqs, predict, exp, ws),
        full_output=True
    )

    print(f"r: {r}\nit: {it}\nnf: {nf}")

    return k
