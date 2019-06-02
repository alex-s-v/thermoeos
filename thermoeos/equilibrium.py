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
    ep, calculated = f(init)
    if not calculated:
        raise Exception("Initial parameters must be solvable.")
    s0 = ep.s
    for _ in range(max_steps):
        init += step_size
        ep, calculated = f(init)
        if calculated:
            if np.sign(s0-1) != np.sign(ep.s-1):
                return sorted([init-step_size, init])
            else:
                if np.abs(ep.s-1) > np.abs(s0-1):
                    step_size *= -1
                else:
                    s0 = ep.s
        else:
            raise Exception(f"Can't solve equation for {init}")
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
        step = 1e3
        max_steps = 1000
        init = P

        def solve(val): return eq.solve(T, val)
    elif predict == Predict.TEMPERATURE:
        step = 2
        max_steps = 100
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


def calc_rmse(calc, exp):
    diff = np.subtract(calc, exp)
    sum_ = np.sum(np.abs(diff)**2)
    return np.sqrt(np.mean(sum_))


def get_prop(mixs, prop):
    if prop == "t":
        return [mix.T for mix in mixs]
    elif prop == "p":
        return [mix.P for mix in mixs]
    elif prop == "x":
        return [mix.x[0][0] for mix in mixs]
    elif prop == "y":
        return [mix.y[0][0] for mix in mixs]
    else:
        raise Exception(f"Unknown poperty: {prop}")


def solve_mult(k, eos, predict):
    for eq in eos:
        eq.mix.set_k([[0, k], [k, 0]], inplace=True)
    eps = [find_equilibrium(eq, eq.mix.T, eq.mix.P, predict) for eq in eos]
    return eps


def optimization_function(k, eos, ems, predict, include, ws):
    """
    Calculate error
    Parameters
    ----------
    eos : list of thermoeos.equations.Equation
        EoS for the optimization.
    ems : list of thermoeos.mixture.Mixture
        Mixture whose parameters are obtained experimentally.
    include : list of str, optional
        Errors of which parameters include in the computation
        of BIP. Any combination of t, p, x and y. All included
        by default.
    Returns
    -------
    rmse : float
        Root of mean square error for the given set of
        experimental mixtures and calculated ones
        using obtained BIP.
    """
    eps = solve_mult(k, eos, predict)
    mixs = [ep[0] for ep in eps if ep[1]]
    if len(mixs) != len(eps):
        raise Exception(f"Not all EoS are solved. ({len(mixs)}/{len(eps)})")
    rmse = 0
    for prop in include:
        rmse += calc_rmse(get_prop(mixs, prop), get_prop(ems, prop)) * ws[prop]
    return rmse


def new_fit(eos, ems, include=["t", "p", "x", "y"], bounds=[-0.2, 0.2],
            ws={"x": 1, "y": 1, "t": 1e-2, "p": 1e-5},
            predict=Predict.PRESSURE):
    """
    Calculate one BIP for the given set of mixtures.
    Parameters
    ----------
    eos : list of thermoeos.equations.Equation
        EoS for the optimization.
    ems : list of thermoeos.mixture.Mixture
        Mixture whose parameters are obtained experimentally.
    include : list of str, optional
        Errors of which parameters include in the computation
        of BIP. Any combination of t, p, x and y. All included
        by default.
    Returns
    -------
    k : float
        Optimized BIP for the given set on mixtures.
    rmse : dict
        Root of mean square error for the given set of
        experimental mixtures and calculated ones
        using obtained BIP.
        rmse = {
            "x": float or None,
            "y": float or None,
            "t": float or None,
            "p": float or None
        }
    """
    k = optimize.fminbound(
        optimization_function, *bounds,
        args=(eos, ems, predict, include, ws)
    )

    eps = solve_mult(k, eos, predict)
    mixs = [ep[0] for ep in eps if ep[1]]

    rmse = {"x": None, "y": None, "t": None, "p": None}

    for prop in include:
        rmse[prop] = calc_rmse(get_prop(mixs, prop), get_prop(ems, prop))

    return k, rmse


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
