"""Python port of the fitting logic in GrowthPredict's fit_model.m (repo HEAD b7829ae).

Ported: the seven growth models, the incidence curve abs([C(1); diff(C)]), the
least-squares and Poisson objectives, the parameter box, the seed, the start
points (Latin hypercube + jittered seed + seed), the unit-box / log10 optimizer
coordinates and the "widen r once" retry.

Not ported: fmincon/MultiStart (scipy least_squares is used instead; Poisson is
minimised through deviance residuals, which has the same minimiser as the NLL).

Box schemes
  'old'  the box of fit_model.m at HEAD
  'new'  proposed: for the models with a free exponent p (GGM, GLM, GRM)
         r in [1e-4, 10*max(y)] / DT; all other models unchanged
  'i0'   optimizer searches the initial incidence rate i0 = r*C0^p*sat(C0)
         instead of r (the "Fix 2" design), box [yref/20, 20*yref]

Gompertz (flag 5) only
  'old'      HEAD: time labels used as they are in r*C*exp(-a*t); box from initialParams
  'gfix'     adopted: time measured from the first point of the curve; HEAD box divided by DT
  'gnew'     tested and not adopted: as 'gfix' with r and a in [1e-4, 10]/DT, both searched in log10
  'gwide'    reference: as 'gnew' with [1e-7, 1e3]/DT
"""
import math
import numpy as np
from numba import njit
from scipy.optimize import least_squares
from scipy.stats import qmc

EXP, GGM, GLM, GRM, LM, RICH, GOM = -1, 0, 1, 2, 3, 4, 5
NAMES = {-1: 'EXP', 0: 'GGM', 1: 'GLM', 2: 'GRM', 3: 'LM', 4: 'RICH', 5: 'GOM'}
FREE_P = (GGM, GLM, GRM)
ELAPSED = ('gfix', 'gnew', 'gwide')      # Gompertz schemes that measure time from the first point


@njit(cache=True)
def _f(flag, t, C, r, p, a, K):
    if C < 0.0:
        C = 0.0
    if flag == -1:
        return r * C
    if flag == 0:
        return r * C ** p
    if flag == 1:
        return r * C ** p * (1.0 - C / K)
    if flag == 2:
        return r * C ** p * (1.0 - (C / K) ** a)
    if flag == 3:
        return r * C * (1.0 - C / K)
    if flag == 4:
        return r * C * (1.0 - (C / K) ** a)
    return r * C * math.exp(-a * t)


@njit(cache=True)
def solve_cum(flag, r, p, a, K, C0, tvec, rtol, atol, maxsteps):
    """Adaptive Dormand-Prince 5(4) for the scalar model; returns (C, ok)."""
    n = tvec.shape[0]
    out = np.empty(n)
    out[0] = C0
    t = tvec[0]
    C = C0
    sat = (flag == 1) or (flag == 2) or (flag == 3) or (flag == 4)
    h = (tvec[1] - tvec[0]) * 1e-3
    steps = 0
    for i in range(1, n):
        tend = tvec[i]
        while t < tend:
            if sat and abs(K - C) <= 1e-13 * K:      # at the plateau: nothing changes any more
                C = K
                t = tend
                break
            if h > tend - t:
                h = tend - t
            k1 = _f(flag, t, C, r, p, a, K)
            k2 = _f(flag, t + h / 5, C + h * k1 / 5, r, p, a, K)
            k3 = _f(flag, t + 3 * h / 10, C + h * (3 * k1 / 40 + 9 * k2 / 40), r, p, a, K)
            k4 = _f(flag, t + 4 * h / 5, C + h * (44 * k1 / 45 - 56 * k2 / 15 + 32 * k3 / 9), r, p, a, K)
            k5 = _f(flag, t + 8 * h / 9, C + h * (19372 * k1 / 6561 - 25360 * k2 / 2187 + 64448 * k3 / 6561 - 212 * k4 / 729), r, p, a, K)
            k6 = _f(flag, t + h, C + h * (9017 * k1 / 3168 - 355 * k2 / 33 + 46732 * k3 / 5247 + 49 * k4 / 176 - 5103 * k5 / 18656), r, p, a, K)
            Cn = C + h * (35 * k1 / 384 + 500 * k3 / 1113 + 125 * k4 / 192 - 2187 * k5 / 6784 + 11 * k6 / 84)
            k7 = _f(flag, t + h, Cn, r, p, a, K)
            err = h * abs(71 * k1 / 57600 - 71 * k3 / 16695 + 71 * k4 / 1920 - 17253 * k5 / 339200 + 22 * k6 / 525 - k7 / 40)
            if not (math.isfinite(Cn) and math.isfinite(err)) or abs(Cn) > 1e200:
                if h < 1e-14 * max(1.0, abs(t)):
                    return out, False
                h *= 0.1
                steps += 1
                if steps > maxsteps:
                    return out, False
                continue
            tol = atol + rtol * max(abs(C), abs(Cn))
            if err <= tol:
                t += h
                C = Cn
                if C < 0.0:
                    C = 0.0
            fac = 0.9 * (tol / (err + 1e-300)) ** 0.2
            if fac > 5.0:
                fac = 5.0
            if fac < 0.2:
                fac = 0.2
            h *= fac
            steps += 1
            if steps > maxsteps or h < 1e-15 * max(1.0, abs(t)):
                return out, False
        out[i] = C
    return out, True


def incidence(z, tvec, flag, rtol=1e-10, atol=1e-10):
    """growthIncidenceCurve.m: abs([C(1); diff(C)]) or None when the solve fails."""
    C, ok = solve_cum(flag, z[0], z[1], z[2], z[3], z[4], np.ascontiguousarray(tvec, dtype=float), rtol, atol, 400000)
    if not ok or not np.all(np.isfinite(C)):
        return None
    return np.abs(np.r_[C[0], np.diff(C)])


def initial_params(y, flag):
    """initialParams.m"""
    s, y1 = float(np.sum(y)), float(y[0])
    if flag == -1:
        return [1, 1, 1, s, y1, 1, 1]
    if flag == 0:
        return [1, 0.9, 1, s, y1, 1, 1]
    if flag in (1, 2):
        return [0.1, 0.9, 1, s, y1, 1, 1]
    if flag == 3:
        return [0.1, 1, 1, s, y1, 1, 1]
    if flag == 4:
        return [0.1, 1, 0.9, s, y1, 1, 1]
    rg = 1 - y1 / s
    ag = rg / math.log(s / y1)
    return [rg + 0.001, 1, ag, s, y1, 1, 1]


def build_box(t, y, flag, fixI0, scheme, params0=None, DT=1.0):
    """Bounds and seed exactly as in fit_model.m ('old'); 'new'/'i0' change the first slot only.
    Vector layout [r p a K I0]; the error-model slots are left out (LS and Poisson have none)."""
    if params0 is None:
        params0 = initial_params(y, flag)
    I0 = float(y[0])
    w = min(7, len(y))
    r0, rlb, rub = params0[0], 1e-3, 10.0
    slope = np.nan
    if w >= 2:
        slope = np.polyfit(t[:w] * DT, np.log(np.maximum(y[:w], 0.5)), 1)[0]
        if slope > 0:
            r0, rlb, rub = slope, slope / 20, slope * 20
    Cobs = float(np.sum(np.abs(y)))
    Klb = max(20.0, Cobs)
    Kub = max(min(1e3 * Cobs, 1e9), Klb)

    if scheme in ('new', 'i0', 'newa') and flag in FREE_P:
        ymax = max(float(np.max(y)), 1.0)
        rlb, rub = 1e-4 / DT, 10.0 * ymax / DT          # <-- the proposed box
        if scheme == 'newa' and flag == GRM:            # variant: a >= 0.05 and r scaled by 1/a_min
            rub = rub / 0.05

    if fixI0 == 1:
        ilb = iub = I0
    else:
        iub = min(5 * max(1.0, I0), Klb)
        ilb = min(max(1.0, 0.2 * I0), iub)
    p_lo = {-1: 1, 0: .01, 1: .01, 2: .01, 3: 1, 4: 1, 5: 1}[flag]
    a_lo, a_hi = {-1: (1, 1), 0: (1, 1), 1: (1, 1), 2: (0, 10), 3: (1, 1), 4: (0, 10)}.get(flag, (None, None))
    if scheme == 'newa' and flag == GRM:
        a_lo = 0.05
    if flag == 5 and scheme in ('gnew', 'gwide'):
        lo, hi = (1e-4, 10.0) if scheme == 'gnew' else (1e-7, 1e3)
        LB = [lo / DT, 1, lo / DT, 1, ilb]
        UB = [hi / DT, 1, hi / DT, 1, iub]
    elif flag == 5:
        s_ = DT if scheme == 'gfix' else 1.0          # rates are per unit of time
        LB = [max(1e-4, params0[0] / 10) / s_, 1, max(1e-6, params0[2] / 10) / s_, 1, ilb]
        UB = [params0[0] * 10 / s_, 1, params0[2] * 10 / s_, 1, iub]
    else:
        K_lo, K_hi = (0, 0) if flag in (-1, 0) else (Klb, Kub)
        LB = [rlb, p_lo, a_lo, K_lo, ilb]
        UB = [rub, 1, a_hi, K_hi, iub]
    z = list(params0[:5])
    z[4] = I0
    if flag != 5:
        z[0] = min(max(r0, LB[0]), UB[0])
    LB, UB = np.array(LB, float), np.array(UB, float)
    z = np.minimum(np.maximum(np.array(z, float), LB), UB)
    info = dict(slope=slope, Cobs=Cobs)

    if scheme == 'i0' and flag in FREE_P:               # slot 1 becomes the initial incidence rate
        yref = max(float(np.mean(y[:3])), 0.5)
        LB[0], UB[0] = yref / 20 / DT, yref * 20 / DT
        z[0] = min(max(r_to_i0(z, flag), LB[0]), UB[0])
    return LB, UB, z, info


def _sat0(z, flag):
    C0, K, a = z[4], z[3], z[2]
    if flag in (1, 3):
        return max(1 - C0 / K, 1e-12)
    if flag in (2, 4):
        return max(1 - (C0 / K) ** a, 1e-12)
    return 1.0


def r_to_i0(z, flag):
    return z[0] * max(z[4], 1e-12) ** z[1] * _sat0(z, flag)


def i0_to_r(z, flag):
    return z[0] / (max(z[4], 1e-12) ** z[1] * _sat0(z, flag))


def objective_value(mu, y, method):
    """growthObjectiveValue.m for method 0 (SSE) and 1 (Poisson NLL without log y!)."""
    mu = np.where(mu == 0, 0.001, mu)
    if method == 0:
        return float(np.sum((y - mu) ** 2))
    return float(-np.sum(y * np.log(mu) - mu))


_BIG = 1e9


def _squash(res):
    """Monotone compression of huge residuals (role of localBound in fit_model.m)."""
    a = np.abs(res)
    return np.where(a > _BIG, np.sign(res) * (_BIG + np.log1p(np.maximum(a - _BIG, 0))), res)


def _residuals(mu, y, method):
    if mu is None:
        return np.full(len(y), 2 * _BIG)
    mu = np.where(mu == 0, 0.001, mu)
    if method == 0:
        return _squash(mu - y)
    with np.errstate(divide='ignore', invalid='ignore'):
        dev = 2 * (np.where(y > 0, y * np.log(y / mu), 0.0) - (y - mu))
    return _squash(np.sign(mu - y) * np.sqrt(np.maximum(dev, 0)))


class Fit:
    __slots__ = ('P', 'fval', 'curve', 'LB', 'UB', 'atLB', 'atUB', 'n_starts', 'n_valid', 'agree', 'widened', 'scheme', 'slot1')


def fit_model(t, y, flag, fixI0=1, method=0, nstart=10, scheme='old', rng=None, params0=None,
              box=None, DT=1.0, _allow_widen=True):
    """One fit, following fit_model.m. Returns a Fit with P = [r p a K I0] in original units."""
    rng = np.random.default_rng(0) if rng is None else rng
    t = np.asarray(t, float)
    y = np.asarray(y, float)
    tv = t * DT
    if flag == GOM and scheme in ELAPSED:
        tv = tv - tv[0]
    if box is None:
        LB, UB, z0, _ = build_box(t, y, flag, fixI0, scheme, params0, DT)
        supplied = False
    else:
        LB, UB = np.array(box[0], float), np.array(box[1], float)
        z0 = np.minimum(np.maximum(np.array(params0[:5], float), LB), UB)
        supplied = True
    use_i0 = scheme == 'i0' and flag in FREE_P
    free = np.where(UB > LB)[0]
    logc = np.zeros(5, bool)
    for j in ((0, 2, 3) if (flag == GOM and scheme in ('gnew', 'gwide')) else (0, 3)):
        logc[j] = UB[j] > LB[j] and LB[j] > 0
    toT = lambda Z: np.where(logc, np.log10(np.maximum(Z, 1e-300)), Z)
    fromT = lambda T: np.where(logc, 10.0 ** T, T)
    LBt, UBt = toT(LB), toT(UB)
    span = UBt - LBt

    def to_model(u):
        T = toT(z0).copy()
        T[free] = LBt[free] + np.clip(u, 0, 1) * span[free]
        Z = fromT(T)
        if use_i0:
            Z = Z.copy()
            Z[0] = i0_to_r(Z, flag)
        return Z

    def res(u):
        return _residuals(incidence(to_model(u), tv, flag), y, method)

    nF = len(free)
    u_seed = (toT(z0)[free] - LBt[free]) / span[free]
    starts = []
    if nstart > 0 and nF > 0:
        starts += list(qmc.LatinHypercube(d=nF, seed=rng).random(nstart))
        for _ in range(min(3, nstart)):
            starts.append(np.clip(u_seed + 0.05 * rng.standard_normal(nF), 0, 1))
    starts.append(np.clip(u_seed, 0, 1))

    sols = []
    for u0 in starts:
        try:
            s = least_squares(res, u0, bounds=(np.zeros(nF), np.ones(nF)), xtol=1e-12, ftol=1e-12,
                              gtol=1e-12, diff_step=1e-6, max_nfev=600)
        except Exception:
            continue
        Z = to_model(s.x)
        mu = incidence(Z, tv, flag)
        if mu is None:
            continue
        fv = objective_value(mu, y, method)
        if np.isfinite(fv):
            sols.append((fv, Z, mu, s.x))
    out = Fit()
    out.scheme, out.n_starts, out.n_valid, out.widened = scheme, len(starts), len(sols), False
    out.LB, out.UB = LB, UB
    if not sols:
        out.P, out.fval, out.curve, out.agree = None, np.inf, None, np.nan
        out.atLB = out.atUB = np.zeros(5, bool)
        out.slot1 = np.nan
        return out
    fvs = np.array([s[0] for s in sols])
    b = int(np.argmin(fvs))
    out.fval, out.P, out.curve = sols[b][0], sols[b][1], sols[b][2]
    out.agree = float(np.mean(fvs <= out.fval + 1e-4 * max(1.0, abs(out.fval))))
    # bound flags are evaluated on the optimizer's own slot-1 quantity (r, or i0 for scheme 'i0')
    Pq = out.P.copy()
    if use_i0:
        Pq[0] = r_to_i0(out.P, flag)
    out.slot1 = Pq[0]
    out.atLB, out.atUB = np.zeros(5, bool), np.zeros(5, bool)
    out.atLB[free] = Pq[free] <= LB[free] + 1e-3 * np.abs(LB[free])
    out.atUB[free] = Pq[free] >= UB[free] - 1e-3 * np.abs(UB[free])

    # fit_model.m: escape from a binding data-derived interval for r
    if (_allow_widen and not supplied and not use_i0 and flag != 5 and (out.atLB[0] or out.atUB[0])
            and (LB[0] > 1e-3 or UB[0] < 10)):
        LB2, UB2 = LB.copy(), UB.copy()
        LB2[0], UB2[0] = min(LB[0], 1e-3), max(UB[0], 10.0)
        f2 = fit_model(t, y, flag, fixI0, method, nstart, scheme, rng, params0=list(out.P) + [1, 1],
                       box=(LB2, UB2), DT=DT, _allow_widen=False)
        if f2.fval < out.fval:
            f2.widened = True
            return f2
    return out


def forecast(P, t, h, flag, DT=1.0, elapsed=False):
    tv = np.arange(t[0], t[-1] + h + 0.5) * DT
    if elapsed:
        tv = tv - tv[0]
    return incidence(P, tv, flag)
