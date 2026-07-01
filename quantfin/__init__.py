import math

def _phi(x):
    return 0.5 * (1 + math.erf(x / math.sqrt(2)))

def _pdf(x):
    return math.exp(-0.5 * x * x) / math.sqrt(2 * math.pi)

def black_scholes(S, K, r, sigma, T, call=True):
    d1 = (math.log(S / K) + (r + 0.5 * sigma ** 2) * T) / (sigma * math.sqrt(T))
    d2 = d1 - sigma * math.sqrt(T)
    if call:
        return S * _phi(d1) - K * math.exp(-r * T) * _phi(d2)
    return K * math.exp(-r * T) * _phi(-d2) - S * _phi(-d1)

def greeks(S, K, r, sigma, T, call=True):
    d1 = (math.log(S / K) + (r + 0.5 * sigma ** 2) * T) / (sigma * math.sqrt(T))
    d2 = d1 - sigma * math.sqrt(T)
    delta = _phi(d1) if call else _phi(d1) - 1
    gamma = _pdf(d1) / (S * sigma * math.sqrt(T))
    vega = S * _pdf(d1) * math.sqrt(T)
    return {"delta": delta, "gamma": gamma, "vega": vega}

def npv(rate, cashflows):
    return sum(cf / (1 + rate) ** t for t, cf in enumerate(cashflows))

def irr(cashflows, guess=0.1, tol=1e-9, max_iter=200):
    r = guess
    for _ in range(max_iter):
        f = npv(r, cashflows)
        df = sum(-t * cf / (1 + r) ** (t + 1) for t, cf in enumerate(cashflows))
        if df == 0:
            break
        step = f / df
        r -= step
        if abs(step) < tol:
            break
    return r
