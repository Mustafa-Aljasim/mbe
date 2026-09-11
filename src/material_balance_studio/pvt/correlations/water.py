"""McCain water FVF as documented by Pengtools, Properties of Petroleum Fluids.

Inputs psia and °F; result rb/STB (=m³/m³). Salinity is not an explicit term in
this approximation; the input's value is retained and its treatment disclosed.
"""


def mccain_bw(p: float, t: float) -> float:
    dvp = -1.95301e-9 * p * t - 1.72834e-13 * p * p * t - 3.58922e-7 * p - 2.25341e-10 * p * p
    dvt = -1.0001e-2 + 1.33391e-4 * t + 5.50654e-7 * t * t
    return (1 + dvp) * (1 + dvt)
