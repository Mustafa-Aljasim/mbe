"""PETEX Appendix C C2.9 piecewise-linear pressure convolution; infinite radial."""
from dataclasses import dataclass, replace
from math import expm1, log, pi, sqrt
from typing import ClassVar

from scipy.integrate import quad

from .veh import VanEverdingenHurstAquifer
from .veh_response import VEH_INFINITE_TABLE, veh_infinite_water_influx_dimensionless as response
from .transient import transient_warnings
from material_balance_studio.domain.validation import finite_value


def integrated_response(low, high):
    """Integral over response ages, split at all shared interpolation boundaries."""
    finite_value("lower response age", low, nonnegative=True)
    finite_value("upper response age", high, nonnegative=True)
    if high < low:
        raise ValueError("Response ages must be ordered.")
    cuts = [low, *(t for t, _ in VEH_INFINITE_TABLE if low < t < high), high]
    total = 0.0
    for left, right in zip(cuts, cuts[1:]):
        if right == left:
            continue
        if right <= .01:
            total += 4/(3*sqrt(pi))*(right**1.5-left**1.5)
        elif left >= VEH_INFINITE_TABLE[-1][0]:
            value, error = quad(response, left, right, epsabs=1e-10, epsrel=1e-10)
            if error > max(1e-8, abs(value)*1e-8):
                raise ValueError("Modified VEH response integration did not converge.")
            total += value
        else:
            # The entire subinterval belongs to one log-log interpolation cell.
            midpoint = (left+right)/2
            index = next(i for i, (t, _) in enumerate(VEH_INFINITE_TABLE) if t >= midpoint)
            t0, w0 = VEH_INFINITE_TABLE[index-1]
            t1, w1 = VEH_INFINITE_TABLE[index]
            power = log(w1/w0)/log(t1/t0)
            wl = w0*(left/t0)**power
            total += wl*left*expm1((power+1)*log(right/left))/(power+1)
    return total


@dataclass(frozen=True)
class ModifiedVanEverdingenHurstAquifer(VanEverdingenHurstAquifer):
    key: ClassVar[str] = "modified_van_everdingen_hurst"

    def compute_step(self, previous, previous_reservoir_pressure, candidate_reservoir_pressure, dt):
        average = self.validate_step(previous, previous_reservoir_pressure, candidate_reservoir_pressure, dt)
        values = dict(previous.model_variables)
        count = int(values.get("segment_count", 0))
        segments = [(values[f"segment_start_{i}"], values[f"segment_end_{i}"],
                     values[f"segment_drop_{i}"]) for i in range(count)]
        elapsed = previous.elapsed_time+dt
        drop = previous_reservoir_pressure-candidate_reservoir_pressure
        segments.append((previous.elapsed_time, elapsed, drop))
        geometry = self.geometry
        contributions, warnings = [], []
        for start, end, dp in segments:
            low, high = geometry.dimensionless_time(elapsed-end), geometry.dimensionless_time(elapsed-start)
            integral = integrated_response(low, high)
            contributions.append(geometry.aquifer_constant*dp/geometry.dimensionless_time(end-start)*integral)
            warnings.extend(transient_warnings(tD=high, response_name="Modified VEH infinite response"))
        variables = [("segment_count", float(len(segments)))]
        for i, (start, end, dp) in enumerate(segments):
            variables.extend(((f"segment_start_{i}", start), (f"segment_end_{i}", end), (f"segment_drop_{i}", dp)))
        variables.extend(("diag_"+name, value) for name, value in {
            "tD": geometry.dimensionless_time(elapsed), "start_pressure": previous_reservoir_pressure,
            "end_pressure": candidate_reservoir_pressure, "average_pressure": average,
            "pressure_slope": -drop/dt, "response_value": response(geometry.dimensionless_time(elapsed)),
            "current_step_contribution": contributions[-1], "historical_contribution": sum(contributions[:-1]),
            "integrated_current_response": integrated_response(0, geometry.dimensionless_time(dt)),
        }.items())
        result = self.step_result(previous, candidate_reservoir_pressure, dt, sum(contributions), None,
                                  average, model_variables=tuple(variables))
        return replace(result, warnings=tuple(dict.fromkeys((*result.warnings, *warnings))))
