"""Real-gas volume ratio; explicit standard conditions (60°F, 1 atm, Zsc=1)."""

STANDARD_PRESSURE = 101325.0
STANDARD_TEMPERATURE = (60 - 32) * 5 / 9 + 273.15


def real_gas_fvf(pressure: float, temperature: float, z: float) -> float:
    """Bg = Z*T*Psc/(P*Tsc*Zsc), canonical reservoir m³/standard m³."""
    return z * temperature * STANDARD_PRESSURE / (pressure * STANDARD_TEMPERATURE)
