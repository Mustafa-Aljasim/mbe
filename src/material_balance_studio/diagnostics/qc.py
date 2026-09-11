"""Presentation QC categories; numerical solver acceptance remains unchanged."""

from math import isfinite

PASS_RELATIVE_TOLERANCE = 1e-8
WARNING_RELATIVE_TOLERANCE = 1e-5


def balance_qc_status(relative_residual: float | None, *, converged: bool = True) -> str:
    """A failed solve never receives PASS, even if its last residual is small."""
    if not converged or relative_residual is None or not isfinite(relative_residual):
        return "FAIL"
    error = abs(relative_residual)
    if error <= PASS_RELATIVE_TOLERANCE:
        return "PASS"
    if error <= WARNING_RELATIVE_TOLERANCE:
        return "WARNING"
    return "FAIL"
