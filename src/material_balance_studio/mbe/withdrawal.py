"""Production and injection voidage, each evaluated at candidate reservoir pressure."""

from material_balance_studio.domain.models import CumulativeVolumes, PVTProperties, WithdrawalTerms


def produced_withdrawal(volumes: CumulativeVolumes, pvt: PVTProperties) -> float:
    """Np Bo + (Gp - Np Rs) Bg + Wp Bw; valid even when Np = 0.

    Do not clip the solution-gas correction: it is part of component accounting.
    """
    return volumes.np * pvt.bo + (volumes.gp - volumes.np * pvt.rs) * pvt.bg + volumes.wp * pvt.bw


def water_injection_support(volumes: CumulativeVolumes, pvt: PVTProperties) -> float:
    """Winj Bwinj; use resident Bw if no separate injection FVF is tabulated."""
    return volumes.winj * (pvt.bw if pvt.bwinj is None else pvt.bwinj)


def gas_injection_support(volumes: CumulativeVolumes, pvt: PVTProperties) -> float:
    """Ginj Bginj; injected gas is represented only by reservoir-volume support."""
    return volumes.ginj * (pvt.bg if pvt.bginj is None else pvt.bginj)


def withdrawal_terms(volumes: CumulativeVolumes, pvt: PVTProperties) -> WithdrawalTerms:
    return WithdrawalTerms(produced_withdrawal(volumes, pvt),
                           water_injection_support(volumes, pvt),
                           gas_injection_support(volumes, pvt))
