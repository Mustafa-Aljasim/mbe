"""Water/gas/total voidage replacement using a common produced-volume basis."""
from .voidage import voidage


def ratio(injected, produced):
    return injected/produced if produced > 1e-12 else None


def vrr(state):
    volumes = voidage(state)
    result = {}
    for basis, values in volumes.items():
        result[basis] = {k:ratio(values[k],values["produced"]) for k in ("water_injected","gas_injected","injected")}
    result["warnings"] = tuple(f"{basis}: VRR undefined because produced voidage is zero, negative or below 1e-12 m³."
                               for basis,values in volumes.items() if values["produced"] <= 1e-12)
    return result
