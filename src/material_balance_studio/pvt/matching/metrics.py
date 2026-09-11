from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class Metrics:
    n: int
    rmse: float
    mae: float
    mape: float | None
    max_ape: float | None
    bias: float
    percentage_n: int


def error_metrics(measured, predicted) -> Metrics:
    y, prediction = np.asarray(measured, dtype=float), np.asarray(predicted, dtype=float)
    if y.ndim != 1 or y.shape != prediction.shape or not y.size or not np.all(np.isfinite(y)) or not np.all(np.isfinite(prediction)):
        raise ValueError("Metrics require equal nonempty finite one-dimensional arrays.")
    error = prediction-y
    nonzero = y != 0
    ape = 100*np.abs(error[nonzero]/y[nonzero])
    return Metrics(len(y), float(np.sqrt(np.mean(error**2))), float(np.mean(np.abs(error))),
                   float(np.mean(ape)) if len(ape) else None, float(np.max(ape)) if len(ape) else None,
                   float(np.mean(error)), int(np.sum(nonzero)))
