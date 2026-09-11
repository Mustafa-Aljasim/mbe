"""Run the bundled synthetic case without importing Streamlit."""

from datetime import date
from pathlib import Path

from material_balance_studio.diagnostics.balance_closure import results_frame
from material_balance_studio.domain.models import ReservoirTank
from material_balance_studio.io.history import history_from_frame, pvt_from_frame, read_table
from material_balance_studio.pvt.table_model import TablePVTModel
from material_balance_studio.solver.simulation import simulate


def main() -> None:
    folder = Path(__file__).resolve().parent
    pvt = TablePVTModel(pvt_from_frame(read_table(folder / "pvt.csv")))
    tank = ReservoirTank(date(2020, 1, 1), 30e6, 1e6, 0.2, 5e-10, 4e-10, 0.0, pvt)
    result = simulate(tank, history_from_frame(read_table(folder / "history.csv")))
    print(results_frame(result)[["date", "pressure_pa", "residual_m3", "relative_residual"]].to_string(index=False))
    if not result.converged:
        raise SystemExit(result.failed_timestep.solution.diagnostics.message)


if __name__ == "__main__":
    main()
