from dataclasses import replace
from datetime import date
from io import BytesIO, StringIO

import pandas as pd
import pytest

from material_balance_studio.domain.models import CUMULATIVE_FIELDS, CumulativeVolumes, HistoryRecord
from material_balance_studio.domain.validation import EngineeringValidationError
from material_balance_studio.io.history import history_from_frame, pvt_from_frame, read_table, validate_history
from material_balance_studio.solver.simulation import simulate


@pytest.mark.parametrize("field", CUMULATIVE_FIELDS)
def test_e_decreasing_cumulatives_rejected(field):
    history = [HistoryRecord(date(2020, month, 1), CumulativeVolumes(**{field: value}))
               for month, value in [(2, 100), (3, 90)]]
    with pytest.raises(EngineeringValidationError, match=f"Cumulative {field} decreases"):
        validate_history(history)


@pytest.mark.parametrize("field", CUMULATIVE_FIELDS)
@pytest.mark.parametrize("value", [-1.0, float("nan"), float("inf")])
def test_negative_and_nonfinite_volumes_rejected(field, value):
    with pytest.raises(EngineeringValidationError):
        CumulativeVolumes(**{field: value})


def test_unsorted_duplicate_empty_and_preinitial_history(tank):
    for history in [[], [HistoryRecord(date(2020, 3, 1)), HistoryRecord(date(2020, 2, 1))],
                    [HistoryRecord(date(2020, 2, 1))] * 2,
                    [HistoryRecord(date(2019, 1, 1))],
                    [HistoryRecord(tank.initial_date, CumulativeVolumes(np=1))]]:
        with pytest.raises(EngineeringValidationError):
            simulate(tank, history)


@pytest.mark.parametrize("field,value", [("swc", -0.1), ("swc", 1), ("m", -1),
                                         ("oil_in_place", 0), ("cf", -1), ("cw", -1),
                                         ("initial_pressure", 5e6)])
def test_reservoir_safeguards(tank, field, value):
    with pytest.raises(EngineeringValidationError):
        replace(tank, **{field: value})


def test_csv_excel_and_missing_observations():
    csv = StringIO("date,np,gp,wp,observed_pressure\n2020-02-01,100,8000,0,\n")
    csv_frame = read_table(csv, "history.csv")
    history = history_from_frame(csv_frame)
    assert history[0].observed_pressure is None
    assert history[0].cumulative.winj == 0
    workbook = BytesIO()
    csv_frame.to_excel(workbook, index=False)
    workbook.seek(0)
    assert history_from_frame(read_table(workbook, "history.xlsx")) == history


def test_pvt_excel_roundtrip(pvt):
    frame = pd.DataFrame([{"pressure": r.pressure, "bo": r.bo, "rs": r.rs, "bg": r.bg, "bw": r.bw}
                          for r in pvt.table.rows])
    workbook = BytesIO()
    frame.to_excel(workbook, index=False)
    workbook.seek(0)
    assert pvt_from_frame(read_table(workbook, "pvt.xlsx")) == pvt.table


@pytest.mark.parametrize("patch", [{"date": "02/01/2020"}, {"np": "bad"}, {"np": None},
                                    {"winj": None}, {"observed_pressure": -1}])
def test_malformed_rows_rejected(patch):
    row = {"date": "2020-02-01", "np": 1, "gp": 80, "wp": 0, **patch}
    with pytest.raises(EngineeringValidationError):
        history_from_frame(pd.DataFrame([row]))


def test_wrong_columns_are_not_silently_ignored():
    with pytest.raises(EngineeringValidationError, match="Invalid columns"):
        history_from_frame(pd.DataFrame([{"Date": "2020-02-01", "np": 1, "gp": 80, "wp": 0}]))
