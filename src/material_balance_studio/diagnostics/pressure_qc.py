"""Row-preserving raw-history screening and optional pressure metadata."""
from dataclasses import dataclass
from math import isfinite
import pandas as pd
from material_balance_studio.io.history import history_from_frame
from material_balance_studio.units.conversions import to_si

SOURCES = ("Average reservoir pressure", "Static well pressure", "PBU-derived pressure",
           "RFT/MDT pressure", "Estimated pressure", "Other")
METADATA = ("pressure_source", "pressure_sigma", "pressure_quality", "pressure_note")


@dataclass(frozen=True)
class PressureObservation:
    date: str
    source: str = ""
    sigma: float | None = None  # Pa, standard deviation; never an optimization weight
    quality: str = ""
    note: str = ""


def number(value):
    try:
        result = float(value)
        return result if isfinite(result) else None
    except (TypeError, ValueError):
        return None


def text(value):
    return "" if pd.isna(value) else str(value)


def pressure_metadata(frame, units="SI"):
    def date_key(value):
        parsed=pd.to_datetime(value,errors="coerce")
        return str(value) if pd.isna(parsed) else parsed.date().isoformat()
    return tuple(PressureObservation(date_key(row.get("date", "")), text(row.get("pressure_source", "")),
                 to_si(number(row.get("pressure_sigma")), "pressure", units) if number(row.get("pressure_sigma")) is not None else None,
                 text(row.get("pressure_quality", "")), text(row.get("pressure_note", "")))
                 for row in frame.to_dict("records"))


def diagnostic_history_from_frame(frame, units="SI"):
    """Strip only documented metadata before the unchanged strict engine parser."""
    return history_from_frame(frame.drop(columns=list(METADATA), errors="ignore"), units), pressure_metadata(frame, units)


def history_qc(frame, initial_pressure, initial_date, units="SI"):
    """No repairs or deletion. Thresholds: 1% depletion, 365-day gaps, 20% Pi jumps."""
    flags = []
    def flag(row, status, rule, note):
        flags.append(dict(Row=row, Status=status, Rule=rule, Note=note))
    if number(initial_pressure) is None or initial_pressure<=0:
        flag(None,"FAIL","Initial pressure","Initial absolute pressure must be positive before screening depletion.")
        return pd.DataFrame(flags)
    if frame.empty or frame.columns.duplicated().any():
        flag(None,"FAIL","History schema","Empty history or duplicate column names.")
        return pd.DataFrame(flags)
    dates = pd.to_datetime(frame.get("date", pd.Series(index=frame.index,dtype=object)), errors="coerce")
    previous_date = pd.Timestamp(initial_date)
    previous = {k:0. for k in ("np","gp","wp","winj","ginj")}
    last_survey, last_pressure, pressures = pd.Timestamp(initial_date), initial_pressure, []
    for i, row in enumerate(frame.to_dict("records")):
        rowno, when = i+2, dates.iloc[i]
        if pd.isna(when):
            flag(rowno,"FAIL","Date","Missing or invalid date.")
        else:
            if when < previous_date:
                flag(rowno,"FAIL","Chronology","Date precedes previous row or initial date.")
            if dates.duplicated(keep=False).iloc[i]:
                flag(rowno,"FAIL","Duplicate date","Duplicate history/pressure date; review both rows.")
        for k in previous:
            value = number(row.get(k,0.))
            if value is None or value < 0:
                flag(rowno,"FAIL",k,"Missing, nonfinite or negative cumulative volume.")
                continue
            if value < previous[k]:
                flag(rowno,"FAIL",k,"Cumulative production/injection decreases.")
            if not pd.isna(when) and when <= previous_date and value > previous[k]:
                flag(rowno,"FAIL",k,"Production/injection increases without positive elapsed time.")
            previous[k] = value
        raw = row.get("observed_pressure")
        p = number(raw)
        if p is None:
            flag(rowno,"CAUTION" if pd.isna(raw) else "FAIL","Pressure","Missing pressure observation." if pd.isna(raw) else "Invalid pressure observation.")
        else:
            p = to_si(p,"pressure",units)
            if p <= 0:
                flag(rowno,"FAIL","Pressure","Observed absolute pressure must be positive.")
            else:
                pressures.append(p)
                if p > initial_pressure:
                    flag(rowno,"CAUTION","Pressure above Pi","May reflect injection or survey inconsistency; retain and review.")
                if abs(p-last_pressure) > .2*initial_pressure:
                    flag(rowno,"CAUTION","Abrupt change","Survey change exceeds 20% of initial pressure.")
                if not pd.isna(when) and (when-last_survey).days > 365:
                    flag(rowno,"CAUTION","Survey gap","More than 365 days since previous survey/reference.")
                if not pd.isna(when):
                    last_survey = when
                last_pressure = p
        sigma_raw = row.get("pressure_sigma")
        if sigma_raw is not None and not pd.isna(sigma_raw):
            sigma = number(sigma_raw)
            if sigma is None or sigma <= 0:
                flag(rowno,"CAUTION","Pressure uncertainty","Sigma must be positive; normalized residual unavailable.")
        source = text(row.get("pressure_source", ""))
        if source and source not in SOURCES:
            flag(rowno,"CAUTION","Pressure source","Unrecognized source label retained.")
        if source and source != SOURCES[0]:
            flag(rowno,"CAUTION","Pressure representativeness","Verify that this measurement represents average tank pressure.")
        quality = text(row.get("pressure_quality", ""))
        if quality and quality.upper() not in ("PASS","GOOD"):
            flag(rowno,"CAUTION","Engineering quality",quality)
        if not pd.isna(when):
            previous_date = when
    if len(pressures) < 5:
        flag(None,"CAUTION","Survey count",f"Only {len(pressures)} valid pressure observations.")
    if pressures and (initial_pressure-min(pressures))/initial_pressure < .01:
        flag(None,"CAUTION","Small depletion","Pressure depletion below 1% of initial pressure; transformed diagnostics may be ill-conditioned.")
    if not frame.empty and (previous_date-last_survey).days > 365:
        flag(None,"CAUTION","Survey gap","Final history extends more than 365 days beyond last pressure survey.")
    if not flags:
        flag(None,"PASS","History screening","No configured screening flags.")
    return pd.DataFrame(flags, columns=["Row","Status","Rule","Note"])
