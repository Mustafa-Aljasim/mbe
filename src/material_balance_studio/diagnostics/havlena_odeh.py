"""Observed-pressure material-balance transformations via the existing balance API."""
import numpy as np
from material_balance_studio.aquifer.none import NoAquifer
from material_balance_studio.mbe.oil_balance import evaluate_balance


def observed_components(tank, history):
    """Replay aquifer at observations independently of the solved pressure history.

    After a missing observation, transient We is unavailable: no silent pressure
    interpolation or borrowing from the calculated trajectory is allowed.
    """
    model = tank.aquifer or NoAquifer()
    aq = model.initial_state(tank.initial_pressure)
    previous_date, previous_pressure, complete = tank.initial_date, tank.initial_pressure, True
    rows = []
    for record in history:
        p = record.observed_pressure
        row = dict(date=str(record.date), pressure=p, F=None, Eo=None, Eg=None, Efw=None, Et=None,
                   We=None, injection=None, adjusted=None, status="PASS", note="")
        if p is None:
            complete = False
            row.update(status="CAUTION",note="Missing pressure; no observed-pressure transform.")
            rows.append(row)
            continue
        try:
            we = 0. if model.key == "none" else None
            if complete:
                dt = (record.date-previous_date).days*86400.
                if dt>0:
                    step = model.compute_step(aq,previous_pressure,p,dt)
                    aq = step.updated_state
                    if step.warnings:
                        row.update(status="CAUTION",note=" | ".join(step.warnings))
                we = aq.cumulative_influx
                previous_date, previous_pressure = record.date,p
            b = evaluate_balance(p,tank,record.cumulative,aquifer_influx=we or 0.)
            row.update(F=b.withdrawal.net,Eo=b.expansion.eo,Eg=b.expansion.eg,Efw=b.expansion.efw,
                       Et=b.expansion.et,We=we,injection=b.withdrawal.water_injection+b.withdrawal.gas_injection,
                       adjusted=b.withdrawal.net-we if we is not None else None)
            if we is None:
                row.update(status="CAUTION",note="Aquifer-adjusted transform unavailable after missing pressure; full observed aquifer history required.")
        except ValueError as exc:
            complete=False
            row.update(status="FAIL",note=str(exc))
        rows.append(row)
    return rows


def linear_diagnostic(x,y, *, infer="slope"):
    """Unweighted line summary only. No reservoir parameter is modified."""
    pairs = [(a,b) for a,b in zip(x,y) if a is not None and b is not None and np.isfinite(a) and np.isfinite(b)]
    result = dict(count=len(pairs),slope=None,intercept=None,r2=None,status="CAUTION",notes=[])
    if len(pairs)<3:
        result["notes"]=["At least three usable observations are required for a trend line."]
        return result
    x,y = np.array(pairs,dtype=float).T
    spread = np.ptp(x)
    if spread <= 1e-12*max(1.,np.max(abs(x))):
        result["notes"]=["Insufficient expansion/axis spread for a stable slope."]
        return result
    z=(x-x.mean())/spread
    scaled,intercept_center = np.polyfit(z,y,1)
    slope,intercept = scaled/spread,intercept_center-scaled*x.mean()/spread
    predicted=slope*x+intercept
    sse=float(np.sum((y-predicted)**2))
    sst=float(np.sum((y-y.mean())**2))
    notes=[]
    if abs(intercept)>.1*max(float(np.max(abs(y))),1e-12) and infer=="slope":
        notes.append("Large intercept (>10% of maximum transformed withdrawal).")
    inferred=slope if infer=="slope" else intercept
    if inferred<=0:
        notes.append("Nonphysical inferred OOIP (nonpositive).")
    if len(x)>=4:
        estimates=[]
        for sl in (slice(1,None),slice(None,-1)):
            if np.ptp(x[sl])>spread*1e-10:
                a,b=np.polyfit((x[sl]-x.mean())/spread,y[sl],1)
                estimates.append(a/spread if infer=="slope" else b-a*x.mean()/spread)
        if any(abs(v-inferred)>.2*max(abs(inferred),1e-12) for v in estimates):
            notes.append("Inferred OOIP changes >20% when first/last point is omitted.")
        quadratic=np.polyval(np.polyfit(z,y,2),z)
        if np.sqrt(sse/len(x))>.05*max(float(np.ptp(y)),1e-12) and np.sum((y-quadratic)**2)<.5*sse:
            notes.append("Strong curvature: quadratic residual reduction exceeds 50% with material line deviations.")
    if sst>0 and 1-sse/sst<.9:
        notes.append("Weak straight-line agreement (R² < 0.9); inspect scatter/curvature.")
    result.update(slope=float(slope),intercept=float(intercept),r2=1-sse/sst if sst>0 else None,
                  status="CAUTION" if notes else "PASS",notes=notes)
    return result


def havlena_odeh(rows, tank, variant="total"):
    gas=variant=="gas_cap"
    points=[]
    for r in rows:
        denominator = (r["Eo"]+r["Efw"]) if gas and r["Eo"] is not None else r["Et"]
        usable=denominator is not None and denominator>1e-12 and r["adjusted"] is not None
        points.append(dict(date=r["date"],x=(r["Eg"]/denominator if gas else denominator) if usable else None,
                           y=(r["adjusted"]/denominator if gas else r["adjusted"]) if usable else None))
    fit=linear_diagnostic([p["x"] for p in points],[p["y"] for p in points],infer="intercept" if gas else "slope")
    pressures=[r["pressure"] for r in rows if r["pressure"] is not None]
    if pressures and (tank.initial_pressure-min(pressures))/tank.initial_pressure<.01:
        fit["notes"].append("Less than 1% pressure depletion: inferred OOIP is poorly constrained.")
        fit["status"]="CAUTION"
    return dict(points=points,fit=fit,
        equation="(Fnet-We)/(Eo+Efw) = N + Nm·Eg/(Eo+Efw)" if gas else "Fnet-We = N·Et; Et=Eo+mEg+Efw",
        x_label="Eg/(Eo+Efw)" if gas else "Et",x_quantity="dimensionless" if gas else "expansion",
        y_label="(Fnet-We)/(Eo+Efw)" if gas else "Fnet-We",y_quantity="oil_volume" if gas else "reservoir_volume",
        inferred_from="intercept" if gas else "slope")
