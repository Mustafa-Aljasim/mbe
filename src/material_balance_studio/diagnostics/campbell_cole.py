"""Oil-reservoir Campbell transformation; dry-gas Cole is not applicable here."""

def campbell(rows):
    points=[]
    for r in rows:
        usable=r["Et"] is not None and r["Et"]>1e-12
        points.append(dict(date=r["date"],x=r["F"],unadjusted=r["F"]/r["Et"] if usable else None,
                           adjusted=r["adjusted"]/r["Et"] if usable and r["adjusted"] is not None else None))
    valid=[p for p in points if p["unadjusted"] is not None]
    cue="Insufficient observations for a shape cue."
    if len(valid)>=3:
        start,end=valid[0]["unadjusted"],valid[-1]["unadjusted"]
        change=(end-start)/max(abs(start),1e-12)
        cue="Approximately horizontal endpoint behavior." if abs(change)<=.05 else "Upward endpoint trend; review all points for systematic behavior." if change>0 else "Downward/abnormal endpoint trend; review PVT, pressure and support assumptions."
    return dict(points=points,equation="Fnet/Et = N + We/Et; adjusted: (Fnet-We)/Et = N",cue=cue)
