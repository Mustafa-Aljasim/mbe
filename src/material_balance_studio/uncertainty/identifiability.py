"""Engineering classification uses local information, bounds and profile evidence."""
import numpy as np


def classify(specs,values,scaled_sensitivity,diagnostics,nobs,pressure_spread,profiles=()):
    result=[]
    sensitivity=np.sqrt(np.mean(np.asarray(scaled_sensitivity)**2,axis=0))
    for j,(s,x) in enumerate(zip(specs,values)):
        notes=[]
        bound=min(x-s.lower,s.upper-x)/(s.upper-s.lower)<1e-4
        corr=diagnostics["correlation"]
        coupling=max((abs(corr[j,k]) for k in range(len(specs)) if k!=j),default=0) if corr is not None else 1.
        deficient=diagnostics["rank"]<len(specs)
        if bound:
            status="BOUND LIMITED"
            notes.append("User bound limits inference; uncertainty outside the bound is not estimated.")
        elif deficient or sensitivity[j]<1e-8:
            status="NON-IDENTIFIABLE"
            notes.append("Insufficient independent local information; no reliable marginal interval.")
        elif coupling>.9 or sensitivity[j]<1e-5 or diagnostics["condition"]>1e4:
            status="POORLY CONSTRAINED"
            notes.append("Low sensitivity and/or strong local coupling/conditioning; inspect compensating profiles.")
        elif coupling>=.7 or nobs<max(5,2*len(specs)+1) or pressure_spread<.01:
            status="MODERATELY CONSTRAINED"
            notes.append("Limited pressure information or moderate parameter coupling.")
        else:
            status="WELL CONSTRAINED"
            notes.append("Local sensitivity, rank, coupling and data-spread checks pass; global uniqueness is not established.")
        for profile in profiles:
            if profile.parameter==s.name and profile.plausible_ranges:
                width=max(b-a for a,b in profile.plausible_ranges)/(s.upper-s.lower)
                if width>.5 and status not in ("BOUND LIMITED","NON-IDENTIFIABLE"):
                    status="POORLY CONSTRAINED"
                    notes.append("Profile permits more than half the allowed parameter span.")
                if len(profile.plausible_ranges)>1:
                    notes.append("Profile contains disconnected grid-resolved plausible regions.")
        result.append((s.name,status," ".join(notes)))
    return tuple(result)
