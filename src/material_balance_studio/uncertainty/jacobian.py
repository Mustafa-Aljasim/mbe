"""SVD and covariance geometry on a parameter-scaled residual Jacobian."""
import numpy as np


def matrix_diagnostics(jacobian,relative_tolerance=1e-6):
    a=np.asarray(jacobian,dtype=float)
    if a.ndim!=2 or not a.size or not np.all(np.isfinite(a)):
        raise ValueError("Jacobian must be a finite nonempty matrix.")
    _,singular,vt=np.linalg.svd(a,full_matrices=False)
    threshold=max(1e-10,relative_tolerance*singular[0])
    rank=int(np.sum(singular>threshold))
    norms=np.linalg.norm(a,axis=0)
    denominator=np.outer(norms,norms)
    coupling=np.divide(a.T@a,denominator,out=np.zeros_like(denominator),where=denominator>0)
    geometry=None
    correlation=None
    if rank==a.shape[1]:
        geometry=(vt.T*(1/singular**2))@vt
        sd=np.sqrt(np.diag(geometry))
        correlation=np.clip(geometry/np.outer(sd,sd),-1,1)
    return dict(singular=singular,rank=rank,condition=float(singular[0]/singular[-1]) if rank==a.shape[1] else None,
                weak=vt[-1],geometry=geometry,correlation=correlation,coupling=np.clip(coupling,-1,1),norms=norms)
