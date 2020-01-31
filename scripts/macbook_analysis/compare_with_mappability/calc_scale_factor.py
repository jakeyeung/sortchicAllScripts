import numpy as np
import warnings

def calc_scale_factor(X, alpha=0.5, return_nbins = False):
    """
    X: (Nbins, Nsamples)
    """
    
    assert X.ndim == 2
    assert (X.sum(axis=0) > 0).all()
    
    w = ((~np.isclose(X, 0)) & (X > 0)).any(axis=1)
    assert w.ndim == 1
    
    assert w.sum() > 0
    
    Xij = X[w]
    S = np.ones(len(Xij), dtype=bool)
    
    Xi = Xij.sum(axis=1)
    sj = Xij[S].sum(axis=0) / Xi[S].sum(axis=0)
    
    assert np.isclose(sj.sum(), 1.)
    
    Eij = np.outer(Xi, sj)
    GOFi = (((Xij - Eij) ** 2) / Eij).sum(axis=1)
    assert len(GOFi) == len(S)
    
    _GOF_low, _GOF_high = np.quantile(GOFi, np.array([0.05, min(1.0 - 0.05, 1 - alpha + 0.05)]))
    
    for _ in range(20):
        S_update = (GOFi >= _GOF_low) & (GOFi < _GOF_high)
        if (S == S_update).all():
            break

        S = S_update.copy()

        sj = Xij[S].sum(axis=0) / Xi[S].sum(axis=0)
        
        Eij = np.outer(Xi, sj)
        GOFi = (((Xij - Eij) ** 2) / Eij).sum(axis=1)
        
        _GOF_low, _GOF_high = np.quantile(GOFi, np.array([0.05, min(1.0 - 0.05, 1 - alpha + 0.05)]))
    else:
        warnings.warn("Did not converge")
    
    sj_final = sj.copy()
    # print("Bins considered:" + str(len(S)))
    if return_nbins:
      return sum(S)
    else:
      return sj_final
