# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 22:51:22 2024

@author: angel
"""

import numpy as np
from scipy.sparse import spdiags, csr_matrix

def discret(aP, aE, aW, aS, aN, aB, aT, Su):
    K = csr_matrix((np.prod(aP.shape), np.prod(aP.shape)))
    aP = aP.flatten()
    tam = spdiags(aP, 0, K.shape[0], K.shape[1])
    K = K + tam

    aE = -aE.flatten()
    aE = aE[1:]
    tam = spdiags(aE, -1, K.shape[0], K.shape[1])
    K = K + tam

    aW = -aW.flatten()
    tam = spdiags(aW, -1, K.shape[0], K.shape[1]).T
    K = K + tam

    aS = -aS.flatten()
    n = aS.shape[0]
    aS = aS[n:]
    tam = spdiags(aS, -n, K.shape[0], K.shape[1])
    K = K + tam

    aN = -aN.flatten()
    n = aN.shape[0]
    tam = spdiags(aN, -n, K.shape[0], K.shape[1]).T
    K = K + tam

    aB = -aB.flatten()
    n = np.prod(aB.shape[:2])
    aB = aB[n:]
    tam = spdiags(aB, -n, K.shape[0], K.shape[1])
    K = K + tam

    aT = -aT.flatten()
    n = np.prod(aT.shape[:2])
    tam = spdiags(aT, -n, K.shape[0], K.shape[1]).T
    K = K + tam

    F = csr_matrix(Su.flatten())

    return K, F
