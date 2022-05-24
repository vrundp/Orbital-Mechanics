# Copyright (c) 2022 Vrund Patel, USA

import numpy as np
from sklearn.metrics import r2_score

def MinimumEnergyLambertSolver(r1_vec, r2_vec, tm, MU):

    r1 = np.linalg.norm(r1_vec)
    r2 = np.linalg.norm(r2_vec)

    cosDeltaNu = np.dot(r1_vec, r2_vec) / (r1 * r2)
    sinDeltaNu = tm * np.sqrt(1 - cosDeltaNu**2)

    c = np.sqrt(r1**2 + r2**2 - 2 * r1 * r2 * cosDeltaNu)
    s = (r1 + r2 + c) / 2

    aMin = s / 2
    pMin = ((r1 * r2) / c) * (1 - cosDeltaNu)
    eMin = np.sqrt(1 - ((2 * pMin) / s))

    v1_vec = ((np.sqrt(MU * pMin)) / (r1 * r2 * sinDeltaNu)) * (r2_vec - (1 - (r2 / pMin) * (1 - cosDeltaNu)) * r1_vec)

    betaMin = 2 * np.arcsin(np.sqrt((s - c) / s))
    tofMin = np.sqrt(aMin**3 / MU) * (np.pi - tm * (betaMin - np.sin(betaMin)))

    return v1_vec



