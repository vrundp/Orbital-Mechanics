# Copyright (c) 2022 Vrund Patel, USA

import numpy as np

def HerrickOrGibbs(r1_vec, r2_vec, r3_vec, t1, t2, t3, MU):

	r1 = np.linalg.norm(r1_vec)
	r2 = np.linalg.norm(r2_vec)
	r3 = np.linalg.norm(r3_vec)

	r1_hat = r1_vec / r1
	r2_hat = r2_vec / r2
	r3_hat = r3_vec / r3

	a12 = np.rad2deg( np.arccos( np.dot(r1_hat, r2_hat) ) )
	a23 = np.rad2deg( np.arccos( np.dot(r2_hat, r3_hat) ) )

	if a12 < 1 or a23 < 1:

		HerrickGibbs(r1_vec, r2_vec, r3_vec, t1, t2, t3, MU)

	elif a12 > 5 and a23 > 5:

		Gibbs(r1_vec, r2_vec, r3_vec, MU)
	
	else:

		return "This IOD method will not work."


def Gibbs(r1_vec, r2_vec, r3_vec, MU):

	r1 = np.linalg.norm(r1_vec)
	r2 = np.linalg.norm(r2_vec)
	r3 = np.linalg.norm(r3_vec)

	Z12 = np.cross(r1_vec, r2_vec)
	Z23 = np.cross(r2_vec, r3_vec)
	Z31 = np.cross(r3_vec, r1_vec)

	N = r1 * Z23 + r2 * Z31 + r3 * Z12

	D = Z12 + Z23 + Z31

	S = (r2 - r3) * r1_vec + (r3 - r1) * r2_vec + (r1 - r2) * r3_vec

	B = np.cross(D, r2_vec)

	Lg = np.sqrt( MU / (np.linalg.norm(N) * np.linalg.norm(D)) )

	v2_vec = (Lg / r2) * B + Lg * S

	return v2_vec

def HerrickGibbs(r1_vec, r2_vec, r3_vec, t1, t2, t3, MU):

	r1 = np.linalg.norm(r1_vec)
	r2 = np.linalg.norm(r2_vec)
	r3 = np.linalg.norm(r3_vec)

	t31 = t3 - t1
	t32 = t3 - t2
	t21 = t2 - t1

	coeff1 = -t32 * ( ( 1 / (t21 * t31) ) + ( MU / (12 * r1**3)) )
	coeff2 = (t32 - t21) * ( ( 1 / (t21 * t32) ) + ( MU / (12 * r2**3)) )
	coeff3 = t21 * ( ( 1 / (t32 * t31) ) + ( MU / (12 * r3**3)) )

	v2_vec = coeff1 * r1_vec + coeff2 * r2_vec + coeff3 * r3_vec

	return v2_vec


