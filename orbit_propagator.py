# Copyright (c) 2021 Vrund Patel, USA
#
import numpy as np
import matplotlib.pyplot as plt
import planetary_data as pd

from scipy.integrate import ode

class OrbitPropagator:

	def __init__(self, r0, v0, t0, tf, dt, cb = pd.earth, J2 = False, CR3BP = False):

		self.r0 = r0
		self.v0 = v0
		self.t0 = t0
		self.tf = tf
		self.dt = dt
		self.tspan = self.tf - self.t0
		self.cb = cb
		self.J2 = J2
		self.CR3BP = CR3BP

	def propagate(self):
		
		self.times = np.linspace(self.t0, self.tf, int((self.tspan / self.dt) + 1))

		self.Y0 = np.zeros(6)
		self.Y0[:3] = self.r0
		self.Y0[3:] = self.v0

		self.rv = ode(self.derivFxn)
		self.rv.set_integrator('dopri5', rtol = 1e-10, atol = 1e-20)
		self.rv.set_initial_value(self.Y0, self.t0)

		self.Y = np.zeros((len(self.times), len(self.Y0)))
		self.Y[0,:] = self.Y0

		for i in range(1, len(self.times)):

			self.Y[i, :] = self.rv.integrate(self.times[i]) 

			if not self.rv.successful():

				raise RuntimeError("Could not integrate.")

		return self.Y

	def derivFxn(self, t, Y):
    
		Y = np.array(Y)
		dY = np.zeros(Y.shape)
		    
		rx, ry, rz, vx, vy, vz = Y
		r_vec = np.array([rx, ry, rz])
		r = np.linalg.norm(r_vec)

		if self.J2 == False:

			ax, ay, az = (-self.cb["mu"] / r**3) * r_vec

		else:

			ax = (-self.cb["mu"] * rx / r**3) * (1 - self.cb["J2"] * 1.5 * (self.cb["radius"]**2 / r**2) * (5 * (rz**2 / r**2) - 1))
			ay = (-self.cb["mu"] * ry / r**3) * (1 - self.cb["J2"] * 1.5 * (self.cb["radius"]**2 / r**2) * (5 * (rz**2 / r**2) - 1))
			az = (-self.cb["mu"] * rz / r**3) * (1 - self.cb["J2"] * 1.5 * (self.cb["radius"]**2 / r**2) * (5 * (rz**2 / r**2) - 3))

		dY = np.array([vx,vy,vz,ax,ay,az])

		return dY
