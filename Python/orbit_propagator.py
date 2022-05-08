# Copyright (c) 2022 Vrund Patel, USA
#
import numpy as np
import matplotlib.pyplot as plt
import planetary_data as pd

from scipy.integrate import ode
from mpl_toolkits.mplot3d import Axes3D

class OrbitPropagator:

	def __init__(self, r0, v0, t0, tf, dt, cb = pd.earth, J2 = False, CR3BP = False, m1 = pd.earth["mass"], m2 = pd.moon["mass"]):

		self.r0 = r0
		self.v0 = v0
		self.t0 = t0
		self.tf = tf
		self.dt = dt
		self.tspan = self.tf - self.t0
		self.cb = cb
		self.J2 = J2
		self.CR3BP = CR3BP
		self.m1 = m1
		self.m2 = m2

	def propagate(self):

		self.times = np.linspace(self.t0, self.tf, int((self.tspan / self.dt) + 1))

		self.Y0 = np.zeros(6)
		self.Y0[:3] = self.r0
		self.Y0[3:] = self.v0

		self.rv = ode(self.derivFunc)
		self.rv.set_integrator('dopri5', rtol = 1e-10, atol = 1e-20)
		self.rv.set_initial_value(self.Y0, self.t0)

		self.Y = np.zeros((len(self.times), len(self.Y0)))
		self.Y[0,:] = self.Y0

		for i in range(1, len(self.times)):

			self.Y[i, :] = self.rv.integrate(self.times[i]) 

			if not self.rv.successful():

				raise RuntimeError("Could not integrate.")

		return self.Y

	def derivFunc(self, t, Y):

		Y = np.array(Y)
		dY = np.zeros(Y.shape)
		    
		rx, ry, rz, vx, vy, vz = Y
		r_vec = np.array([rx, ry, rz])
		r = np.linalg.norm(r_vec)

		if self.J2 == True:

			self.CR3BP = False

			ax = (-self.cb["mu"] * rx / r**3) * (1 - self.cb["J2"] * 1.5 * (self.cb["radius"]**2 / r**2) * (5 * (rz**2 / r**2) - 1))
			ay = (-self.cb["mu"] * ry / r**3) * (1 - self.cb["J2"] * 1.5 * (self.cb["radius"]**2 / r**2) * (5 * (rz**2 / r**2) - 1))
			az = (-self.cb["mu"] * rz / r**3) * (1 - self.cb["J2"] * 1.5 * (self.cb["radius"]**2 / r**2) * (5 * (rz**2 / r**2) - 3))

		if self.CR3BP == True:

			self.J2 = False
			
			mu = self.m2 / (self.m1 + self.m2)
			r1 = np.sqrt((rx + mu)**2 + ry**2 + rz**2)
			r2 = np.sqrt((rx + mu - 1)**2 + ry**2 + rz**2)


			ax = 2 * vy + rx - (1 - mu) * (rx + mu) / (r1**3) - mu * (rx + mu - 1) / (r2**3)
			ay = -2 * vx + ry - (1 - mu) * (ry) / (r1**3) - mu * (ry) / (r2**3)
			az = -(1 - mu) * (rz) / (r1**3) - mu * (rz) / (r2**3)

		else:

			ax, ay, az = (-self.cb["mu"] / r**3) * r_vec

		dY = np.array([vx, vy, vz, ax, ay, az])

		return dY

	def plot(self):
		
		pos = np.zeros((len(self.times), 1))
		for i in range(len(pos)):

			pos[i] = np.sqrt(self.Y[i, 0]**2 + self.Y[i, 1]**2 + self.Y[i, 2]**2)

		posx = np.zeros((len(self.times), 1))
		for i in range(len(posx)):

			posx[i] = self.Y[i, 0]

		posy = np.zeros((len(self.times), 1))
		for i in range(len(posy)):

			posy[i] = self.Y[i, 1]

		posz = np.zeros((len(self.times), 1))
		for i in range(len(posz)):

			posz[i] = self.Y[i, 2]

		vel = np.zeros((len(self.times), 1))
		for i in range(len(vel)):

			vel[i] = np.sqrt(self.Y[i, 3]**2 + self.Y[i, 4]**2 + self.Y[i, 5]**2)

		acc_vec = np.zeros((len(self.times), 3))
		for i in range(len(self.times)):

			temp = self.derivFunc(self.times[i], self.Y[i, :])
			acc_vec[i, :] = temp[3:]
		
		acc = np.zeros((len(acc_vec), 1))
		for i in range(len(acc)):

			acc[i] = np.sqrt(acc_vec[i, 0]**2 + acc_vec[i, 1]**2 + acc_vec[i, 2]**2)
		
		energy = np.zeros((len(self.times), 1))
		for i in range(len(self.times)):

			energy[i] = (vel[i]**2 / 2) - self.cb["mu"] / pos[i]


		fig = plt.figure(figsize = (10, 10))
		ax = fig.add_subplot(111, projection = '3d')

		u = np.linspace(0, 2 * np.pi, 100)
		v = np.linspace(0, np.pi, 100)

		x = self.cb["radius"] * np.outer(np.cos(u), np.sin(v))
		y = self.cb["radius"] * np.outer(np.sin(u), np.sin(v))
		z = self.cb["radius"] * np.outer(np.ones(np.size(u)), np.cos(v))

		ax.plot_surface(x, y, z, color = 'blue', label = self.cb["name"])
		ax.scatter3D(posx, posy, posz, color = "black", label = "Trajectory")
		ax.set_xlabel("x position", fontsize = 14)
		ax.set_ylabel("y position", fontsize = 14)
		ax.set_zlabel("z position", fontsize = 14)

		plt.show()
