# Copyright (c) 2022 Vrund Patel, USA
#
# This file contains constsants and enumerations.
# Data is extracted from Appendix A of "Orbital Mechanics for Engineering Students".

'''

Units:

mu: km^3/s^2
mass: kg
radius: km
sma: km
period: s
soi: km

'''

sun = {
	
	"name": "Sun",
	"mu": 132.712e9,
	"mass": 1.989e30,
	"radius": 696000

}

mercury = {
	
	"name": "Mercury",
	"mu": 22030,
	"mass": 330.2e21,
	"radius": 2440,
	"sma": 57.91e6,
	"ecc": 0.2056,
	"period": 87.97*24*3600,
	"soi": 112000

}

venus = {
	
	"name": "Venus",
	"mu": 324900,
	"mass": 4.869e24,
	"radius": 6052,
	"sma": 108.2e6,
	"ecc": 0.0067,
	"period": 224.7*24*3600,
	"soi": 616000

}

earth = {
	
	"name": "Earth",
	"mu": 398600.4415,
	"mass": 5.974e24,
	"radius": 6378.1363,
	"sma": 149.6e6,
	"ecc": 0.0167,
	"period": 365.256*24*3600,
	"soi": 925000,
	"J2": 0.00108248

}

moon = {
	
	"name": "Moon",
	"mu": 4903,
	"mass": 73.48e21,
	"radius": 1737,
	"sma": 384.4e3,
	"ecc": 0.0549,
	"period": 27.322*24*3600,
	"soi": 66200

}

mars = {
	
	"name": "Mars",
	"mu": 42828,
	"mass": 641.9e21,
	"radius": 3396,
	"sma": 227.9e6,
	"ecc": 0.0935,
	"period": 1.881*365.256*24*3600,
	"soi": 577000

}

jupiter = {
	
	"name": "Jupiter",
	"mu": 126.686e6,
	"mass": 1.899e27,
	"radius": 71490,
	"sma": 778.6e6,
	"ecc": 0.0489,
	"period": 11.86*365.256*24*3600,
	"soi": 48.2e6

}

saturn = {
	
	"name": "Saturn",
	"mu": 37.931e6,
	"mass": 568.5e24,
	"radius": 60270,
	"sma": 1.433e9,
	"ecc": 0.0565,
	"period": 29.46*365.256*24*3600,
	"soi": 54.8e6

}

uranus = {
	
	"name": "Uranus",
	"mu": 5.794e6,
	"mass": 86.83e24,
	"radius": 25560,
	"sma": 2.872e9,
	"ecc": 0.0457,
	"period": 84.01*365.256*24*3600,
	"soi": 51.8e6

}

neptune = {
	
	"name": "Neptune",
	"mu": 6.8351e6,
	"mass": 102.4e24,
	"radius": 24760,
	"sma": 4.495e9,
	"ecc": 0.0113,
	"period": 164.8*365.256*24*3600,
	"soi": 86.6e6

}

pluto = {
	
	"name": "Pluto",
	"mu": 830,
	"mass": 12.5e21,
	"radius": 1195,
	"sma": 5.870e9,
	"ecc": 0.2444,
	"period": 247.7*365.256*24*3600,
	"soi": 3.08e6

}

