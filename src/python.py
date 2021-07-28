import numpy as np
import matplotlib.pyplot as plt
import numba as nb

# AU, years
μ = 39.4234021    # gravitational parameter G*M_sun. AU^3/Years^2
a = 0.38709893    # semimajor axis
α = 1.09778201e-8 # units of AU^2, is 3L^2/(mc)^2, where m is reduced mass, L angular momentum
e = 0.20563069    # eccentricity


@nb.njit("f8(f8[:])",fastmath=True)
def angle(point):
	return np.mod(np.arctan2(point[1],point[0]),2*np.pi)

@nb.njit("f8[:](f8[:])",fastmath=True)
def acceleration(r_vec):
	r = np.linalg.norm(r_vec)
	β = -μ/r**3 * (1+α/r**2)
	return r_vec*β

@nb.njit("none(f8,f8,f8)",fastmath=True)
def verlet(t_max,dt_max,dt_min):
	#results = TypedList()
	
	r_past = np.array([a*(1+e),0.])
	v_past = np.array([0.,np.sqrt((μ/a)*(1-e)/(1+e))])
	a_past = acceleration(r_past)

	r_now = np.zeros(2)
	v_now = np.zeros(2)
	a_now = np.zeros(2)
	
	event_past = np.dot(r_past,v_past)
	event_now = 0.0
	
	dt = dt_max
	t = 0.
	while t < t_max:
		r_now = r_past + v_past*dt + a_past*dt**2/2
		a_now = acceleration(r_now)
		v_now = v_past + dt*(a_now+a_past)/2
		event_now = np.dot(r_now,v_now)
		
		if event_past<0 and event_now>0:
			if dt<dt_min:
				print(t,np.rad2deg(angle((r_now+r_past)/2)-np.pi)*3600)
				dt = dt_max
			else:
				# rewind time
				t -= dt
				r_now = r_past
				v_now = v_past
				a_now = a_past
				event_now = event_past
				# adjust dt
				dt /= 2
		
		r_past = r_now
		v_past = v_now
		a_past = a_now
		event_past = event_now
		t += dt

verlet(10.0,1e-6,1e-12)

