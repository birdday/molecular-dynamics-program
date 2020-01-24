"""
Test of a 1-D Molecuar Dynamic Problem with a Verlet Integration Scheme.
"""

import numpy as np
import matplotlib.pyplot as plt

# r_vect = np.linspace(2.5,5.5,51)
# V_vect = 4*eps*k_Bp*( (sigma/r_vect)**12 - (sigma/r_vect)**6)
# F_vect = 4*eps*k_Bp*( -12*(sigma**12/r_vect**13) + 6*(sigma**6/r_vect**7))
# plt.plot(r_vect, V_vect)


# Initial Positions & Velocities
r1 = 4.0
r2 = 7.0
v1 = 0
v2 = 0

# Molecule Parameters - Helium
sigma = 2.628
eps = 5.465
m = 2

# Calculate Potential and Force
dt = 1e-15
k_Bp = 8.314459920816467e23

algorithm = 'Verlet'
for i in range(2000):
	if algorithm == 'Verlet':
		r = r1-r2
		F1 = -4*eps*k_Bp*( -12*(sigma**12/r**13) + 6*(sigma**6/r**7))
		F2 = -F1
		a1 = F1/m
		a2 = F2/m
		v1_half = v1 + a1*dt/2
		v2_half = v2 + a2*dt/2
		r1 = r1 + v1_half*dt
		r2 = r2 + v2_half*dt
		r = r1-r2
		V = 4*eps*k_Bp*( (sigma/r)**12 - (sigma/r)**6 )
		F1 = -4*eps*k_Bp*( -12*(sigma**12/r**13) + 6*(sigma**6/r**7))
		F2 = -F1
		a1 = F1/m
		a2 = F2/m
		v1 = v1_half + a1*dt/2
		v2 = v2_half + a2*dt/2

	if algorithm == 'Simple':
		if r1 < 0:
			r1 = r1+10
		if r2 > 10:
			r1 = r1-10
		if r2 < 0:
			r2 = r2+10
		if r2 > 10:
			r2 = r2-10
		r = r1-r2
		F1 = -4*eps*k_Bp*( -12*(sigma**12/r**13) + 6*(sigma**6/r**7))
		F2 = -F1
		a1 = F1/m
		a2 = F2/m
		V = 4*eps*k_Bp*( (sigma/r)**12 - (sigma/r)**6 )
		v1 = v1 + a1*dt
		v2 = v2 + a2*dt
		r1 = r1 + v1*dt
		r2 = r2 + v2*dt
		# v1 = v1_temp
		# v2 = v2_temp

	if i%20 == 0:
		KE = 0.5*m*v1**2 + 0.5*m*v2**2
		TE = V + KE
		# print('TE ='+np.format_float_scientific(TE,6))
		print('Cycle: '+str(i)+'\tr = '+str(np.round(r1,3))+'\tTE ='+np.format_float_scientific(TE,4))
		# print('r = '+str(np.round(r,3))+'\tV = '+np.format_float_scientific(V,precision=3)+'\tKE = '+np.format_float_scientific(KE,precision=3))
	# print(V)
