# Example Molecular Dynamics Code
# For understanding MD / Practice with Python, Github, etc.
# Written by: Brian Day

# Import necessary packages
import copy
import csv
import numpy as np
import random as rand

# Internal System of Units
# Time (t') = 1 s
# Length (l') = 1x10^-10 m = 1 Angstrom
# Mass (m') = 1.660539x10^-27 kg = 1 Da = 1 u
	# This is simply inverse of Avogadro's number.
# Temperature = Kelvin

# New Boltmann Coefficient
k_B = 1.38064852e-23 	# [m^2 kg s^-2 K^-1]
k_Bp = 8.314459920816467e23 	# [A^2 Da s^-2 K^-1]
	# Universal gas constant with differnet units

# Define functions
def load_positions(filename):
	"""
	Load the particle types and positions as a dictionary by atom number.
	"""

	# Load the positions file.
	with open(filename,newline='') as csvfile:
		output_data = csv.reader(csvfile, delimiter="\t")
		output_data = list(output_data)

	# Initialize empty positions list, parse from file.
		positions_list = []
		for i in range(1,len(output_data)):
			row = output_data[i]
			temp_dict = {}
			temp_dict['num'] = i
			temp_dict['type'] = row[0]
			temp_dict['pos_vect'] = np.array([float(item) for item in row[1:4]])
			temp_dict['vel_vect'] = np.array([0.0, 0.0, 0.0])
			temp_dict['force_vect'] = np.array([0.0, 0.0, 0.0])
			temp_dict['KE'] = 0.0
			positions_list.append(temp_dict)

		return positions_list


def load_parameters(filename):
	"""
	Load the Lennard-Jones (and other relevant parameters as updated) to be used in calculating the
	forces.
	"""

	# Load the parameters file.
	with open(filename,newline='') as csvfile:
		output_data = csv.DictReader(csvfile, delimiter="\t")
		output_data = list(output_data)

	# Initialize empty parameters dictionary.
		params_dict = {}
		for row in output_data:
			type = row['type']
			del row['type']
			keys = row.keys()
			for key in keys:
				row[key] = float(row[key])
			params_dict[type] = dict(row)

		return params_dict


def create_random_unit_vector():
	"""
	Creates a random unit vector.
	"""
	vector = np.array([rand.uniform(-1,1), rand.uniform(-1,1), rand.uniform(-1,1)])
	vector_mag = np.linalg.norm(vector)
	unit_vector = np.array([val/vector_mag for val in vector])
	return unit_vector


def initialize_MD(temperature, positions, parameters, distribution='Boltzmann'):
	"""
	Takes in the temperture to define a set of inital velocity vectors. Velcities can either be
	assigned to have an (approximately) uniform distribtuion or an (approximately) Boltzmann
	distribution. Simulations should still include initialization cycles before any properties are
	calculated.
	Velocities are normalized such that the total momentum is zero (i.e. No externl forces).
	"""

	# Assign initial normalized velocities of a given distribution type.
	# It can be shown that the Boltzmann Distribution is simply the sum of three normal
	# distributions of the momentum with variance mkT, where m i the molecular mass, k is the
	# Boltzmann constant, and T is the temperature. This property can be leverage during
	# initialization.
	k_Bp = 8.314459920816467e23
	momentum_total = np.array([0.0,0.0,0.0])

	if distribution == 'Uniform':
		for particle in positions:
			type = particle['type']
			mass = parameters[type]['mass']
			particle['vel_vect'] = rand.uniform(0,1) * create_random_unit_vector()
			momentum_total += [ mass*vel for vel in particle['vel_vect'] ]

	elif distribution == 'Boltzmann':
		for particle in positions:
			type = particle['type']
			mass = parameters[type]['mass']
			sigma = np.sqrt(mass*k_Bp*temperature)
			particle['vel_vect'] = np.array([val/mass for val in np.random.normal(0,sigma,3)])
			momentum_total += [ mass*vel for vel in particle['vel_vect'] ]

	else:
		print('Invalid Distribution Type!')

	# Shift All Vectors so that Total Momentum is 0 (i.e. No External Forces)
	num_particles = len(positions)
	velocity_total = np.array([0.0,0.0,0.0])
	for particle in positions:
		type = particle['type']
		mass = parameters[type]['mass']
		particle['vel_vect'] += -1*momentum_total/mass/num_particles
		particle['KE'] = 0.5*mass*np.linalg.norm(particle['vel_vect'])**2
		velocity_total += particle['vel_vect']

	# Scale All Vectors so that the initial Temperature is Correct
	# Does the total momentum in x need to equal the total momentum in y and z, or just overall total is 0?
	return positions, velocity_total


def apply_mixing_rules(sigma1, sigma2, eps1, eps2, mixing_rules=None):
	"""
	Apply mixing rules on the sigma/epsilon parameters. Written as a separate function rather than
	part of the calculate_forces code so that multiple mixing_rules can be stored.
	"""
	if mixing_rules == None or mixing_rules == 'Lorentz-Berthelot' or mixing_rules == 'LB':
		sigma_mix = 0.5*(sigma1+sigma2)
		epsilon_mix = np.sqrt(eps1*eps2)
	else:
		print('Invalid Mixing Rule!')
		sigma_mix = 0
		epsilon_mix = 0

	return sigma_mix, epsilon_mix


def replicate_cell(positions, box_size):

	if box_size != None:
		box_x = box_size[0]
		box_y = box_size[1]
		box_z = box_size[2]
	else:
		print('Box size needed to replicate cells!')

	x_trans = [0, box_x, -1*box_x]
	y_trans = [0, box_y, -1*box_y]
	z_trans = [0, box_z, -1*box_z]

	positions_with_replicas = []

	for x_scalar in x_trans:
		for y_scalar in y_trans:
			for z_scalar in z_trans:
				for position in positions:
					position_new = copy.deepcopy(position)
					x_pos_new = position_new['pos_vect'][0]+x_scalar
					y_pos_new = position_new['pos_vect'][1]+y_scalar
					z_pos_new = position_new['pos_vect'][2]+z_scalar
					position_new['pos_vect'] = [x_pos_new, y_pos_new, z_pos_new]
					positions_with_replicas.extend([position_new])

	return positions_with_replicas


def calculate_forces(positions, positions_with_replicas, parameters, r_cutoff):
	"""
    Evaluates Lennard-Jones potentials for all particles and assigns forces.
    Update later to include charge-charge interactions.
    """
	k_Bp = 8.314459920816467e23

	for i in range(len(positions)):
		positions[i]['force_vect'] = np.array([0.0,0.0,0.0])

	PE_total = 0.0
	for i in range(len(positions)):
		for j in range(i+1,len(positions_with_replicas)):
			r_vect = positions[i]['pos_vect'] - positions_with_replicas[j]['pos_vect']
			r_mag = np.linalg.norm(r_vect)
			if r_mag <= r_cutoff:
				type1 = positions[i]['type']
				type2 = positions_with_replicas[j]['type']
				sigma1 = parameters[type1]['sigma']
				sigma2 = parameters[type2]['sigma']
				eps1 = parameters[type1]['epsilon']
				eps2 = parameters[type2]['epsilon']
				sigma_mix, epsilon_mix = apply_mixing_rules(sigma1, sigma2, eps1, eps2, mixing_rules='LB')

				LJ_potential = 4*epsilon_mix*k_Bp*( (sigma_mix/r_mag)**12 - (sigma_mix/r_mag)**6 )
				LJ_force = -4*epsilon_mix*k_Bp*( -12*(sigma_mix**12/r_mag**13) + 6*(sigma_mix**6/r_mag**7))
				PE_total += LJ_potential

				positions[i]['force_vect'] = positions[i]['force_vect'] + [LJ_force*val for val in r_vect/r_mag]
				if j < len(positions):
					positions[j]['force_vect'] = positions[j]['force_vect'] + [-1*LJ_force*val for val in r_vect/r_mag]

	return PE_total


def integrate_forces(positions, parameters, timestep):
	"""
	Uses the calculate forces and the time step to update the positions of the particles in the
	system.
	"""
	for particle in positions:

		type = particle['type']
		mass = parameters[type]['mass']
		vel_init = particle['vel_vect']
		pos_init = particle['pos_vect']

		accel = [force/mass for force in particle['force_vect']]
		vel_final = [vel_init[i] + accel[i]*timestep for i in range(len(vel_init))]
		pos_final = [vel_final[i]*timestep + pos_init[i] for i in range(len(vel_init))]
		particle['vel_vect'] = np.array(vel_final)
		particle['pos_vect'] = np.array(pos_final)
		particle['KE'] = 0.5*mass*np.linalg.norm(particle['vel_vect'])**2


def integrate_forces_verlet(positions, parameters, timestep, box_size, r_cutoff):

	for particle in positions:

		type = particle['type']
		mass = parameters[type]['mass']
		vel_init = particle['vel_vect']
		pos_init = particle['pos_vect']

		accel = [force/mass for force in particle['force_vect']]
		vel_half  = [vel_init[i] + accel[i]*timestep/2 for i in range(len(vel_init))]
		pos_final = [pos_init[i] + vel_half[i]*timestep for i in range(len(vel_half))]
		particle['vel_half'] = np.array(vel_half)
		particle['pos_vect'] = np.array(pos_final)

	check_boundaries(positions, box_size)
	# positions_with_replicas = replicate_cell(positions, box_size)
	PE_total = calculate_forces(positions, positions, parameters, r_cutoff)

	for particle in positions:

		type = particle['type']
		mass = parameters[type]['mass']
		vel_half = particle['vel_half']

		accel_final = [force/mass for force in particle['force_vect']]
		vel_final = [vel_half[i] + accel_final[i]*timestep/2 for i in range(len(vel_half))]
		particle['vel_vect'] = np.array(vel_final)
		particle['KE'] = 0.5*mass*np.linalg.norm(particle['vel_vect'])**2

	return positions


def check_boundaries(positions, box_size):
	"""
	Check if any particle has gone beyond the boundary, and if so update position through periodic
	boundary conditions.
	"""
	if box_size != None:
		box_x = box_size[0]
		box_y = box_size[1]
		box_z = box_size[2]
	else:
		print('Box size needed to apply periodic boundary conditions!')

	for particle in positions:
		# Check x-positiions
		if particle['pos_vect'][0] < 0:
			particle['pos_vect'][0] = box_x + particle['pos_vect'][0];
		if particle['pos_vect'][0] > box_x:
			particle['pos_vect'][0] = particle['pos_vect'][0] - box_x;

		# Check x-positiions
		if particle['pos_vect'][1] < 0:
			particle['pos_vect'][1] = box_y + particle['pos_vect'][1];
		if particle['pos_vect'][1] > box_y:
			particle['pos_vect'][1] = particle['pos_vect'][1] - box_y;

		# Check x-positiions
		if particle['pos_vect'][2] < 0:
			particle['pos_vect'][2] = box_z + particle['pos_vect'][2];
		if particle['pos_vect'][2] > box_z:
			particle['pos_vect'][2] = particle['pos_vect'][2] - box_z;

	return positions


def modulo(A,B):
	return int(A) % int(B)


def log_positions(positions, cycle_num, filename):
	if cycle_num == 0:
		csvfile = open(filename,'w', newline='')
	else:
		csvfile = open(filename,'a', newline='')
	writer = csv.writer(csvfile, delimiter="\t")
	writer.writerow([str(len(positions))])
	writer.writerow(['Cycle Number = '+str(cycle_num)])
	for particle in positions:
		type = particle['type']
		xyz = ' '.join(map(str, particle['pos_vect']))
		writer.writerow([type, xyz])


def log_properties(positions, cycle_num, total_PE, filename):
	if cycle_num == -1:
		csvfile = open(filename,'w', newline='')
	else:
		csvfile = open(filename,'a', newline='')
	writer = csv.writer(csvfile, delimiter="\t")

	total_KE = 0
	for particle in positions:
		total_KE += particle['KE']
	total_E = total_KE + total_PE

	writer.writerow(['Cycle Number = '+str(cycle_num)])
	writer.writerow(['     Kinetic Energy = '+str(np.format_float_scientific(total_KE, precision=4))])
	writer.writerow(['     Potential Energy = '+str(np.format_float_scientific(total_PE, precision=4))])
	writer.writerow(['     Total Energy = '+str(np.format_float_scientific(total_E, precision=4))])


def run_MD(positions_file, parameters_file, log_filename, init_cycles, production_cycles, temperature, box_size, r_cutoff, time_step, log_freq, log_prop_filename):
	positions = load_positions(positions_file)
	parameters = load_parameters(parameters_file)
	initialize_MD(temperature, positions, parameters)
	positions_with_replicas = replicate_cell(positions, box_size)
	PE_total = calculate_forces(positions, positions_with_replicas, parameters, r_cutoff)

	log_properties(positions, -1, PE_total, log_prop_filename)

	for i in range(init_cycles):
		cycle_num = i
		# integrate_forces(positions, parameters, timestep)
		integrate_forces_verlet(positions, parameters, timestep, box_size, r_cutoff)
		check_boundaries(positions, box_size)
		positions_with_replicas = replicate_cell(positions, box_size)
		PE_total = calculate_forces(positions, positions_with_replicas, parameters, r_cutoff)
		if modulo(cycle_num, log_freq) == 0:
			log_positions(positions, cycle_num, log_filename)
			log_properties(positions, cycle_num, PE_total, log_prop_filename)


# Exeute Program
pos_file = 'test_pos.xyz'
params_file = 'LJ_params.def'
log_file = '/Users/brian_day/Desktop/log_file_text.xyz'
log_prop_file = '/Users/brian_day/Desktop/log_prop_text.txt'
init_cycles = 2000
prod_cycles = 200
temperature = 300 # Kelvin
box_size = [20,20,20] # Angstrom
r_cutoff = 20 # Angstrom
timestep = 1e-16 #Seconds
log_freq = 20

run_MD(pos_file, params_file, log_file, init_cycles, prod_cycles, temperature, box_size, r_cutoff, timestep, log_freq, log_prop_file)
