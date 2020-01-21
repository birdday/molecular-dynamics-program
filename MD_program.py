# Example Molecular Dynamics Code
# For understanding MD / Practice with Python, Github, etc.
# Written by: Brian brian_day

# Import necessary packages
import copy
import csv
import numpy as np
import random as rand


# Define functions
def load_positions(filename, box_size=None):
	"""
	Load the particle types and positions as a dictionary by atom number.
	"""
	with open(filename,newline='') as csvfile:
		output_data = csv.reader(csvfile, delimiter="\t")
		output_data = list(output_data)

		positions_list = []

		if box_size != None:
			box_x = box_size[0]
			box_y = box_size[1]
			box_z = box_size[2]

		for i in range(1,len(output_data)):
			row = output_data[i]
			temp_dict = {}
			temp_dict['num'] = i
			temp_dict['type'] = row[0]
			if len(row) > 1:
				temp_dict['pos_vect'] = np.array([float(item) for item in row[1:4]])
			elif box_size != None:
				temp_dict['pos_vect'] = np.array([rand.uniform(0,box_x), rand.uniform(0,box_y), rand.uniform(0,box_z)])
			else:
				print('Box size needed to apply random coordinates!')

			temp_dict['vel_vect'] = ''
			temp_dict['force_vect'] = np.array([0,0,0])
			positions_list.append(temp_dict)
		return positions_list


def load_parameters(filename):
	"""
	Load the Lennard-Jones (and other relevant parameters as updated) to be used in calculating the
	forces.
	"""
	with open(filename,newline='') as csvfile:
		output_data = csv.DictReader(csvfile, delimiter="\t")
		output_data = list(output_data)

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
	unit_vector = [val/vector_mag for val in vector]
	return unit_vector


def initialize_MD(temperature, temp_range, positions, parameters):
	"""
	Takes in the temperture to define a set of inital velocity vectors. Speeds of the particles are
	assigned by creating a list of temperatures whose average equals the input temperature. Current
	approch does NOT guarantee a Boltzmann distribution from the start. CHECK UNITS!!
	"""
	list_of_temperatures = []
	list_of_speeds = []

	for particle in positions:
		type = particle['type']
		mass = parameters[type]['mass']
		temp = rand.uniform(temperature-temp_range, temperature+temp_range)
		list_of_temperatures.extend([temp])
		speed = np.sqrt(2*temp/mass)
		list_of_speeds.extend([speed])
		vel_vect = np.array([speed*i for i in create_random_unit_vector()])
		particle['vel_vect'] = vel_vect

	return positions


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
				#LJ_potential = 4*epsilon_mix*((sigma_mix/r_mag)**12 - (sigma_mix/r_mag)**6)
				LF_force = 4*epsilon_mix*( -12*(sigma_mix**12/r_mag**13) + 6*(sigma_mix**6/r_mag**7))

				positions[i]['force_vect'] = positions[i]['force_vect'] + [LJ_force*val for val in r_vect/r_mag]
				if j <= len(positions)-1:
					positions[j]['force_vect'] = positions[j]['force_vect'] + [-1*LJ_potential*val for val in r_vect/r_mag]

	return positions


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
		pos_final = [0.5*(vel_final[i]+vel_init[i])*timestep + pos_init[i] for i in range(len(vel_init))]
		particle['vel_vect'] = np.array(vel_final)
		particle['pos_vect'] = np.array(pos_final)

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


def log_positions(positions, cycle_num, filename):
	if cycle_num == 0:
		csvfile = open(filename,'w', newline='')
	else:
		csvfile = open(filename,'a', newline='')
	writer = csv.writer(csvfile, delimiter="\t")
	writer.writerow(str(len(positions)))
	writer.writerow(['Cycle Number = '+str(cycle_num)])
	for particle in positions:
		type = particle['type']
		xyz = ' '.join(map(str, particle['pos_vect']))
		writer.writerow([type, xyz])


def run_MD(positions_file, parameters_file, log_filename, init_cycles, production_cycles, temperature, temp_range, box_size, r_cutoff, time_step):
	positions = load_positions(positions_file, box_size)
	parameters = load_parameters(parameters_file)
	initialize_MD(temperature, temp_range, positions, parameters)

	for i in range(init_cycles):
		cycle_num = i
		positions_with_replicas = replicate_cell(positions, box_size)
		calculate_forces(positions, positions_with_replicas, parameters, r_cutoff)
		integrate_forces(positions, parameters, timestep)
		check_boundaries(positions, box_size)
		log_positions(positions, cycle_num, log_filename)
	# for i in range(production_cycles):
	# 	cycle_num = init_cycles+i
	# 	positions_with_replicas = replicate_cell(positions, box_size)
	# 	calculate_forces(positions, positions_with_replicas, parameters, r_cutoff)
	# 	integrate_forces(positions, parameters, timestep)
	# 	check_boundaries(positions, box_size)
	# 	log_positions(positions, cycle_num, log_filename)


# Exeute Program
pos_file = 'test_pos.xyz'
params_file = 'LJ_params.def'
log_file = '/Users/brian_day/Desktop/log_file_text.xyz'
init_cycles = 100
prod_cycles = 200
temperature = 300 # Kelvin
temp_range = 30
box_size = [10,10,10] # Angstrom
r_cutoff = 8
timestep = 1e-3

run_MD(pos_file, params_file, log_file, init_cycles, prod_cycles, temperature, temp_range, box_size, r_cutoff, timestep
