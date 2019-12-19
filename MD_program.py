# Example Molecular Dynamics Code
# For understanding MD / Practice with Python, Github, etc.
# Written by: Brian brian_day

# Import necessary packages
import csv
import numpy as np
from random import random

# Define functions
def load_positions(filename):
	"""
	Load the particle types and positions as a dictionary by atom number.
	"""
	with open(filename,newline='') as csvfile:
		output_data = csv.reader(csvfile, delimiter="\t")
		output_data = list(output_data)

		positions_list = []
		for i in range(1,len(output_data)):
			row = output_data[i]
			temp_dict = {}
			temp_dict['num'] = i
			temp_dict['type'] = row[0]
			temp_dict['pos_vect'] = np.array([float(item) for item in row[1:4]])
			temp_dict['vel_vect'] = ''
			positions_list.append(temp_dict)
		return positions_list


def load_parameters(filename):
	"""
	Load the Lennard-Jones (and other relevant parameters as updated) to be used in calculating the
	forces.
	"""
	with open(filename,newline='') as csvfile:
		output_data = csv.DictReader(csvfile, delimiter="\t")
		return list(output_data)


def initialize_MD(temperature, positions, parameters):
	"""
	Takes in the temperture to define a set of inital velocity vectors.
	"""
	for particle in positions:
		speed = 5
		vel_vect = speed* np.array([random(), random(), random()])
		# Need to fix the above def. Will affect velocity since not a unit vector.
		particle['vel_vect'] = vel_vect
	return positions


def calculate_forces():
    """
    Evaluates Lennard-Jones potentials for all particles and assigns forces.
    Update later to include charge-charge interactions.
    """
    return


def integrate_forces():
    """
    Uses the calculate forces and the time step to update the positions of the particles in the
    system.
    """
    return


def check_boundaries():
    """
    Check if any particle has gone beyond the boundary, and if so update position through periodic
    boundary conditions.
    """
	for particle in positions:
		# Check x-positiions
		if particle['pos_vect'][0] < 0:
			['pos_vect'][0] = box_x + ['pos_vect'][0];
		if particle['pos_vect'][0] > box_x:
			['pos_vect'][0] = ['pos_vect'][0] - box_x;

		# Check x-positiions
		if particle['pos_vect'][1] < 0:
			['pos_vect'][1] = box_y + ['pos_vect'][1];
		if particle['pos_vect'][1] > box_y:
			['pos_vect'][1] = ['pos_vect'][1] - box_y;

		# Check x-positiions
		if particle['pos_vect'][2] < 0:
			['pos_vect'][2] = box_z + ['pos_vect'][2];
		if particle['pos_vect'][2] > box_z:
			['pos_vect'][2] = ['pos_vect'][2] - box_z;

    return


def run_MD(init_cycles, prod_cycles, positions, parameters, temperature, r_cutoff):
    particles = load_positions
    parameters = load_parameters
    initialize_MD()

    for i in range(init_cycles):
        forces = calculate_forces()
        pos_vects_temp, vel_vects = integrate_forces()
        pos_vects = check_boundaries()
    for i in range(production_cycles):
        forces = calculate_forces()
        pos_vects_temp, vel_vects = integrate_forces()
        pos_vects = check_boundaries()

    return


# Exeute Program
params_file = 'LJ_params.def'
pos_file = 'test_pos.xyz'
init_cycles = 100
prod_cycles = 200
r_cutoff = 10
temperature = 300

var = load_positions(pos_file)
var = initialize_MD(5, var, 5)
