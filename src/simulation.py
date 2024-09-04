import numpy as np

import csv
import logging
import argparse

from atoms import Atom, Atoms
from forcefield import Forcefield


class MolecularDynamics:
    """_summary_

    Attributes
    ----------

    Notes
    -----
    - Want to keep the notion of atoms and forcefields as separate as possible, so we do not make potential or kinetic
    energy properties of the Atoms class, but rather here where we first combine Atoms with Forcefields. Could
    reconsider making forcefield a property of the atoms themselves, and MD class a "what-to-do-with-that" simulator.
    Probably a good idea since MC could use same atoms/forcefields but do something different. 
    """
    def __init__(self, **kwargs):
        default_kwargs = {
            'timestep': 0.1,
            'start_time': 0,
            'stop_time': 100
        }
        kwargs = default_kwargs | kwargs

        self.system = Atoms()
        self.start_time = kwargs['start_time']
        self.stop_time = kwargs['stop_time']
        self.dt = kwargs['timestep']

        # TODO: make this optional, also implement it lol
        self.nearest_neighbor_recalculation = 10

        self.logger = logging.getLogger(__name__)
        logging.basicConfig(filename='example.log', encoding='utf-8', level=logging.DEBUG)

    def load_input_file(file):
        pass

    def update_kwargs():
        pass

    # TODO: make utility function in another file
    @classmethod
    def create_random_unit_vector():
        """
        Creates a random unit vector.
        """
        vector = np.array([rand.uniform(-1,1), rand.uniform(-1,1), rand.uniform(-1,1)])
        vector_mag = np.linalg.norm(vector)
        unit_vector = np.array([val/vector_mag for val in vector])
        return unit_vector
    
    def initialize_MD(self, temperature, positions, parameters, distribution='Boltzmann'):
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
        
    def integrate_forces(self):
        """
        Uses the calculate forces and the time step to update the positions of the particles in the
        system.
        """
        positions = self.positions
        parameters = self.forcefield
        timestep = self.timestep
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

    def integrate_forces_verlet(self, box_size, r_cutoff):
        positions = self.positions
        parameters = self.forcefield
        timestep = self.timestep
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

        self.cell.apply_pbc()

        for particle in positions:

            type = particle['type']
            mass = parameters[type]['mass']
            vel_half = particle['vel_half']

            accel_final = [force/mass for force in particle['force_vect']]
            vel_final = [vel_half[i] + accel_final[i]*timestep/2 for i in range(len(vel_half))]
            particle['vel_vect'] = np.array(vel_final)
            particle['KE'] = 0.5*mass*np.linalg.norm(particle['vel_vect'])**2

        return positions

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

    def run(self, positions_file, parameters_file, log_filename, temperature, box_size, r_cutoff, time_step, log_freq, log_prop_filename):
        positions_with_replicas = replicate_cell(positions, box_size)
        PE_total = calculate_forces(positions, positions_with_replicas, parameters, r_cutoff)

        log_properties(positions, -1, PE_total, log_prop_filename)

        times = np.linspace(self.start_time, self.stop_time, self.timestep)
        for i, t in enumerate(times):
            integrate_forces_verlet(positions, parameters, timestep, box_size, r_cutoff)
            self.cell.apply_pbc()
            positions_with_replicas = replicate_cell(positions, box_size)
            PE_total = calculate_forces(positions, positions_with_replicas, parameters, r_cutoff)
            if i % log_freq == 0:
                log_properties(self.positions, i, PE_total, log_prop_filename)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'Say hello')
    parser.add_argument('name', help='your name, enter it')
    args = parser.parse_args()
    main(args.name)