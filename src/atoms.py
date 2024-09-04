import numpy as np

import csv
from typing import TypeAlias, NewType

from forcefield import Forcefield


class Atom:
    """_summary_

    Attributes
    ----------
    id: int
    molecule_id: int
    type: str
    mass: float
    charge: float
    position: np.ndarray(float)
    velocty: np.ndarray(float)
    acceleration: np.ndarray(float)
    force: np.ndarray(float)
    """

    def __init__(self) -> None:
        self.id
        self.molecule_id
        self.type
        self.mass
        self.charge
        self.position
        self.velocity
        self.acceleration
        self.force


class Atoms:
    """Container class for Atoms.

    Attributes
    ----------
    atoms: np.ndarray of Atom
    bonds: np.ndarray of tuple(int, int)
    angles: np.ndarray of tuple(int, int, int)
    impropers: np.ndarray of tuple(int, int, int, int)
    dihedrals: np.ndarray of tuple(int, int, int, int)
    nearest_neighbors: np.ndarray of list[int]
        jagged array of the atom ids neighboring the atom of index i
    forcefield: Forcefield | None
    potential_energy: float | None
    kinetic_energy: float | None

    Notes
    ------
    Add methods for bulk/avg atomic properties, like mass, charge, velocity, acceleration, center of mass, etc.

    """
    def __init__(self, molecules) -> None:
        self.atoms = np.Array()
        self.bonds = np.Array()
        self.angles = np.Array()
        self.dihedrals = np.Array()
        self.impropers = np.Array()
        self.nearest_neighbors = np.Array()

        self.forcefield = None
        self.potential_energy = None
        self.kinetic_energy = None

    def load_positions(filename):
        """
        Load the particle types and positions as a dictionary by atom number.
        """

        # Load the positions file.
        with open(filename, newline='') as csvfile:
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

    def apply_pbc(self):
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

    def calculate_nearest_neighbors():
        pass

    def validate_forcefield_for_atoms():
        pass
