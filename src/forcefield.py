import pandas as pd
import numpy as np

class Forcefield:
    def __init__(self) -> None:
        self.cutoff = 12
        self.ecutoff = 15

    def load_parameters(csvfile):
        """
        Load the Lennard-Jones (and other relevant parameters as updated) to be used in calculating the
        forces.
        """

        # Load the parameters file.
        output_data = pd.read_csv(csvfile)

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

class Bond:
    def __init__(self) -> None:
        pass

class Angle:
    def __init__(self) -> None:
        pass

class Dihedral:
    def __init__(self) -> None:
        pass

class Improper:
    def __init__(self) -> None:
        pass
