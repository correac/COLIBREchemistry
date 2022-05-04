import numpy as np
from .unitilies.helper_functions import cosmic_time_approx_Gyr, calculate_angmomentum
from swiftsimio.visualisation.rotation import rotation_matrix_from_vector


class GasParticleData:
    """
    General class containing gas particle properties
    """

    def __init__(self, sim_info, halo_id):

        mask_gas = sim_info.make_mask_gas(halo_id=halo_id)

        self.mass = sim_info.snapshot.gas.masses[mask_gas].value * sim_info.to_Msun_units
        self.n_parts = len(self.mass)

        self.xcoordinate = sim_info.snapshot.gas.coordinates[mask_gas, 0].value * sim_info.a * sim_info.to_kpc_units
        self.ycoordinate = sim_info.snapshot.gas.coordinates[mask_gas, 1].value * sim_info.a * sim_info.to_kpc_units
        self.zcoordinate = sim_info.snapshot.gas.coordinates[mask_gas, 2].value * sim_info.a * sim_info.to_kpc_units

        self.xvelocity = sim_info.snapshot.gas.velocities[mask_gas, 0].value  # km/s
        self.yvelocity = sim_info.snapshot.gas.velocities[mask_gas, 1].value  # km/s
        self.zvelocity = sim_info.snapshot.gas.velocities[mask_gas, 2].value  # km/s

        self.smoothing_length = sim_info.snapshot.gas.smoothing_lengths[mask_gas].value * sim_info.a * sim_info.to_kpc_units

        XH = sim_info.snapshot.gas.element_mass_fractions.hydrogen[mask_gas].value
        gas_HI = sim_info.snapshot.gas.species_fractions.HI[mask_gas].value
        gas_H2 = sim_info.snapshot.gas.species_fractions.H2[mask_gas].value * 2.0

        self.HI_mass = gas_HI * XH * self.mass
        self.H2_mass = gas_H2 * XH * self.mass

        self.star_formation_rates = sim_info.snapshot.gas.star_formation_rates[mask_gas].value * sim_info.to_Msun_units / sim_info.to_yr_units
        self.densities = sim_info.snapshot.gas.densities[mask_gas].value * (sim_info.a * sim_info.to_Msun_units / sim_info.to_kpc_units) ** 3
        self.metal_mass_fractions = sim_info.snapshot.gas.metal_mass_fractions[mask_gas].value / sim_info.Zsolar


class StarParticleData:
    """
    General class containing gas particle properties
    """

    def __init__(self, sim_info, halo_id):
        mask_stars = sim_info.make_mask_stars(halo_id=halo_id)

        self.mass = sim_info.snapshot.stars.masses[mask_stars].value * sim_info.to_Msun_units
        stars_birthz = (
                1.0 / sim_info.snapshot.stars.birth_scale_factors[mask_stars].value - 1.0
        )

        if len(stars_birthz) > 1:
            self.age = cosmic_time_approx_Gyr(
                z=0.0, Omega_L=sim_info.Omega_l, Hubble_time=sim_info.hubble_time_Gyr
            ) - cosmic_time_approx_Gyr(
                z=stars_birthz, Omega_L=sim_info.Omega_l, Hubble_time=sim_info.hubble_time_Gyr
            )
        else:
            self.age = 0.0

        self.metal_mass_fractions = sim_info.snapshot.stars.metal_mass_fractions[mask_stars].value
        self.initmass = sim_info.snapshot.stars.initial_masses[mask_stars].value * sim_info.to_Msun_units

        self.n_parts = len(self.mass)

        self.xcoordinate = sim_info.snapshot.stars.coordinates[mask_stars, 0].value * sim_info.a * sim_info.to_kpc_units
        self.ycoordinate = sim_info.snapshot.stars.coordinates[mask_stars, 1].value * sim_info.a * sim_info.to_kpc_units
        self.zcoordinate = sim_info.snapshot.stars.coordinates[mask_stars, 2].value * sim_info.a * sim_info.to_kpc_units

        self.xvelocity = sim_info.snapshot.stars.velocities[mask_stars, 0].value  # km/s
        self.yvelocity = sim_info.snapshot.stars.velocities[mask_stars, 1].value  # km/s
        self.zvelocity = sim_info.snapshot.stars.velocities[mask_stars, 2].value  # km/s

        self.baryon_max_soft = 0.5 * sim_info.baryon_max_soft * np.ones(self.n_parts)
        self.weighted_mass = self.mass * (1.2348 / self.baryon_max_soft) ** 3

        self.oxygen = sim_info.snapshot.stars.element_mass_fractions.oxygen[mask_stars].value
        self.iron = sim_info.snapshot.stars.element_mass_fractions.iron[mask_stars].value
        self.magnesium = sim_info.snapshot.stars.element_mass_fractions.magnesium[mask_stars].value
        self.hydrogen = sim_info.snapshot.stars.element_mass_fractions.hydrogen[mask_stars].value
        self.carbon = sim_info.snapshot.stars.element_mass_fractions.carbon[mask_stars].value
        self.silicon = sim_info.snapshot.stars.element_mass_fractions.silicon[mask_stars].value
        self.europium = sim_info.snapshot.stars.element_mass_fractions.europium[mask_stars].value
        self.nitrogen = sim_info.snapshot.stars.element_mass_fractions.nitrogen[mask_stars].value
        self.neon = sim_info.snapshot.stars.element_mass_fractions.neon[mask_stars].value
        self.iron_SNIa_fraction = sim_info.snapshot.stars.iron_mass_fractions_from_snia[mask_stars].value

        if (hasattr(sim_info.snapshot.stars.element_mass_fractions, 'barium')):
            self.barium = sim_info.snapshot.stars.element_mass_fractions.barium[mask_stars].value
            self.strontium = sim_info.snapshot.stars.element_mass_fractions.strontium[mask_stars].value
        else:
            self.barium = np.zeros(self.n_parts)
            self.strontium = np.zeros(self.n_parts)

        self.in_halo = np.zeros(self.n_parts)

        indx = sim_info.halo_data.halo_ids == halo_id

        particles_pos = np.zeros((self.n_parts,3))
        particles_pos[:,0] = self.xcoordinate - sim_info.halo_data.xminpot[indx]
        particles_pos[:,1] = self.ycoordinate - sim_info.halo_data.yminpot[indx]
        particles_pos[:,2] = self.zcoordinate - sim_info.halo_data.zminpot[indx]

        particles_vel = np.zeros((self.n_parts,3))
        particles_vel[:,0] = self.xvelocity - sim_info.halo_data.vxminpot[indx]
        particles_vel[:,1] = self.yvelocity - sim_info.halo_data.vyminpot[indx]
        particles_vel[:,2] = self.zvelocity - sim_info.halo_data.vzminpot[indx]

        momentum, momentum_inner_5kpc, distance = \
            calculate_angmomentum(sim_info, particles_pos, particles_vel, self.mass )

        self.distance = distance

        face_on_rotation_matrix = rotation_matrix_from_vector(momentum, axis="z")

        particles_pos_face_on = np.matmul(face_on_rotation_matrix, particles_pos.T)
        particles_pos_face_on = particles_pos_face_on.T

        self.R = np.sqrt( particles_pos_face_on[:,0]**2 + particles_pos_face_on[:,1]**2)
        self.z = particles_pos_face_on[:, 2]

        #self.luminosity = sim_info.snapshot.stars.luminosity.rband[mask_stars].value
