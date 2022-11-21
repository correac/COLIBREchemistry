import numpy as np

def calculate_kappa_co(sim_info, halo_index):
    # subhalo contain subhalo data and is strutured as follow
    # [ (0:3)CentreOfPotential[kpc]: (0)X | (1)Y | (2)Z  | (3:6)Velocity[km/s]: (3)Vx | (4)Vy | (5)Vz  | (6)R200c[kpc]]

    # let's create super array (partsDATA) that contains particles data and is structured as follow
    # [ (:3)Position[kpc]: (0)X | (1)Y | (2)Z  | (3)Mass[Msun] | (4:7)Velocity[km/s]: (4)Vx | (5)Vy | (6)Vz | (7)hsml]

    len_partsDATA = sim_info.stars.n_parts

    particlesDATA = np.zeros( (len_partsDATA, 7) )
    particlesDATA[:, 0] = sim_info.stars.xcoordinate.copy()
    particlesDATA[:, 1] = sim_info.stars.ycoordinate.copy()
    particlesDATA[:, 2] = sim_info.stars.zcoordinate.copy()
    particlesDATA[:, 3] = sim_info.stars.mass.copy()

    particlesDATA[:, 4] = sim_info.stars.xvelocity.copy()
    particlesDATA[:, 5] = sim_info.stars.yvelocity.copy()
    particlesDATA[:, 6] = sim_info.stars.zvelocity.copy()


    # Centering onto subhalo CoP
    particlesDATA[:, 0] -= sim_info.halo_data.xminpot[halo_index]
    particlesDATA[:, 1] -= sim_info.halo_data.yminpot[halo_index]
    particlesDATA[:, 2] -= sim_info.halo_data.zminpot[halo_index]
    particlesDATA[:, :3] += sim_info.boxSize / 2
    particlesDATA[:, :3] %= sim_info.boxSize
    particlesDATA[:, :3] -= sim_info.boxSize / 2  # end the unwrap

    # Center velocities on the subhalo CoM velocity
    particlesDATA[:, 4] -= sim_info.halo_data.vxminpot[halo_index]
    particlesDATA[:, 5] -= sim_info.halo_data.vyminpot[halo_index]
    particlesDATA[:, 6] -= sim_info.halo_data.vzminpot[halo_index]

    # Compute distances
    distancesDATA = np.linalg.norm(particlesDATA[:, :3], axis=1)

    # Restrict particles
    extract = distancesDATA < 30.0

    particlesDATA = particlesDATA[extract, :]
    distancesDATA = distancesDATA[extract]

    Mstar = np.sum(particlesDATA[:, 3])  # compute total in-aperture stellar mass
    # Compute 30kpc CoM to Sub CoM velocty offset & recenter
    dvVmass = (
        np.sum(particlesDATA[:, 3][:, np.newaxis] * particlesDATA[:, 4:7], axis=0)
        / Mstar
    )
    particlesDATA[:, 4:7] -= dvVmass

    # Compute momentum
    smomentums = np.cross(particlesDATA[:, :3], particlesDATA[:, 4:7])
    momentum = np.sum(particlesDATA[:, 3][:, np.newaxis] * smomentums, axis=0)

    # Compute rotational velocities
    smomentumz = np.sum(momentum * smomentums / np.linalg.norm(momentum), axis=1)
    cyldistances = (
        distancesDATA ** 2
        - np.sum(momentum * particlesDATA[:, :3] / np.linalg.norm(momentum), axis=1)
        ** 2
    )
    cyldistances = np.sqrt(np.abs(cyldistances))

    if len(cyldistances[cyldistances > 0]) > 0:
        cylmin = np.min(cyldistances[cyldistances > 0])
        cyldistances[cyldistances == 0] = cylmin
        vrots = smomentumz / cyldistances
    else:
        vrots = smomentumz

    # Compute kappa_co
    Mvrot2 = np.sum((particlesDATA[:, 3] * vrots ** 2)[vrots > 0])

    kappa_co = Mvrot2 / np.sum(
        particlesDATA[:, 3] * (np.linalg.norm(particlesDATA[:, 4:7], axis=1)) ** 2
    )

    # Return
    return kappa_co