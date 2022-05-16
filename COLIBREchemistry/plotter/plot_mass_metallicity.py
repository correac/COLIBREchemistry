import matplotlib.pylab as plt
from matplotlib.pylab import rcParams
import numpy as np
from tqdm import tqdm
import scipy.stats as stat


# def compute_median_ratios(Oxygen_fraction, Iron_fraction,
#                    Magnesium_fraction,Hydrogen_fraction):
#
#     mp_in_cgs = 1.6737236e-24
#     mH_in_cgs = 1.00784 * mp_in_cgs
#     mFe_in_cgs = 55.845 * mp_in_cgs
#     mO_in_cgs = 15.999 * mp_in_cgs
#     mMg_in_cgs = 24.305 * mp_in_cgs
#
#     # Asplund et al. (2009)
#     Fe_H_Sun = 7.5
#     O_H_Sun = 8.69
#     Mg_H_Sun = 7.6
#
#     O_Fe_Sun = O_H_Sun - Fe_H_Sun - np.log10(mFe_in_cgs / mO_in_cgs)
#     Mg_Fe_Sun = Mg_H_Sun - Fe_H_Sun - np.log10(mFe_in_cgs / mMg_in_cgs)
#     Fe_H_Sun = Fe_H_Sun - 12.0 - np.log10(mH_in_cgs / mFe_in_cgs)
#
#     Fe_H = np.log10(Iron_fraction / Hydrogen_fraction) - Fe_H_Sun
#     O_Fe = np.log10(Oxygen_fraction / Iron_fraction) - O_Fe_Sun
#     Mg_Fe = np.log10(Magnesium_fraction / Iron_fraction) - Mg_Fe_Sun
#
#     # Let's set lower and upper limits:
#     Fe_H[Iron_fraction == 0] = -7  # set lower limit
#     Fe_H[Fe_H < -7] = -7  # set lower limit
#     Mg_Fe[Iron_fraction == 0] = -2  # set lower limit
#     Mg_Fe[Magnesium_fraction == 0] = -2  # set lower limit
#     Mg_Fe[Mg_Fe < -2] = -2  # set lower limit
#     O_Fe[Iron_fraction == 0] = -2  # set lower limit
#     O_Fe[Oxygen_fraction == 0] = -2  # set lower limit
#     O_Fe[O_Fe < -2] = -2  # set lower limit
#
#
#     select = np.where(Fe_H >= -3)[0]
#     Fe_H = Fe_H[select]
#     O_Fe = O_Fe[select]
#     Mg_Fe = Mg_Fe[select]
#
#     Fe_H = np.median(Fe_H)
#     O_Fe = np.median(O_Fe)
#     Mg_Fe = np.median(Mg_Fe)
#
#     return {'Fe_H': Fe_H, 'O_Fe': O_Fe, 'Mg_Fe': Mg_Fe}

# def compute_mass_weighted_ratios(Oxygen_fraction, Iron_fraction, Magnesium_fraction,
#                                  Hydrogen_fraction, stellar_mass):
#
#     mp_in_cgs = 1.6737236e-24
#     mH_in_cgs = 1.00784 * mp_in_cgs
#     mFe_in_cgs = 55.845 * mp_in_cgs
#     mO_in_cgs = 15.999 * mp_in_cgs
#     mMg_in_cgs = 24.305 * mp_in_cgs
#
#     # Asplund et al. (2009)
#     Fe_H_Sun = 7.5
#     O_H_Sun = 8.69
#     Mg_H_Sun = 7.6
#
#     O_Fe_Sun = O_H_Sun - Fe_H_Sun - np.log10(mFe_in_cgs / mO_in_cgs)
#     Mg_Fe_Sun = Mg_H_Sun - Fe_H_Sun - np.log10(mFe_in_cgs / mMg_in_cgs)
#     Fe_H_Sun = Fe_H_Sun - 12.0 - np.log10(mH_in_cgs / mFe_in_cgs)
#
#     Fe_H = Iron_fraction / Hydrogen_fraction
#     O_Fe = Oxygen_fraction / Iron_fraction
#     Mg_Fe = Magnesium_fraction / Iron_fraction
#
#     # Let's set lower and upper limits:
#     Mg_Fe[Iron_fraction == 0] = 0.  # set lower limit
#     O_Fe[Iron_fraction == 0] = 0.   # set lower limit
#
#     Fe_H_mass_weighted = np.sum( stellar_mass * Iron_fraction) / np.sum( Hydrogen_fraction * stellar_mass)
#     O_Fe_mass_weighted = np.sum( stellar_mass * Oxygen_fraction) / np.sum( Iron_fraction * stellar_mass)
#     Mg_Fe_mass_weighted = np.sum( stellar_mass * Magnesium_fraction) / np.sum( Iron_fraction * stellar_mass)
#     Fe_H_mass_weighted = np.log10(Fe_H_mass_weighted) - Fe_H_Sun
#     O_Fe_mass_weighted = np.log10(O_Fe_mass_weighted) - O_Fe_Sun
#     Mg_Fe_mass_weighted = np.log10(Mg_Fe_mass_weighted) - Mg_Fe_Sun
#
#     Fe_H_mass_weighted_ratio = np.sum( stellar_mass * Fe_H) / np.sum( stellar_mass)
#     O_Fe_mass_weighted_ratio = np.sum( stellar_mass * O_Fe) / np.sum( stellar_mass)
#     Mg_Fe_mass_weighted_ratio = np.sum( stellar_mass * Mg_Fe) / np.sum( stellar_mass)
#     Fe_H_mass_weighted_ratio = np.log10(Fe_H_mass_weighted_ratio) - Fe_H_Sun
#     O_Fe_mass_weighted_ratio = np.log10(O_Fe_mass_weighted_ratio) - O_Fe_Sun
#     Mg_Fe_mass_weighted_ratio = np.log10(Mg_Fe_mass_weighted_ratio) - Mg_Fe_Sun
#
#     return {'Fe_H': Fe_H_mass_weighted,
#             'O_Fe': O_Fe_mass_weighted,
#             'Mg_Fe': Mg_Fe_mass_weighted,
#             'Fe_H_weighted_ratio': Fe_H_mass_weighted_ratio,
#             'O_Fe_weighted_ratio': O_Fe_mass_weighted_ratio,
#             'Mg_Fe_weighted_ratio': Mg_Fe_mass_weighted_ratio}
#
# def compute_light_weighted_ratios(Z, Fe, H, mass, L_r_band, L_i_band, L_Z_band):
#
#     mp_in_cgs = 1.6737236e-24
#     mH_in_cgs = 1.00784 * mp_in_cgs
#     mFe_in_cgs = 55.845 * mp_in_cgs
#
#     # Asplund et al. (2009)
#     Fe_H_Sun = 7.5
#     Z_Sun = 0.0134
#     Fe_H_Sun = Fe_H_Sun - 12.0 - np.log10(mH_in_cgs / mFe_in_cgs)
#
#     # Ratio
#     Fe_H = Fe / H
#
#     Z_weighted_r_band = np.sum(Z * L_r_band) / np.sum(L_r_band)
#     Z_weighted_r_band /= Z_Sun
#
#     Z_weighted_i_band = np.sum(Z * L_i_band) / np.sum(L_i_band)
#     Z_weighted_i_band /= Z_Sun
#
#     Z_weighted_Z_band = np.sum(Z * L_Z_band) / np.sum(L_Z_band)
#     Z_weighted_Z_band /= Z_Sun
#
#     Fe_H_weighted_r_band = np.sum(Fe_H * L_r_band) / np.sum(L_r_band)
#     Fe_H_weighted_i_band = np.sum(Fe_H * L_i_band) / np.sum(L_i_band)
#     Fe_H_weighted_Z_band = np.sum(Fe_H * L_Z_band) / np.sum(L_Z_band)
#     Fe_H_weighted_r_band = np.log10(Fe_H_weighted_r_band) - Fe_H_Sun
#     Fe_H_weighted_i_band = np.log10(Fe_H_weighted_i_band) - Fe_H_Sun
#     Fe_H_weighted_Z_band = np.log10(Fe_H_weighted_Z_band) - Fe_H_Sun
#
#
#     Fe_weighted = np.sum(Fe * mass * L_r_band) / np.sum(L_r_band)
#     H_weighted = np.sum(H * mass * L_r_band) / np.sum(L_r_band)
#     Fe_H_weighted_r_band_with_mass = np.log10(Fe_weighted / H_weighted) - Fe_H_Sun
#     Fe_weighted = np.sum(Fe * mass * L_i_band) / np.sum(L_i_band)
#     H_weighted = np.sum(H * mass * L_i_band) / np.sum(L_i_band)
#     Fe_H_weighted_i_band_with_mass = np.log10(Fe_weighted / H_weighted) - Fe_H_Sun
#     Fe_weighted = np.sum(Fe * mass * L_Z_band) / np.sum(L_Z_band)
#     H_weighted = np.sum(H * mass * L_Z_band) / np.sum(L_Z_band)
#     Fe_H_weighted_Z_band_with_mass = np.log10(Fe_weighted / H_weighted) - Fe_H_Sun
#
#     return {'Fe_H_weighted_r_band': Fe_H_weighted_r_band,
#             'Fe_H_weighted_i_band': Fe_H_weighted_i_band,
#             'Fe_H_weighted_Z_band': Fe_H_weighted_Z_band,
#             'Z_weighted_r_band': Z_weighted_r_band,
#             'Z_weighted_i_band': Z_weighted_i_band,
#             'Z_weighted_Z_band': Z_weighted_Z_band,
#             'Fe_H_weighted_r_band_with_mass': Fe_H_weighted_r_band_with_mass,
#             'Fe_H_weighted_i_band_with_mass': Fe_H_weighted_i_band_with_mass,
#             'Fe_H_weighted_Z_band_with_mass': Fe_H_weighted_Z_band_with_mass}

# def compute_alpha_light_weighted_ratios(Fe, O, Mg, mass, L_r_band):
#
#     mp_in_cgs = 1.6737236e-24
#     mH_in_cgs = 1.00784 * mp_in_cgs
#     mFe_in_cgs = 55.845 * mp_in_cgs
#     mO_in_cgs = 15.999 * mp_in_cgs
#     mMg_in_cgs = 24.305 * mp_in_cgs
#
#     # Asplund et al. (2009)
#     Fe_H_Sun = 7.5
#     O_H_Sun = 8.69
#     Mg_H_Sun = 7.6
#
#     O_Fe_Sun = O_H_Sun - Fe_H_Sun - np.log10(mFe_in_cgs / mO_in_cgs)
#     Mg_Fe_Sun = Mg_H_Sun - Fe_H_Sun - np.log10(mFe_in_cgs / mMg_in_cgs)
#
#     # Ratios
#     O_Fe = O / Fe
#     Mg_Fe = Mg / Fe
#     Mg_Fe[Fe == 0] = 0.  # set lower limit
#     O_Fe[Fe == 0] = 0.  # set lower limit
#
#     O_Fe_weighted_r_band = np.sum(O_Fe * L_r_band) / np.sum(L_r_band)
#     Mg_Fe_weighted_r_band = np.sum(Mg_Fe * L_r_band) / np.sum(L_r_band)
#     O_Fe_weighted_r_band = np.log10(O_Fe_weighted_r_band) - O_Fe_Sun
#     Mg_Fe_weighted_r_band = np.log10(Mg_Fe_weighted_r_band) - Mg_Fe_Sun
#
#     Fe_weighted = np.sum(Fe * mass * L_r_band) / np.sum(L_r_band)
#     O_weighted = np.sum(O * mass * L_r_band) / np.sum(L_r_band)
#     O_Fe_weighted_r_band_with_mass = np.log10(O_weighted / Fe_weighted) - O_Fe_Sun
#     Mg_weighted = np.sum(Mg * mass * L_r_band) / np.sum(L_r_band)
#     Mg_Fe_weighted_r_band_with_mass = np.log10(Mg_weighted / Fe_weighted) - Mg_Fe_Sun
#
#     return {'O_Fe_weighted_r_band': O_Fe_weighted_r_band,
#             'Mg_Fe_weighted_r_band': Mg_Fe_weighted_r_band,
#             'O_Fe_weighted_r_band_with_mass': O_Fe_weighted_r_band_with_mass,
#             'Mg_Fe_weighted_r_band_with_mass': Mg_Fe_weighted_r_band_with_mass}

# def compute_metallicity_weighted_ratios(Z, Fe, H, O, Mg, L_r_band):
#
#     mp_in_cgs = 1.6737236e-24
#     mH_in_cgs = 1.00784 * mp_in_cgs
#     mFe_in_cgs = 55.845 * mp_in_cgs
#     mO_in_cgs = 15.999 * mp_in_cgs
#     mMg_in_cgs = 24.305 * mp_in_cgs
#
#     # Asplund et al. (2009)
#     Fe_H_Sun = 7.5
#     O_H_Sun = 8.69
#     Mg_H_Sun = 7.6
#
#     O_Fe_Sun = O_H_Sun - Fe_H_Sun - np.log10(mFe_in_cgs / mO_in_cgs)
#     Mg_Fe_Sun = Mg_H_Sun - Fe_H_Sun - np.log10(mFe_in_cgs / mMg_in_cgs)
#     Fe_H_Sun = Fe_H_Sun - 12.0 - np.log10(mH_in_cgs / mFe_in_cgs)
#
#     # Ratios
#     Fe_H = Fe / H
#     O_Fe = O / Fe
#     Mg_Fe = Mg / Fe
#     Mg_Fe[Fe == 0] = 0.  # set lower limit
#     O_Fe[Fe == 0] = 0.  # set lower limit
#
#     Fe_H_Z_weighted_r_band = np.sum(Fe_H * Z * L_r_band) / np.sum(Z * L_r_band)
#     Fe_H_Z_weighted_r_band = np.log10(Fe_H_Z_weighted_r_band) - Fe_H_Sun
#     Fe_H_Z_weighted = np.sum(Fe_H * Z) / np.sum(Z)
#     Fe_H_Z_weighted = np.log10(Fe_H_Z_weighted) - Fe_H_Sun
#
#     O_Fe_Z_weighted_r_band = np.sum(O_Fe * Z * L_r_band) / np.sum(Z * L_r_band)
#     O_Fe_Z_weighted_r_band = np.log10(O_Fe_Z_weighted_r_band) - O_Fe_Sun
#     O_Fe_Z_weighted = np.sum(O_Fe * Z) / np.sum(Z)
#     O_Fe_Z_weighted = np.log10(O_Fe_Z_weighted) - O_Fe_Sun
#
#     Mg_Fe_Z_weighted_r_band = np.sum(Mg_Fe * Z * L_r_band) / np.sum(Z * L_r_band)
#     Mg_Fe_Z_weighted_r_band = np.log10(Mg_Fe_Z_weighted_r_band) - Mg_Fe_Sun
#     Mg_Fe_Z_weighted = np.sum(Mg_Fe * Z) / np.sum(Z)
#     Mg_Fe_Z_weighted = np.log10(Mg_Fe_Z_weighted) - Mg_Fe_Sun
#
#     return {'Fe_H_Z_weighted_r_band': Fe_H_Z_weighted_r_band,
#             'Fe_H_Z_weighted': Fe_H_Z_weighted,
#             'O_Fe_Z_weighted_r_band': O_Fe_Z_weighted_r_band,
#             'O_Fe_Z_weighted': O_Fe_Z_weighted,
#             'Mg_Fe_Z_weighted_r_band': Mg_Fe_Z_weighted_r_band,
#             'Mg_Fe_Z_weighted': Mg_Fe_Z_weighted}


def compute_weighted_FeH(X_Fe, X_H, weight):

    mp_in_cgs = 1.6737236e-24
    mH_in_cgs = 1.00784 * mp_in_cgs
    mFe_in_cgs = 55.845 * mp_in_cgs

    # Asplund et al. (2009)
    Fe_H_Sun = 7.5
    Fe_H_Sun = Fe_H_Sun - 12.0 - np.log10(mH_in_cgs / mFe_in_cgs)

    log_of_weight = np.sum ( X_Fe * weight ) / np.sum ( X_H * weight )
    log_of_weight = np.log10( log_of_weight ) - Fe_H_Sun
    log_of_weighted_ratio = np.sum ( (X_Fe/X_H) * weight ) / np.sum( weight )
    log_of_weighted_ratio = np.log10(log_of_weighted_ratio) - Fe_H_Sun
    weight_of_log = np.log10( X_Fe / X_H)

    # Let's not include Fe=0 stars, since that introduces too low values for weight of log,
    # Also, setting an arbitrary lower limit complicates comparison.
    nonzero = X_Fe != 0
    weight_of_log = weight_of_log[nonzero]
    weight = weight[nonzero]
    weight_of_log = np.sum( weight_of_log * weight) / np.sum(weight) - Fe_H_Sun

    return {'log_of_weight': log_of_weight,
            'log_of_weighted_ratio': log_of_weighted_ratio,
            'weight_of_log': weight_of_log}

def compute_weighted_OFe(X_O, X_Fe, weight):

    mp_in_cgs = 1.6737236e-24
    mFe_in_cgs = 55.845 * mp_in_cgs
    mO_in_cgs = 15.999 * mp_in_cgs

    # Asplund et al. (2009)
    Fe_H_Sun = 7.5
    O_H_Sun = 8.69
    O_Fe_Sun = O_H_Sun - Fe_H_Sun - np.log10(mFe_in_cgs / mO_in_cgs)

    ratio = X_O/X_Fe
    ratio[X_Fe==0] = 0

    log_of_weight = np.sum ( X_O * weight ) / np.sum ( X_Fe * weight )
    log_of_weight = np.log10( log_of_weight ) - O_Fe_Sun

    log_of_weighted_ratio = np.sum ( ratio * weight ) / np.sum( weight )
    log_of_weighted_ratio = np.log10(log_of_weighted_ratio) - O_Fe_Sun

    weight_of_log = np.log10( X_O / X_Fe)
    weight_of_log[X_O==0] = -2
    weight_of_log[X_Fe==0] = -2
    weight_of_log = np.sum( weight_of_log * weight) / np.sum(weight) - O_Fe_Sun

    return {'log_of_weight': log_of_weight,
            'log_of_weighted_ratio': log_of_weighted_ratio,
            'weight_of_log': weight_of_log}


def compute_weighted_MgFe(X_Mg, X_Fe, weight):

    mp_in_cgs = 1.6737236e-24
    mFe_in_cgs = 55.845 * mp_in_cgs
    mMg_in_cgs = 24.305 * mp_in_cgs

    # Asplund et al. (2009)
    Fe_H_Sun = 7.5
    Mg_H_Sun = 7.6

    Mg_Fe_Sun = Mg_H_Sun - Fe_H_Sun - np.log10(mFe_in_cgs / mMg_in_cgs)

    ratio = X_Mg / X_Fe
    ratio[X_Fe == 0] = 0

    log_of_weight = np.sum(X_Mg * weight) / np.sum(X_Fe * weight)
    log_of_weight = np.log10(log_of_weight) - Mg_Fe_Sun

    log_of_weighted_ratio = np.sum(ratio * weight) / np.sum(weight)
    log_of_weighted_ratio = np.log10(log_of_weighted_ratio) - Mg_Fe_Sun

    weight_of_log = np.log10(X_Mg / X_Fe)
    weight_of_log[X_Mg == 0] = -2
    weight_of_log[X_Fe == 0] = -2
    weight_of_log = np.sum(weight_of_log * weight) / np.sum(weight) - Mg_Fe_Sun

    return {'log_of_weight': log_of_weight,
            'log_of_weighted_ratio': log_of_weighted_ratio,
            'weight_of_log': weight_of_log}


def calculate_galaxies_metallicity(sim_info):

    Zsun = 0.0134

    select_sample = np.where((sim_info.halo_data.log10_stellar_mass >= 7.0) &
                              (sim_info.halo_data.log10_stellar_mass <= 12))[0]

    select_centrals = np.where(sim_info.halo_data.type[select_sample] == 10)[0]

    sample = select_sample[select_centrals]
    num_sample = len(sample)

    Z_mass_weighted = np.zeros(num_sample)
    Z_light_weighted = np.zeros(num_sample)
    #Z_light_weighted_i_band = np.zeros(num_sample)
    #Z_light_weighted_Z_band = np.zeros(num_sample)

    #FeH_median = np.zeros(num_sample)
    #MgFe_median = np.zeros(num_sample)
    #OFe_median = np.zeros(num_sample)
    FeH_log_mass_weighted = np.zeros(num_sample)
    FeH_log_mass_weighted_ratio = np.zeros(num_sample)
    FeH_mass_weighted_log = np.zeros(num_sample)

    FeH_log_light_weighted = np.zeros(num_sample)
    FeH_log_light_weighted_ratio = np.zeros(num_sample)
    FeH_light_weighted_log = np.zeros(num_sample)

    FeH_log_metallicity_weighted = np.zeros(num_sample)
    FeH_log_metallicity_weighted_ratio = np.zeros(num_sample)
    FeH_metallicity_weighted_log = np.zeros(num_sample)

    OFe_log_mass_weighted = np.zeros(num_sample)
    OFe_log_mass_weighted_ratio = np.zeros(num_sample)
    OFe_mass_weighted_log = np.zeros(num_sample)

    OFe_log_light_weighted = np.zeros(num_sample)
    OFe_log_light_weighted_ratio = np.zeros(num_sample)
    OFe_light_weighted_log = np.zeros(num_sample)

    OFe_log_metallicity_weighted = np.zeros(num_sample)
    OFe_log_metallicity_weighted_ratio = np.zeros(num_sample)
    OFe_metallicity_weighted_log = np.zeros(num_sample)

    MgFe_log_mass_weighted = np.zeros(num_sample)
    MgFe_log_mass_weighted_ratio = np.zeros(num_sample)
    MgFe_mass_weighted_log = np.zeros(num_sample)

    MgFe_log_light_weighted = np.zeros(num_sample)
    MgFe_log_light_weighted_ratio = np.zeros(num_sample)
    MgFe_light_weighted_log = np.zeros(num_sample)

    MgFe_log_metallicity_weighted = np.zeros(num_sample)
    MgFe_log_metallicity_weighted_ratio = np.zeros(num_sample)
    MgFe_metallicity_weighted_log = np.zeros(num_sample)

    stellar_mass = sim_info.halo_data.log10_stellar_mass[sample]

    for i in tqdm(range(len(sample))):
        sim_info.make_particle_data(halo_id=sim_info.halo_data.halo_ids[sample[i]])

        M_stars = sim_info.stars.mass.copy()
        Z_stars = sim_info.stars.metal_mass_fractions.copy()
        stars_mag = sim_info.calculate_luminosities(sim_info.stars.age, Z_stars, sim_info.stars.initmass)
        L_stars_r_band = pow(10.0, -0.4 * stars_mag["r"]) * 3631
        #L_stars_i_band = pow(10.0, -0.4 * stars_mag["i"]) * 3631
        #L_stars_z_band = pow(10.0, -0.4 * stars_mag["z"]) * 3631

        # First metallicity
        Z_mass_weighted[i] = (np.sum(Z_stars * M_stars) / np.sum(M_stars)) / Zsun
        Z_light_weighted[i] = (np.sum(Z_stars * L_stars_r_band) / np.sum(L_stars_r_band)) / Zsun

        Fe_stars = sim_info.stars.iron.copy()
        H_stars = sim_info.stars.hydrogen.copy()

        # Now [Fe/H], variations of mass weighted:
        ratios = compute_weighted_FeH(Fe_stars, H_stars, M_stars)
        FeH_log_mass_weighted[i] = ratios['log_of_weight']
        FeH_log_mass_weighted_ratio[i] = ratios['log_of_weighted_ratio']
        FeH_mass_weighted_log[i] = ratios['weight_of_log']

        # Now [Fe/H], variations of light weighted:
        ratios = compute_weighted_FeH(Fe_stars, H_stars, L_stars_r_band)
        FeH_log_light_weighted[i] = ratios['log_of_weight']
        FeH_log_light_weighted_ratio[i] = ratios['log_of_weighted_ratio']
        FeH_light_weighted_log[i] = ratios['weight_of_log']

        # Now [Fe/H], variations of metallicity weighted:
        ratios = compute_weighted_FeH(Fe_stars, H_stars, Z_stars)
        FeH_log_metallicity_weighted[i] = ratios['log_of_weight']
        FeH_log_metallicity_weighted_ratio[i] = ratios['log_of_weighted_ratio']
        FeH_metallicity_weighted_log[i] = ratios['weight_of_log']

        #####
        O_stars = sim_info.stars.oxygen.copy()
        Mg_stars = sim_info.stars.magnesium.copy()

        # Now [O/Fe], variations of mass weighted:
        ratios = compute_weighted_OFe(O_stars, Fe_stars, M_stars)
        OFe_log_mass_weighted[i] = ratios['log_of_weight']
        OFe_log_mass_weighted_ratio[i] = ratios['log_of_weighted_ratio']
        OFe_mass_weighted_log[i] = ratios['weight_of_log']

        # Now [O/Fe], variations of light weighted:
        ratios = compute_weighted_OFe(O_stars, Fe_stars, L_stars_r_band)
        OFe_log_light_weighted[i] = ratios['log_of_weight']
        OFe_log_light_weighted_ratio[i] = ratios['log_of_weighted_ratio']
        OFe_light_weighted_log[i] = ratios['weight_of_log']

        # Now [O/Fe], variations of metallicity weighted:
        ratios = compute_weighted_OFe(O_stars, Fe_stars, Z_stars)
        OFe_log_metallicity_weighted[i] = ratios['log_of_weight']
        OFe_log_metallicity_weighted_ratio[i] = ratios['log_of_weighted_ratio']
        OFe_metallicity_weighted_log[i] = ratios['weight_of_log']

        # Now [Mg/Fe], variations of mass weighted:
        ratios = compute_weighted_MgFe(Mg_stars, Fe_stars, M_stars)
        MgFe_log_mass_weighted[i] = ratios['log_of_weight']
        MgFe_log_mass_weighted_ratio[i] = ratios['log_of_weighted_ratio']
        MgFe_mass_weighted_log[i] = ratios['weight_of_log']

        # Now [Mg/Fe], variations of light weighted:
        ratios = compute_weighted_MgFe(Mg_stars, Fe_stars, L_stars_r_band)
        MgFe_log_light_weighted[i] = ratios['log_of_weight']
        MgFe_log_light_weighted_ratio[i] = ratios['log_of_weighted_ratio']
        MgFe_light_weighted_log[i] = ratios['weight_of_log']

        # Now [Mg/Fe], variations of metallicity weighted:
        ratios = compute_weighted_MgFe(Mg_stars, Fe_stars, Z_stars)
        MgFe_log_metallicity_weighted[i] = ratios['log_of_weight']
        MgFe_log_metallicity_weighted_ratio[i] = ratios['log_of_weighted_ratio']
        MgFe_metallicity_weighted_log[i] = ratios['weight_of_log']

        #
        #
        # ratios = compute_alpha_light_weighted_ratios(Fe_stars, O_stars, Mg_stars, M_stars, L_stars_r_band)
        #
        # OFe_light_weighted_r_band[i] = ratios['O_Fe_weighted_r_band']
        # MgFe_light_weighted_r_band[i] = ratios['Mg_Fe_weighted_r_band']
        #
        # OFe_light_weighted_r_band_with_mass[i] = ratios['O_Fe_weighted_r_band_with_mass']
        # MgFe_light_weighted_r_band_with_mass[i] = ratios['Mg_Fe_weighted_r_band_with_mass']
        #
        # ratios = compute_light_weighted_ratios(Z_stars, Fe_stars, H_stars, M_stars,
        #                                        L_stars_r_band, L_stars_i_band, L_stars_z_band)
        # Z_light_weighted_r_band[i] = ratios['Z_weighted_r_band']
        # Z_light_weighted_i_band[i] = ratios['Z_weighted_i_band']
        # Z_light_weighted_Z_band[i] = ratios['Z_weighted_Z_band']
        #
        # FeH_light_weighted_r_band[i] = ratios['Fe_H_weighted_r_band']
        # FeH_light_weighted_i_band[i] = ratios['Fe_H_weighted_i_band']
        # FeH_light_weighted_Z_band[i] = ratios['Fe_H_weighted_Z_band']
        #
        # FeH_light_weighted_r_band_with_mass[i] = ratios['Fe_H_weighted_r_band_with_mass']
        # FeH_light_weighted_i_band_with_mass[i] = ratios['Fe_H_weighted_i_band_with_mass']
        # FeH_light_weighted_Z_band_with_mass[i] = ratios['Fe_H_weighted_Z_band_with_mass']
        #
        #
        # ratios = compute_median_ratios(O_stars, Fe_stars, Mg_stars, H_stars)
        # FeH_median[i] = ratios['Fe_H']
        # MgFe_median[i] = ratios['Mg_Fe']
        # OFe_median[i] = ratios['O_Fe']
        #
        #
        # ratios = compute_mass_weighted_ratios(O_stars, Fe_stars, Mg_stars, H_stars, M_stars)
        # FeH_total_mass_weighted[i] = ratios['Fe_H']
        # MgFe_total_mass_weighted[i] = ratios['Mg_Fe']
        # OFe_total_mass_weighted[i] = ratios['O_Fe']
        # FeH_ratio_weighted[i] = ratios['Fe_H_weighted_ratio']
        # MgFe_ratio_weighted[i] = ratios['Mg_Fe_weighted_ratio']
        # OFe_ratio_weighted[i] = ratios['O_Fe_weighted_ratio']
        #
        # ratios = compute_metallicity_weighted_ratios(Z_stars, Fe_stars, H_stars, O_stars, Mg_stars, L_stars_r_band)
        # FeH_Z_light_weighted[i] = ratios['Fe_H_Z_weighted_r_band']
        # FeH_Z_weighted[i] = ratios['Fe_H_Z_weighted']
        # OFe_Z_light_weighted[i] = ratios['O_Fe_Z_weighted_r_band']
        # OFe_Z_weighted[i] = ratios['O_Fe_Z_weighted']
        # MgFe_Z_light_weighted[i] = ratios['Mg_Fe_Z_weighted_r_band']
        # MgFe_Z_weighted[i] = ratios['Mg_Fe_Z_weighted']

    return{
        'Z_mass_weighted': Z_mass_weighted,
        'Z_light_weighted': Z_light_weighted,
        'FeH_log_mass_weighted': FeH_log_mass_weighted,
        'FeH_log_mass_weighted_ratio': FeH_log_mass_weighted_ratio,
        'FeH_mass_weighted_log': FeH_mass_weighted_log,
        'FeH_log_light_weighted': FeH_log_light_weighted,
        'FeH_log_light_weighted_ratio': FeH_log_light_weighted_ratio,
        'FeH_light_weighted_log': FeH_light_weighted_log,
        'FeH_log_metallicity_weighted': FeH_log_metallicity_weighted,
        'FeH_log_metallicity_weighted_ratio': FeH_log_metallicity_weighted_ratio,
        'FeH_metallicity_weighted_log':FeH_metallicity_weighted_log,

        'OFe_log_mass_weighted': OFe_log_mass_weighted,
        'OFe_log_mass_weighted_ratio': OFe_log_mass_weighted_ratio,
        'OFe_mass_weighted_log': OFe_mass_weighted_log,
        'OFe_log_light_weighted': OFe_log_light_weighted,
        'OFe_log_light_weighted_ratio': OFe_log_light_weighted_ratio,
        'OFe_light_weighted_log': OFe_light_weighted_log,
        'OFe_log_metallicity_weighted': OFe_log_metallicity_weighted,
        'OFe_log_metallicity_weighted_ratio': OFe_log_metallicity_weighted_ratio,
        'OFe_metallicity_weighted_log': OFe_metallicity_weighted_log,

        'MgFe_log_mass_weighted': MgFe_log_mass_weighted,
        'MgFe_log_mass_weighted_ratio': MgFe_log_mass_weighted_ratio,
        'MgFe_mass_weighted_log': MgFe_mass_weighted_log,
        'MgFe_log_light_weighted': MgFe_log_light_weighted,
        'MgFe_log_light_weighted_ratio': MgFe_log_light_weighted_ratio,
        'MgFe_light_weighted_log': MgFe_light_weighted_log,
        'MgFe_log_metallicity_weighted': MgFe_log_metallicity_weighted,
        'MgFe_log_metallicity_weighted_ratio': MgFe_log_metallicity_weighted_ratio,
        'MgFe_metallicity_weighted_log': MgFe_metallicity_weighted_log,

        'Mstellar': stellar_mass
    }
    # return {'Fe_H_median': FeH_median,
    #         'O_Fe_median': OFe_median,
    #         'Mg_Fe_median': MgFe_median,
    #
    #         'Fe_H_mass_weighted': FeH_total_mass_weighted,
    #         'O_Fe_mass_weighted': OFe_total_mass_weighted,
    #         'Mg_Fe_mass_weighted': MgFe_total_mass_weighted,
    #
    #         'Fe_H_ratio_weighted': FeH_ratio_weighted,
    #         'O_Fe_ratio_weighted': OFe_ratio_weighted,
    #         'Mg_Fe_ratio_weighted': MgFe_ratio_weighted,
    #
    #         'Z_mass_weighted': Z_mass_weighted,
    #         'Z_light_weighted_r_band': Z_light_weighted_r_band,
    #         'Z_light_weighted_i_band': Z_light_weighted_i_band,
    #         'Z_light_weighted_Z_band': Z_light_weighted_Z_band,
    #         'Fe_H_light_weighted_r_band': FeH_light_weighted_r_band,
    #         'Fe_H_light_weighted_i_band': FeH_light_weighted_i_band,
    #         'Fe_H_light_weighted_Z_band': FeH_light_weighted_Z_band,
    #
    #         'Fe_H_light_weighted_r_band_with_mass': FeH_light_weighted_r_band_with_mass,
    #         'Fe_H_light_weighted_i_band_with_mass': FeH_light_weighted_i_band_with_mass,
    #         'Fe_H_light_weighted_Z_band_with_mass': FeH_light_weighted_Z_band_with_mass,
    #
    #         'O_Fe_light_weighted_r_band': OFe_light_weighted_r_band,
    #         'Mg_Fe_light_weighted_r_band': MgFe_light_weighted_r_band,
    #
    #         'O_Fe_light_weighted_r_band_with_mass': OFe_light_weighted_r_band_with_mass,
    #         'Mg_Fe_light_weighted_r_band_with_mass': MgFe_light_weighted_r_band_with_mass,
    #
    #         'Fe_H_Z_weighted_r_band': FeH_Z_light_weighted,
    #         'Fe_H_Z_weighted': FeH_Z_weighted,
    #         'O_Fe_Z_weighted_r_band': OFe_Z_light_weighted,
    #         'O_Fe_Z_weighted': OFe_Z_weighted,
    #         'Mg_Fe_Z_weighted_r_band': MgFe_Z_light_weighted,
    #         'Mg_Fe_Z_weighted': MgFe_Z_weighted,
    #
    #         'Mstellar':stellar_mass}

def plot_Kirby_distributions(output_path):

    kroupa_to_chabrier_mass = 0.912
    delimiter = "\t"

    # Read the data
    input_filename = "./plotter/obs_data/Kirby_2013_data.txt"
    galaxy_name = np.loadtxt(input_filename, usecols=[0], dtype=str)
    galaxy_mass = np.loadtxt(input_filename, usecols=[1])

    input_filename = "./plotter/obs_data/Kirby_2013_individual_stars.txt"
    raw_galaxy_list = np.loadtxt(input_filename, usecols=[0], dtype=str)
    raw = np.loadtxt(input_filename, usecols=[2, 3])

    num_galaxies = len(galaxy_mass)

    # Plot parameters
    params = {
        "font.size": 12,
        "font.family": "STIXGeneral",
        #"font.family": "Times",
        "text.usetex": False,
        "mathtext.fontset": "stix",
        "figure.figsize": (4, 3),
        "figure.subplot.left": 0.18,
        "figure.subplot.right": 0.95,
        "figure.subplot.bottom": 0.15,
        "figure.subplot.top": 0.95,
        "figure.subplot.wspace": 0.3,
        "figure.subplot.hspace": 0.3,
        "lines.markersize": 0.5,
        "lines.linewidth": 0.2,
    }
    rcParams.update(params)

    plt.figure()
    ax = plt.subplot(1, 1, 1)
    plt.grid(True)

    for i in range(num_galaxies):

        pick = np.where(raw_galaxy_list == galaxy_name[i])[0]
        if len(pick) == 0: continue

        Fe_H_data = 10 ** raw[pick, 0]

        bins = np.arange(-7.2,3,0.2)
        bins = 10 ** bins
        hist, _, _ = stat.binned_statistic(Fe_H_data , values=np.ones(len(Fe_H_data)), statistic="sum", bins=bins)
        bin_centers = 0.5*(bins[1:] + bins[:-1])
        hist = hist/np.sum(hist)
        plt.plot(bin_centers, hist, '-', lw=1, label=galaxy_name[i])

    plt.xlabel("[Fe/H]", labelpad=2)
    plt.ylabel("PDF", labelpad=2)
    plt.xscale('log')
    plt.yscale('log')
    ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    plt.axis([2e-4, 2e0, 1e-3, 1e0])
    plt.legend(loc='upper left', labelspacing=0.1, handlelength=1.5,
               handletextpad=0.1, frameon=False, ncol=2,
               fontsize=9, columnspacing=0.02)

    plt.savefig(f"{output_path}/Kirby_distribution.png", dpi=200)


def plot_Kirby_analysed():

    kroupa_to_chabrier_mass = 0.912
    delimiter = "\t"

    # Read the data
    input_filename = "./plotter/obs_data/Kirby_2013_data.txt"
    galaxy_name = np.loadtxt(input_filename, usecols=[0], dtype=str)
    galaxy_mass = np.loadtxt(input_filename, usecols=[1])

    input_filename = "./plotter/obs_data/Kirby_2013_individual_stars.txt"
    raw_galaxy_list = np.loadtxt(input_filename, usecols=[0], dtype=str)
    raw = np.loadtxt(input_filename, usecols=[2, 3])

    num_galaxies = len(galaxy_mass)
    M_star = []
    Z_star = []
    Z_star_lo = []
    Z_star_hi = []

    for i in range(num_galaxies):

        pick = np.where(raw_galaxy_list == galaxy_name[i])[0]
        if len(pick) == 0: continue

        M_star = np.append(M_star, 10 ** galaxy_mass[i] * kroupa_to_chabrier_mass)

        Fe_H_data = 10 ** raw[pick, 0]

        Fe_H = np.median(Fe_H_data)
        Fe_H_lo = np.percentile(Fe_H_data, 16)
        Fe_H_hi = np.percentile(Fe_H_data, 84)

        Z_star = np.append(Z_star, Fe_H)
        Z_star_lo = np.append(Z_star_lo, Fe_H_lo)
        Z_star_hi = np.append(Z_star_hi, Fe_H_hi)

    # Define the scatter as offset from the mean value
    y_scatter = np.array((Z_star - Z_star_lo, Z_star_hi - Z_star))
    plt.errorbar(M_star, Z_star, yerr= y_scatter, fmt='*',
                 ls='none',color='black',lw=1,ms=2.5,label='Kirby et al. (2013)/Analysed')


def plot_Kirby_data():

    kroupa_to_chabrier_mass = 0.912

    input_filename = "./plotter/obs_data/Kirby_2013_ascii.dat"
    delimiter = "\t"

    # Read the data
    raw = np.loadtxt(input_filename, delimiter=delimiter)
    M_star = 10 ** raw[:, 0] * kroupa_to_chabrier_mass
    M_star_lo = 10 ** (raw[:, 0] - raw[:, 1]) * kroupa_to_chabrier_mass
    M_star_hi = 10 ** (raw[:, 0] + raw[:, 2]) * kroupa_to_chabrier_mass

    Z_star = 10 ** raw[:, 3]  # Z/Z_sun
    Z_star_lo = 10 ** (raw[:, 3] - raw[:, 4]) # Z/Z_sun
    Z_star_hi = 10 ** (raw[:, 3] + raw[:, 5]) # Z/Z_sun

    # Define the scatter as offset from the mean value
    x_scatter = np.array((M_star - M_star_lo, M_star_hi - M_star))
    y_scatter = np.array((Z_star - Z_star_lo, Z_star_hi - Z_star))
    plt.errorbar(M_star, Z_star, yerr= y_scatter, xerr=x_scatter, fmt='',
                 ls='none',color='lightgreen',lw=1,ms=1,label='Kirby et al. (2013)')

def plot_gallazzi(element):
    # Cosmology
    h_sim = 0.6777
    h_obs = 0.7 # WMAP 7
    Z_solar_obs = 0.02

    input_filename = "./plotter/obs_data/gallazzi_2021_ascii.txt"
    delimiter = "\t"

    # Read the data
    raw = np.loadtxt(input_filename, delimiter=delimiter)
    M_star = (
#            10 ** raw[:, 0] * h_obs ** 2
            10 ** raw[:, 0] * h_sim ** -2
#            10 ** raw[:, 0] * (h_sim / h_obs) ** -2
    )

    # Correction factor due to the difference in (X_O/X_Fe)_Sun
    # From Grevesse, Noels & Sauval (1996) to Asplund+ (2009)

    O_over_H_Gr96 = 8.87
    Fe_over_H_Gr96 = 7.50
    Mg_over_H_Gr96 = 7.58

    O_over_H_Asplund09 = 8.69
    Fe_over_H_Asplund09 = 7.50
    Mg_over_H_Asplund09 = 7.6

    O_over_Fe_solar_Gr96 = O_over_H_Gr96 - Fe_over_H_Gr96
    Mg_over_Fe_solar_Gr96 = Mg_over_H_Gr96 - Fe_over_H_Gr96
    O_over_Fe_solar_Asplund09 = O_over_H_Asplund09 - Fe_over_H_Asplund09
    Mg_over_Fe_solar_Asplund09 = Mg_over_H_Asplund09 - Fe_over_H_Asplund09

    if element == 'O':
        mH_in_cgs = 1
        mFe_in_cgs = 55.845  # * mp (times proton mass)
        mO_in_cgs = 15.999  # * mp

        O_Fe_Sun = O_over_H_Asplund09 - Fe_over_H_Asplund09 - np.log10(mFe_in_cgs / mO_in_cgs)

        # From Grevesse, Noels & Sauval (1996):
        O_Fe_Sun_Gr96 = np.log10(0.0096 / 0.0016)

        correction = O_Fe_Sun_Gr96 - O_Fe_Sun

        #correction = O_over_Fe_solar_Gr96 - O_over_Fe_solar_Asplund09

    if element == 'Mg':

        correction = Mg_over_Fe_solar_Gr96 - Mg_over_Fe_solar_Asplund09

    Z_median = (raw[:, 1] + correction)
    Z_lo = (raw[:, 2] + correction)
    Z_hi = (raw[:, 3] + correction)

    # Define the scatter as offset from the mean value
    y_scatter = np.array((Z_median - Z_lo, Z_hi - Z_median))
    plt.errorbar(M_star, Z_median, yerr= y_scatter, fmt='o',
                 ls='none',color='purple',lw=1,ms=1.5,label='Gallazzi et al. (2021)', zorder=23)

def plot_thomas_2010():

    mH_in_cgs = 1
    mFe_in_cgs = 55.845 # * mp (times proton mass)
    mO_in_cgs = 15.999 # * mp

    # Asplund et al. (2009)
    Fe_H_Sun_Asplund = 7.5
    O_H_Sun_Asplund = 8.69

    O_Fe_Sun = O_H_Sun_Asplund - Fe_H_Sun_Asplund - np.log10(mFe_in_cgs / mO_in_cgs)
    Fe_H_Sun = Fe_H_Sun_Asplund - 12.0 - np.log10(mH_in_cgs / mFe_in_cgs)

    # From Grevesse, Noels & Sauval (1996):
    alpha_Fe_Sun_Gr96 = np.log10(0.0148 / 0.0016)
    O_Fe_Sun_Gr96 = np.log10(0.0096 / 0.0016)

    correction_Sun_O_over_Fe = O_Fe_Sun_Gr96 - O_Fe_Sun

    # IMF correction factor (from M* with Kroupa IMF to M* with Chabrier IMF)
    # (Lacey+ 2016, https://arxiv.org/pdf/1509.08473, table B2)
    kroupa_to_chabrier_mass = 0.912

    # Conversion from log velocity dispersion to [alpha/Fe]
    # (eq. 1 in Thomas et al. 2010)
    func_alpha_over_Fe = (
        lambda log_sigma: -0.55 + 0.33 * log_sigma + correction_Sun_O_over_Fe
    )

    # Conversion from log Mstar to log velocity dispersion (eq. 2 in
    # Thomas et al. 2005, https://arxiv.org/pdf/astro-ph/0410209)
    func_log_sigma = lambda log_Mstar: (log_Mstar - 0.63) / 4.52

    # X values of the final data
    M_star = np.arange(np.log10(3e10), np.log10(3e11), 0.1)

    # Conversions
    M_star = 10**M_star * kroupa_to_chabrier_mass
    sigma_log = func_log_sigma(log_Mstar=np.log10(M_star))

    # Y values of the final data
    alpha_over_Fe = func_alpha_over_Fe(log_sigma=sigma_log)

    plt.plot(M_star, alpha_over_Fe,'-', color='grey',lw=1.5, label='Thomas et al. (2010)', zorder=25)


def plot_conroy_2014():

    M_star = np.array([9.8, 10.08, 10.55, 10.70, 11.07])
    O_Fe = np.array([0.10, 0.17, 0.20, 0.25, 0.28])

    # IMF correction factor (from M* with Kroupa IMF to M* with Chabrier IMF)
    # (Lacey+ 2016, https://arxiv.org/pdf/1509.08473, table B2)
    kroupa_to_chabrier_mass = 0.912

    M_star = 10**M_star * kroupa_to_chabrier_mass

    plt.plot(M_star, O_Fe,'-', color='lightgreen',lw=1.5, label='Conroy et al. (2014)', zorder=24)


def plot_gallazzi_2005():
    # Cosmology
    h_sim = 0.6777
    h_obs = 0.7
    Z_solar_obs = 0.02
    kroupa_to_chabrier_mass = 0.912
    solar_metallicity = 0.0126

    input_filename = "./plotter/obs_data/Gallazzi_2005_ascii.txt"
    delimiter = "\t"

    # Read the data
    raw = np.loadtxt(input_filename, delimiter=delimiter)
    M_star = (
            10 ** raw[:, 0] * (h_sim / h_obs) ** -2 * kroupa_to_chabrier_mass
    )
    Z_median = 10 ** raw[:, 1] * Z_solar_obs / solar_metallicity
    Z_lo = 10 ** raw[:, 2] * Z_solar_obs / solar_metallicity
    Z_hi = 10 ** raw[:, 3] * Z_solar_obs / solar_metallicity

    # Define the scatter as offset from the mean value
    y_scatter = np.array((Z_median - Z_lo, Z_hi - Z_median))

    plt.errorbar(M_star, Z_median, yerr= y_scatter, fmt='o',
                 ls='none',color='grey',lw=1,ms=1,label='Gallazzi et al. (2005)', zorder=10)


def plot_Zahid_2017():

    mass_bins = np.array([8.55, 8.65, 8.75, 8.85, 8.95, 9.05, 9.15, 9.25, 9.35, 9.45, 9.55, 9.65, 9.75, 9.85,
                          9.95, 10.05, 10.15, 10.25, 10.35, 10.45, 10.55, 10.65, 10.75, 10.85])

    input_filename = "./plotter/obs_data/Zahid_2017.txt"

    # Read the data
    raw = np.loadtxt(input_filename)
    M_star = raw[:, 1]
    Z_star = 10 ** raw[:, 3]

    Z_median = np.zeros(len(mass_bins))
    for i in range(len(mass_bins)):
        select = np.where(M_star == mass_bins[i])[0]
        Z_median[i] = np.median(Z_star[select])

    plt.plot(10**mass_bins, Z_median, '+',color='purple',ms=5,label='Zahid et al. (2017)')


def plot_Kudritzki_2016():

    input_filename = "./plotter/obs_data/Kudritzki_2016.txt"

    # Read the data
    raw = np.loadtxt(input_filename)
    M_star = 10 ** raw[:, 0]
    Z_median = 10 ** raw[:, 1]

    plt.plot(M_star, Z_median, 'v',color='black',ms=3,label='Kudritzki et al. (2016)')

def compute_metallicity_relation(sim_info, metallicity_data):

    # Look for abundance ratios from COLIBRE snaps:
    data = calculate_galaxies_metallicity(sim_info)

    Mstellar_median = data['Mstellar']
    Z_mass_weighted = data['Z_mass_weighted']
    Z_light_weighted = data['Z_light_weighted']

    FeH_log_mass_weighted = data['FeH_log_mass_weighted']
    FeH_log_mass_weighted_ratio = data['FeH_log_mass_weighted_ratio']
    FeH_mass_weighted_log = data['FeH_mass_weighted_log']
    FeH_log_light_weighted = data['FeH_log_light_weighted']
    FeH_log_light_weighted_ratio = data['FeH_log_light_weighted_ratio']
    FeH_light_weighted_log = data['FeH_light_weighted_log']
    FeH_log_metallicity_weighted = data['FeH_log_metallicity_weighted']
    FeH_log_metallicity_weighted_ratio = data['FeH_log_metallicity_weighted_ratio']
    FeH_metallicity_weighted_log = data['FeH_metallicity_weighted_log']

    OFe_log_mass_weighted = data['OFe_log_mass_weighted']
    OFe_log_mass_weighted_ratio = data['OFe_log_mass_weighted_ratio']
    OFe_mass_weighted_log = data['OFe_mass_weighted_log']
    OFe_log_light_weighted = data['OFe_log_light_weighted']
    OFe_log_light_weighted_ratio = data['OFe_log_light_weighted_ratio']
    OFe_light_weighted_log = data['OFe_light_weighted_log']
    OFe_log_metallicity_weighted = data['OFe_log_metallicity_weighted']
    OFe_log_metallicity_weighted_ratio = data['OFe_log_metallicity_weighted_ratio']
    OFe_metallicity_weighted_log = data['OFe_metallicity_weighted_log']

    MgFe_log_mass_weighted = data['MgFe_log_mass_weighted']
    MgFe_log_mass_weighted_ratio = data['MgFe_log_mass_weighted_ratio']
    MgFe_mass_weighted_log = data['MgFe_mass_weighted_log']
    MgFe_log_light_weighted = data['MgFe_log_light_weighted']
    MgFe_log_light_weighted_ratio = data['MgFe_log_light_weighted_ratio']
    MgFe_light_weighted_log = data['MgFe_light_weighted_log']
    MgFe_log_metallicity_weighted = data['MgFe_log_metallicity_weighted']
    MgFe_log_metallicity_weighted_ratio = data['MgFe_log_metallicity_weighted_ratio']
    MgFe_metallicity_weighted_log = data['MgFe_metallicity_weighted_log']

    counter = np.array([len(Mstellar_median)])

    if metallicity_data == None:

        metallicity_data = {
            'Z_mass_weighted': Z_mass_weighted,
            'Z_light_weighted': Z_light_weighted,

            'FeH_log_mass_weighted': FeH_log_mass_weighted,
            'FeH_log_mass_weighted_ratio': FeH_log_mass_weighted_ratio,
            'FeH_mass_weighted_log': FeH_mass_weighted_log,
            'FeH_log_light_weighted': FeH_log_light_weighted,
            'FeH_log_light_weighted_ratio': FeH_log_light_weighted_ratio,
            'FeH_light_weighted_log': FeH_light_weighted_log,
            'FeH_log_metallicity_weighted': FeH_log_metallicity_weighted,
            'FeH_log_metallicity_weighted_ratio': FeH_log_metallicity_weighted_ratio,
            'FeH_metallicity_weighted_log': FeH_metallicity_weighted_log,

            'OFe_log_mass_weighted': OFe_log_mass_weighted,
            'OFe_log_mass_weighted_ratio': OFe_log_mass_weighted_ratio,
            'OFe_mass_weighted_log': OFe_mass_weighted_log,
            'OFe_log_light_weighted': OFe_log_light_weighted,
            'OFe_log_light_weighted_ratio': OFe_log_light_weighted_ratio,
            'OFe_light_weighted_log': OFe_light_weighted_log,
            'OFe_log_metallicity_weighted': OFe_log_metallicity_weighted,
            'OFe_log_metallicity_weighted_ratio': OFe_log_metallicity_weighted_ratio,
            'OFe_metallicity_weighted_log': OFe_metallicity_weighted_log,

            'MgFe_log_mass_weighted': MgFe_log_mass_weighted,
            'MgFe_log_mass_weighted_ratio': MgFe_log_mass_weighted_ratio,
            'MgFe_mass_weighted_log': MgFe_mass_weighted_log,
            'MgFe_log_light_weighted': MgFe_log_light_weighted,
            'MgFe_log_light_weighted_ratio': MgFe_log_light_weighted_ratio,
            'MgFe_light_weighted_log': MgFe_light_weighted_log,
            'MgFe_log_metallicity_weighted': MgFe_log_metallicity_weighted,
            'MgFe_log_metallicity_weighted_ratio': MgFe_log_metallicity_weighted_ratio,
            'MgFe_metallicity_weighted_log': MgFe_metallicity_weighted_log,

            'Mstellar': Mstellar_median,
            'counter':counter}

    else:

        Mstellar = metallicity_data['Mstellar']
        Mstellar = np.append(Mstellar, Mstellar_median)
        counter_sim = metallicity_data['counter']
        counter_sim = np.append(counter_sim, counter)

        Z_mw = metallicity_data['Z_mass_weighted']
        Z_mw = np.append(Z_mw, Z_mass_weighted)
        Z_lw = metallicity_data['Z_light_weighted']
        Z_lw = np.append(Z_lw, Z_light_weighted)

        FeH_lw = metallicity_data['FeH_log_light_weighted']
        FeH_lw = np.append(FeH_lw, FeH_log_light_weighted)
        FeH_lw_r = metallicity_data['FeH_log_light_weighted_ratio']
        FeH_lw_r = np.append(FeH_lw_r, FeH_log_light_weighted_ratio)
        FeH_lw_2 = metallicity_data['FeH_light_weighted_log']
        FeH_lw_2 = np.append(FeH_lw_2, FeH_light_weighted_log)

        FeH_mw = metallicity_data['FeH_log_mass_weighted']
        FeH_mw = np.append(FeH_mw, FeH_log_mass_weighted)
        FeH_mw_r = metallicity_data['FeH_log_mass_weighted_ratio']
        FeH_mw_r = np.append(FeH_mw_r, FeH_log_mass_weighted_ratio)
        FeH_mw_2 = metallicity_data['FeH_mass_weighted_log']
        FeH_mw_2 = np.append(FeH_mw_2, FeH_mass_weighted_log)

        FeH_zw = metallicity_data['FeH_log_metallicity_weighted']
        FeH_zw = np.append(FeH_zw, FeH_log_metallicity_weighted)
        FeH_zw_r = metallicity_data['FeH_log_metallicity_weighted_ratio']
        FeH_zw_r = np.append(FeH_zw_r, FeH_log_metallicity_weighted_ratio)
        FeH_zw_2 = metallicity_data['FeH_metallicity_weighted_log']
        FeH_zw_2 = np.append(FeH_zw_2, FeH_metallicity_weighted_log)

        OFe_lw = metallicity_data['OFe_log_light_weighted']
        OFe_lw = np.append(OFe_lw, OFe_log_light_weighted)
        OFe_lw_r = metallicity_data['OFe_log_light_weighted_ratio']
        OFe_lw_r = np.append(OFe_lw_r, OFe_log_light_weighted_ratio)
        OFe_lw_2 = metallicity_data['OFe_light_weighted_log']
        OFe_lw_2 = np.append(OFe_lw_2, OFe_light_weighted_log)

        OFe_mw = metallicity_data['OFe_log_mass_weighted']
        OFe_mw = np.append(OFe_mw, OFe_log_mass_weighted)
        OFe_mw_r = metallicity_data['OFe_log_mass_weighted_ratio']
        OFe_mw_r = np.append(OFe_mw_r, OFe_log_mass_weighted_ratio)
        OFe_mw_2 = metallicity_data['OFe_mass_weighted_log']
        OFe_mw_2 = np.append(OFe_mw_2, OFe_mass_weighted_log)

        OFe_zw = metallicity_data['OFe_log_metallicity_weighted']
        OFe_zw = np.append(OFe_zw, OFe_log_metallicity_weighted)
        OFe_zw_r = metallicity_data['OFe_log_metallicity_weighted_ratio']
        OFe_zw_r = np.append(OFe_zw_r, OFe_log_metallicity_weighted_ratio)
        OFe_zw_2 = metallicity_data['OFe_metallicity_weighted_log']
        OFe_zw_2 = np.append(OFe_zw_2, OFe_metallicity_weighted_log)

        MgFe_lw = metallicity_data['MgFe_log_light_weighted']
        MgFe_lw = np.append(MgFe_lw, MgFe_log_light_weighted)
        MgFe_lw_r = metallicity_data['MgFe_log_light_weighted_ratio']
        MgFe_lw_r = np.append(MgFe_lw_r, MgFe_log_light_weighted_ratio)
        MgFe_lw_2 = metallicity_data['MgFe_light_weighted_log']
        MgFe_lw_2 = np.append(MgFe_lw_2, MgFe_light_weighted_log)

        MgFe_mw = metallicity_data['MgFe_log_mass_weighted']
        MgFe_mw = np.append(MgFe_mw, MgFe_log_mass_weighted)
        MgFe_mw_r = metallicity_data['MgFe_log_mass_weighted_ratio']
        MgFe_mw_r = np.append(MgFe_mw_r, MgFe_log_mass_weighted_ratio)
        MgFe_mw_2 = metallicity_data['MgFe_mass_weighted_log']
        MgFe_mw_2 = np.append(MgFe_mw_2, MgFe_mass_weighted_log)

        MgFe_zw = metallicity_data['MgFe_log_metallicity_weighted']
        MgFe_zw = np.append(MgFe_zw, MgFe_log_metallicity_weighted)
        MgFe_zw_r = metallicity_data['MgFe_log_metallicity_weighted_ratio']
        MgFe_zw_r = np.append(MgFe_zw_r, MgFe_log_metallicity_weighted_ratio)
        MgFe_zw_2 = metallicity_data['MgFe_metallicity_weighted_log']
        MgFe_zw_2 = np.append(MgFe_zw_2, MgFe_metallicity_weighted_log)

        metallicity_data = {
            'Z_mass_weighted': Z_mw,
            'Z_light_weighted': Z_lw,
            'FeH_log_mass_weighted': FeH_mw,
            'FeH_log_mass_weighted_ratio': FeH_mw_r,
            'FeH_mass_weighted_log': FeH_mw_2,
            'FeH_log_light_weighted': FeH_lw,
            'FeH_log_light_weighted_ratio': FeH_lw_r,
            'FeH_light_weighted_log': FeH_lw_2,
            'FeH_log_metallicity_weighted': FeH_zw,
            'FeH_log_metallicity_weighted_ratio': FeH_zw_r,
            'FeH_metallicity_weighted_log': FeH_zw_2,

            'OFe_log_mass_weighted': OFe_mw,
            'OFe_log_mass_weighted_ratio': OFe_mw_r,
            'OFe_mass_weighted_log': OFe_mw_2,
            'OFe_log_light_weighted': OFe_lw,
            'OFe_log_light_weighted_ratio': OFe_lw_r,
            'OFe_light_weighted_log': OFe_lw_2,
            'OFe_log_metallicity_weighted': OFe_zw,
            'OFe_log_metallicity_weighted_ratio': OFe_zw_r,
            'OFe_metallicity_weighted_log': OFe_zw_2,

            'MgFe_log_mass_weighted': MgFe_mw,
            'MgFe_log_mass_weighted_ratio': MgFe_mw_r,
            'MgFe_mass_weighted_log': MgFe_mw_2,
            'MgFe_log_light_weighted': MgFe_lw,
            'MgFe_log_light_weighted_ratio': MgFe_lw_r,
            'MgFe_light_weighted_log': MgFe_lw_2,
            'MgFe_log_metallicity_weighted': MgFe_zw,
            'MgFe_log_metallicity_weighted_ratio': MgFe_zw_r,
            'MgFe_metallicity_weighted_log': MgFe_zw_2,


            'Mstellar': Mstellar,
            'counter': counter_sim}

    return metallicity_data
