import matplotlib
import matplotlib.pylab as plt
from matplotlib.pylab import rcParams
import numpy as np
from tqdm import tqdm
from .loadObservationalData import plot_MW_Satellites, plot_MW_data, \
    plot_GALAH_data, plot_StrontiumObsData, plot_APOGEE_data
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


def compute_ratios(Hydrogen_fraction, Magnesium_fraction, Oxygen_fraction, Iron_fraction,
                   Carbon_fraction, Silicon_fraction, Europium_fraction, Barium_fraction,
                   Strontium_fraction, Nitrogen_fraction, Neon_fraction, Iron_fraction_SNIa):
    
    mp_in_cgs = 1.6737236e-24
    mH_in_cgs = 1.00784 * mp_in_cgs
    mFe_in_cgs = 55.845 * mp_in_cgs
    mO_in_cgs = 15.999 * mp_in_cgs
    mMg_in_cgs = 24.305 * mp_in_cgs
    mC_in_cgs =  12.0107 * mp_in_cgs
    mSi_in_cgs = 28.0855 * mp_in_cgs
    mEu_in_cgs = 151.964 * mp_in_cgs
    mBa_in_cgs = 137.327 * mp_in_cgs
    mSr_in_cgs = 87.62 * mp_in_cgs
    mN_in_cgs = 14.0067 * mp_in_cgs
    mNe_in_cgs = 20.1797 * mp_in_cgs

    # Asplund et al. (2009)
    Fe_H_Sun = 7.5
    O_H_Sun = 8.69
    Mg_H_Sun = 7.6
    C_H_Sun = 8.43
    Si_H_Sun = 7.51
    Eu_H_Sun = 0.52
    Ba_H_Sun = 2.18
    Sr_H_Sun = 2.87
    N_H_Sun = 7.83
    Ne_H_Sun = 7.93

    O_Fe_Sun = O_H_Sun - Fe_H_Sun - np.log10(mFe_in_cgs / mO_in_cgs)
    Mg_Fe_Sun = Mg_H_Sun - Fe_H_Sun - np.log10(mFe_in_cgs / mMg_in_cgs)
    C_Fe_Sun = C_H_Sun - Fe_H_Sun - np.log10(mFe_in_cgs / mC_in_cgs)
    Si_Fe_Sun = Si_H_Sun - Fe_H_Sun - np.log10(mFe_in_cgs / mSi_in_cgs)
    Eu_Fe_Sun = Eu_H_Sun - Fe_H_Sun - np.log10(mFe_in_cgs / mEu_in_cgs)
    Ba_Fe_Sun = Ba_H_Sun - Fe_H_Sun - np.log10(mFe_in_cgs / mBa_in_cgs)
    Sr_Fe_Sun = Sr_H_Sun - Fe_H_Sun - np.log10(mFe_in_cgs / mSr_in_cgs)
    N_Fe_Sun = N_H_Sun - Fe_H_Sun - np.log10(mFe_in_cgs / mN_in_cgs)
    Ne_Fe_Sun = Ne_H_Sun - Fe_H_Sun - np.log10(mFe_in_cgs / mNe_in_cgs)
    Fe_H_Sun = Fe_H_Sun - 12.0 - np.log10(mH_in_cgs / mFe_in_cgs)

    Fe_H = np.log10(Iron_fraction / Hydrogen_fraction) - Fe_H_Sun
    O_Fe = np.log10(Oxygen_fraction / Iron_fraction) - O_Fe_Sun
    Mg_Fe = np.log10(Magnesium_fraction / Iron_fraction) - Mg_Fe_Sun
    C_Fe = np.log10(Carbon_fraction / Iron_fraction) - C_Fe_Sun
    Si_Fe = np.log10(Silicon_fraction / Iron_fraction) - Si_Fe_Sun
    Eu_Fe = np.log10(Europium_fraction / Iron_fraction) - Eu_Fe_Sun
    Ba_Fe = np.log10(Barium_fraction / Iron_fraction) - Ba_Fe_Sun
    Sr_Fe = np.log10(Strontium_fraction / Iron_fraction) - Sr_Fe_Sun
    N_Fe = np.log10(Nitrogen_fraction / Iron_fraction) - N_Fe_Sun
    Ne_Fe = np.log10(Neon_fraction / Iron_fraction) - Ne_Fe_Sun
    FeSNIa_Fe = Iron_fraction_SNIa / Iron_fraction

    # Let's set lower and upper limits:
    Fe_H[Iron_fraction == 0] = -7  # set lower limit
    Fe_H[Fe_H < -7] = -7  # set lower limit
    Mg_Fe[Iron_fraction == 0] = -2  # set lower limit
    Mg_Fe[Magnesium_fraction == 0] = -2  # set lower limit
    Mg_Fe[Mg_Fe < -2] = -2  # set lower limit
    O_Fe[Iron_fraction == 0] = -2  # set lower limit
    O_Fe[Oxygen_fraction == 0] = -2  # set lower limit
    O_Fe[O_Fe < -2] = -2  # set lower limit
    C_Fe[Iron_fraction == 0] = -2  # set lower limit
    C_Fe[Carbon_fraction == 0] = -2  # set lower limit
    C_Fe[C_Fe < -2] = -2  # set lower limit
    Si_Fe[Iron_fraction == 0] = -2  # set lower limit
    Si_Fe[Silicon_fraction == 0] = -2  # set lower limit
    Si_Fe[Si_Fe < -2] = -2  # set lower limit
    Eu_Fe[Iron_fraction == 0] = -2  # set lower limit
    Eu_Fe[Europium_fraction == 0] = -2  # set lower limit
    Eu_Fe[Eu_Fe < -2] = -2  # set lower limit
    Ba_Fe[Iron_fraction == 0] = -2  # set lower limit
    Ba_Fe[Barium_fraction == 0] = -2  # set lower limit
    Ba_Fe[Ba_Fe < -2] = -2  # set lower limit
    Sr_Fe[Iron_fraction == 0] = -2  # set lower limit
    Sr_Fe[Strontium_fraction == 0] = -2  # set lower limit
    Sr_Fe[Sr_Fe < -2] = -2  # set lower limit
    N_Fe[Iron_fraction == 0] = -2  # set lower limit
    N_Fe[Nitrogen_fraction == 0] = -2  # set lower limit
    N_Fe[N_Fe < -2] = -2  # set lower limit
    Ne_Fe[Iron_fraction == 0] = -2  # set lower limit
    Ne_Fe[Neon_fraction == 0] = -2  # set lower limit
    Ne_Fe[Ne_Fe < -2] = -2  # set lower limit

    return {'Fe_H':Fe_H, 'O_Fe':O_Fe, 'Mg_Fe':Mg_Fe,
            'C_Fe':C_Fe, 'Si_Fe':Si_Fe, 'Eu_Fe':Eu_Fe,
            'Ba_Fe':Ba_Fe, 'Sr_Fe':Sr_Fe, 'N_Fe':N_Fe,
            'Ne_Fe':Ne_Fe, 'FeSNIa_Fe': FeSNIa_Fe}

def calculate_abundaces_from_MW_type_galaxies(sim_info):

    select_MW_mass = np.where((sim_info.halo_data.log10_stellar_mass >= 10.0) &
                              (sim_info.halo_data.log10_stellar_mass <= 11.0))[0]
    select_centrals = np.where(sim_info.halo_data.type[select_MW_mass] == 10)[0]

    sample = select_MW_mass[select_centrals]

    Oxygen_fraction = []
    Iron_fraction = []
    Magnesium_fraction = []
    Hydrogen_fraction = []
    Carbon_fraction = []
    Silicon_fraction = []
    Europium_fraction = []
    Strontium_fraction = []
    Barium_fraction = []
    Nitrogen_fraction = []
    Neon_fraction = []
    Iron_fraction_SNIa = []

    stars_R = []
    stars_z = []
    stars_age = []

    for i in tqdm(range(len(sample))):

        sim_info.make_particle_data(halo_id=sim_info.halo_data.halo_ids[sample[i]])

        O_stars = sim_info.stars.oxygen
        Fe_stars = sim_info.stars.iron
        Mg_stars = sim_info.stars.magnesium
        H_stars = sim_info.stars.hydrogen
        C_stars = sim_info.stars.carbon
        Si_stars = sim_info.stars.silicon
        Eu_stars = sim_info.stars.europium
        Ba_stars = sim_info.stars.barium
        Sr_stars = sim_info.stars.strontium
        N_stars = sim_info.stars.nitrogen
        Ne_stars = sim_info.stars.neon
        FeSNIa_Fe = sim_info.stars.iron_SNIa_fraction

        Oxygen_fraction = np.append(Oxygen_fraction, O_stars)
        Iron_fraction = np.append(Iron_fraction, Fe_stars)
        Magnesium_fraction = np.append(Magnesium_fraction, Mg_stars)
        Hydrogen_fraction = np.append(Hydrogen_fraction, H_stars)
        Carbon_fraction = np.append(Carbon_fraction, C_stars)
        Silicon_fraction = np.append(Silicon_fraction, Si_stars)
        Europium_fraction = np.append(Europium_fraction, Eu_stars)
        Barium_fraction = np.append(Barium_fraction, Ba_stars)
        Strontium_fraction = np.append(Strontium_fraction, Sr_stars)
        Nitrogen_fraction = np.append(Nitrogen_fraction, N_stars)
        Neon_fraction = np.append(Neon_fraction, Ne_stars)
        Iron_fraction_SNIa = np.append(Iron_fraction_SNIa, FeSNIa_Fe)

        stars_R = np.append(stars_R, sim_info.stars.R)
        stars_z = np.append(stars_z, sim_info.stars.z)
        stars_age = np.append(stars_age, sim_info.stars.age)

    ratios = compute_ratios(Hydrogen_fraction, Magnesium_fraction, 
                            Oxygen_fraction, Iron_fraction,
                            Carbon_fraction, Silicon_fraction,
                            Europium_fraction, Barium_fraction,
                            Strontium_fraction, Nitrogen_fraction,
                            Neon_fraction, Iron_fraction_SNIa)
    
    return ratios, stars_age, stars_R, stars_z


def calculate_abundaces_from_satellite_galaxies(sim_info):

    select_mass = np.where((sim_info.halo_data.log10_stellar_mass >= 6.) &
                           (sim_info.halo_data.log10_stellar_mass <= 11.))[0]

    select_satellites = np.where(sim_info.halo_data.type[select_mass] > 10)[0]

    sample = select_mass[select_satellites]

    Oxygen_fraction = []
    Iron_fraction = []
    Magnesium_fraction = []
    Hydrogen_fraction = []
    Carbon_fraction = []
    Silicon_fraction = []
    Europium_fraction = []
    Strontium_fraction = []
    Barium_fraction = []
    Nitrogen_fraction = []
    Neon_fraction = []
    Iron_fraction_SNIa = []

    for i in tqdm(range(len(sample))):
        sim_info.make_particle_data(halo_id=sim_info.halo_data.halo_ids[sample[i]])

        O_stars = sim_info.stars.oxygen
        Fe_stars = sim_info.stars.iron
        Mg_stars = sim_info.stars.magnesium
        H_stars = sim_info.stars.hydrogen
        C_stars = sim_info.stars.carbon
        Si_stars = sim_info.stars.silicon
        Eu_stars = sim_info.stars.europium
        Ba_stars = sim_info.stars.barium
        Sr_stars = sim_info.stars.strontium
        N_stars = sim_info.stars.nitrogen
        Ne_stars = sim_info.stars.neon
        FeSNIa_Fe = sim_info.stars.iron_SNIa_fraction

        Oxygen_fraction = np.append(Oxygen_fraction, O_stars)
        Iron_fraction = np.append(Iron_fraction, Fe_stars)
        Magnesium_fraction = np.append(Magnesium_fraction, Mg_stars)
        Hydrogen_fraction = np.append(Hydrogen_fraction, H_stars)
        Carbon_fraction = np.append(Carbon_fraction, C_stars)
        Silicon_fraction = np.append(Silicon_fraction, Si_stars)
        Europium_fraction = np.append(Europium_fraction, Eu_stars)
        Barium_fraction = np.append(Barium_fraction, Ba_stars)
        Strontium_fraction = np.append(Strontium_fraction, Sr_stars)
        Nitrogen_fraction = np.append(Nitrogen_fraction, N_stars)
        Neon_fraction = np.append(Neon_fraction, Ne_stars)
        Iron_fraction_SNIa = np.append(Iron_fraction_SNIa, FeSNIa_Fe)

    ratios = compute_ratios(Hydrogen_fraction, Magnesium_fraction, 
                            Oxygen_fraction, Iron_fraction,
                            Carbon_fraction, Silicon_fraction,
                            Europium_fraction, Barium_fraction,
                            Strontium_fraction, Nitrogen_fraction,
                            Neon_fraction, Iron_fraction_SNIa)
    return ratios


def plot_stellar_abundances(sim_info, output_path, abundance_data):

    # Look for abundance ratios from COLIBRE snaps:
    ratios_MW, stars_age, stars_R, stars_z = \
        calculate_abundaces_from_MW_type_galaxies(sim_info)
    O_Fe = ratios_MW['O_Fe']
    Mg_Fe = ratios_MW['Mg_Fe']
    Fe_H = ratios_MW['Fe_H']
    FeSNIa_Fe = ratios_MW['FeSNIa_Fe']

    # Load Satellite data:
    ratios_sat = calculate_abundaces_from_satellite_galaxies(sim_info)
    Fe_H_sat = ratios_sat['Fe_H']
    O_Fe_sat = ratios_sat['O_Fe']
    Mg_Fe_sat = ratios_sat['Mg_Fe']

    # Plot parameters
    params = {
        "font.size": 11,
        "font.family": "STIXGeneral",
        "text.usetex": False,
        "mathtext.fontset": "stix",
        "figure.figsize": (7, 3),
        "figure.subplot.left": 0.1,
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

    # Box stellar abundance --------------------------------
    ax = plt.subplot(1, 2, 1)
    plt.grid("True")

    plt.plot(Fe_H, O_Fe, 'o', ms=0.2, color='grey', alpha=0.2)

    plot_MW_data('O')
    plot_GALAH_data('O')
    
    bins = np.arange(-7.2, 1, 0.2)
    ind = np.digitize(Fe_H, bins)
    xm = [np.median(Fe_H[ind == i]) for i in range(1, len(bins)) if len(Fe_H[ind == i]) > 10]
    ym = [np.median(O_Fe[ind == i]) for i in range(1, len(bins)) if len(O_Fe[ind == i]) > 10]
    plt.plot(xm, ym, '-', lw=0.5, color='tab:blue', label='GALAH DR3')
    plt.plot(xm, ym, '-', lw=1.5, color='black', label=sim_info.simulation_name)
    Fe_H_median = xm.copy()
    O_Fe_median = ym.copy()
    counter = len(xm)

    plt.text(-3.8, 1.3, "MW-type galaxies")
    plt.xlabel("[Fe/H]", labelpad=2)
    plt.ylabel("[O/Fe]", labelpad=2)
    plt.axis([-4, 1, -1, 1.5])
    plt.legend(loc=[0.0, 0.02], labelspacing=0.1, handlelength=1.5, handletextpad=0.1, frameon=False, ncol=1,
               columnspacing=0.02)

    ########################
    ax = plt.subplot(1, 2, 2)
    plt.grid("True")

    plt.plot(Fe_H_sat, O_Fe_sat, 'o', ms=0.2, color='grey')

    plot_MW_Satellites('O')
    plot_GALAH_data('O')
    
    bins = np.arange(-7.2, 1, 0.2)
    ind = np.digitize(Fe_H_sat, bins)
    xm = [np.median(Fe_H_sat[ind == i]) for i in range(1, len(bins)) if len(Fe_H_sat[ind == i]) > 10]
    ym = [np.median(O_Fe_sat[ind == i]) for i in range(1, len(bins)) if len(O_Fe_sat[ind == i]) > 10]
    plt.plot(xm, ym, '-', lw=0.5, color='tab:blue',label='GALAH DR3')
    plt.plot(xm, ym, '-', lw=1.5, color='black')

    plt.text(-3.8, 1.3, "Satellite galaxies")
    plt.xlabel("[Fe/H]", labelpad=2)
    plt.ylabel("[O/Fe]", labelpad=2)
    plt.axis([-4, 1, -1, 1.5])
    plt.legend(loc=[0.0, 0.02], labelspacing=0.1, handlelength=1.5, handletextpad=0.1, frameon=False, ncol=2,
               columnspacing=0.02)

    plt.savefig(f"{output_path}/O_Fe_" + sim_info.simulation_name + ".png", dpi=200)

    ## Let's remake the plot with APOGEE data
    plt.figure()

    # Box stellar abundance --------------------------------
    ax = plt.subplot(1, 2, 1)
    plt.grid("True")

    plt.plot(Fe_H, O_Fe, 'o', ms=0.2, color='grey', alpha=0.2)

    plot_MW_data('O')
    plot_GALAH_data('O')

    bins = np.arange(-7.2, 1, 0.2)
    ind = np.digitize(Fe_H, bins)
    xm = [np.median(Fe_H[ind == i]) for i in range(1, len(bins)) if len(Fe_H[ind == i]) > 10]
    ym = [np.median(O_Fe[ind == i]) for i in range(1, len(bins)) if len(O_Fe[ind == i]) > 10]
    plt.plot(xm, ym, '-', lw=0.5, color='tab:blue', label='GALAH DR3')
    plt.plot(xm, ym, '-', lw=1.5, color='black', label=sim_info.simulation_name)
    Fe_H_median = xm.copy()
    O_Fe_median = ym.copy()
    counter = len(xm)

    plt.text(-3.8, 1.3, "MW-type galaxies")
    plt.xlabel("[Fe/H]", labelpad=2)
    plt.ylabel("[O/Fe]", labelpad=2)
    plt.axis([-4, 1, -1, 1.5])
    plt.legend(loc=[0.0, 0.02], labelspacing=0.1, handlelength=1.5, handletextpad=0.1, frameon=False, ncol=1,
               columnspacing=0.02)

    ########################
    ax = plt.subplot(1, 2, 2)
    plt.grid("True")

    plt.plot(Fe_H, O_Fe, 'o', ms=0.2, color='grey', alpha=0.2)

    plot_MW_data('O')
    plot_APOGEE_data('O')

    bins = np.arange(-7.2, 1, 0.2)
    ind = np.digitize(Fe_H, bins)
    xm = [np.median(Fe_H[ind == i]) for i in range(1, len(bins)) if len(Fe_H[ind == i]) > 10]
    ym = [np.median(O_Fe[ind == i]) for i in range(1, len(bins)) if len(O_Fe[ind == i]) > 10]
    plt.plot(xm, ym, '-', lw=0.5, color='blue', label='APOGEE data')
    plt.plot(xm, ym, '-', lw=1.5, color='black', label=sim_info.simulation_name)

    plt.text(-3.8, 1.3, "MW-type galaxies")
    plt.xlabel("[Fe/H]", labelpad=2)
    plt.ylabel("[O/Fe]", labelpad=2)
    plt.axis([-4, 1, -1, 1.5])
    plt.legend(loc=[0.0, 0.02], labelspacing=0.1, handlelength=1.5, handletextpad=0.1, frameon=False, ncol=1,
               columnspacing=0.02)

    plt.savefig(f"{output_path}/O_Fe_" + sim_info.simulation_name + "_2.png", dpi=200)

    # New plot. Now turn Magnesium / Fe --------------------------------

    # Box stellar abundance --------------------------------
    fig = plt.figure()
    ax = plt.subplot(1, 2, 1)
    plt.grid("True")

    plt.plot(Fe_H, Mg_Fe, 'o', ms=0.5, color='grey',alpha=0.5)

    plot_MW_data('Mg')
    plot_GALAH_data('Mg')
    
    bins = np.arange(-7.2, 1, 0.2)
    ind = np.digitize(Fe_H, bins)
    xm = [np.median(Fe_H[ind == i]) for i in range(1, len(bins)) if len(Fe_H[ind == i]) > 10]
    ym = [np.median(Mg_Fe[ind == i]) for i in range(1, len(bins)) if len(Mg_Fe[ind == i]) > 10]
    plt.plot(xm, ym, '-', lw=0.5, color='tab:blue', label='GALAH DR3')
    plt.plot(xm, ym, '-', lw=1.5, color='black',label=sim_info.simulation_name)
    Mg_Fe_median = ym.copy()

    plt.xlabel("[Fe/H]", labelpad=2)
    plt.ylabel("[Mg/Fe]", labelpad=2)
    plt.text(-3.8, 1.2, "MW-type galaxies")
    plt.axis([-4, 1, -2, 1.5])
    plt.legend(loc=[0, 0.02], labelspacing=0.1, handlelength=1.5, handletextpad=0.1, frameon=False, ncol=1,
               columnspacing=0.02)

    ###############
    ax = plt.subplot(1, 2, 2)
    plt.grid("True")

    plt.plot(Fe_H_sat, Mg_Fe_sat, 'o', ms=0.5, color='grey')

    plot_MW_Satellites('Mg')
    plot_GALAH_data('Mg')
    
    bins = np.arange(-7.2, 1, 0.2)
    ind = np.digitize(Fe_H_sat, bins)
    xm = [np.median(Fe_H_sat[ind == i]) for i in range(1, len(bins)) if len(Fe_H_sat[ind == i]) > 10]
    ym = [np.median(Mg_Fe_sat[ind == i]) for i in range(1, len(bins)) if len(Mg_Fe_sat[ind == i]) > 10]
    plt.plot(xm, ym, '-', lw=0.5, color='tab:blue', label='GALAH DR3')
    plt.plot(xm, ym, '-', lw=1.5, color='black')

    plt.text(-3.8,1.2,"Satellite galaxies")
    plt.xlabel("[Fe/H]", labelpad=2)
    plt.ylabel("[Mg/Fe]", labelpad=2)
    plt.axis([-4, 1, -2, 1.5])
    plt.legend(loc=[0, 0.02], labelspacing=0.1, handlelength=1.5, handletextpad=0.1, frameon=False, ncol=2,
               columnspacing=0.02)

    plt.savefig(f"{output_path}/Mg_Fe_" + sim_info.simulation_name + ".png", dpi=200)

    # Box stellar abundance --------------------------------
    fig = plt.figure()
    ax = plt.subplot(1, 2, 1)
    plt.grid("True")

    plt.plot(Fe_H, Mg_Fe, 'o', ms=0.5, color='grey', alpha=0.5)

    plot_MW_data('Mg')
    plot_GALAH_data('Mg')

    bins = np.arange(-7.2, 1, 0.2)
    ind = np.digitize(Fe_H, bins)
    xm = [np.median(Fe_H[ind == i]) for i in range(1, len(bins)) if len(Fe_H[ind == i]) > 10]
    ym = [np.median(Mg_Fe[ind == i]) for i in range(1, len(bins)) if len(Mg_Fe[ind == i]) > 10]
    plt.plot(xm, ym, '-', lw=0.5, color='tab:blue', label='GALAH DR3')
    plt.plot(xm, ym, '-', lw=1.5, color='black', label=sim_info.simulation_name)

    plt.xlabel("[Fe/H]", labelpad=2)
    plt.ylabel("[Mg/Fe]", labelpad=2)
    plt.text(-3.8, 1.2, "MW-type galaxies")
    plt.axis([-4, 1, -2, 1.5])
    plt.legend(loc=[0, 0.02], labelspacing=0.1, handlelength=1.5, handletextpad=0.1, frameon=False, ncol=1,
               columnspacing=0.02)

    ###############
    ax = plt.subplot(1, 2, 2)
    plt.grid("True")

    plt.plot(Fe_H, Mg_Fe, 'o', ms=0.5, color='grey', alpha=0.5)

    plot_MW_data('Mg')
    plot_APOGEE_data('MG')

    bins = np.arange(-7.2, 1, 0.2)
    ind = np.digitize(Fe_H, bins)
    xm = [np.median(Fe_H[ind == i]) for i in range(1, len(bins)) if len(Fe_H[ind == i]) > 10]
    ym = [np.median(Mg_Fe[ind == i]) for i in range(1, len(bins)) if len(Mg_Fe[ind == i]) > 10]
    plt.plot(xm, ym, '-', lw=0.5, color='blue', label='APOGEE data')
    plt.plot(xm, ym, '-', lw=1.5, color='black', label=sim_info.simulation_name)

    plt.xlabel("[Fe/H]", labelpad=2)
    plt.ylabel("[Mg/Fe]", labelpad=2)
    plt.text(-3.8, 1.2, "MW-type galaxies")
    plt.axis([-4, 1, -2, 1.5])
    plt.legend(loc=[0, 0.02], labelspacing=0.1, handlelength=1.5, handletextpad=0.1, frameon=False, ncol=1,
               columnspacing=0.02)

    plt.savefig(f"{output_path}/Mg_Fe_" + sim_info.simulation_name + "_2.png", dpi=200)

    # make remaining plots with just GALAH data (Buder+21)
    for el in ['C','Si','Eu','Ba']:
        fig = plt.figure(figsize=(3.8,3))
        ax = plt.subplot(1, 1, 1)
        plt.grid("True")

        plt.plot(Fe_H, ratios_MW[f'{el}_Fe'], 'o', ms=0.5, color='grey',alpha=0.5)
        plot_GALAH_data(el)

        bins = np.arange(-7.2, 1, 0.2)
        ind = np.digitize(Fe_H, bins)
        xm = [np.median(Fe_H[ind == i]) for i in range(1, len(bins)) if len(Fe_H[ind == i]) > 10]
        ym = [np.median(ratios_MW[f'{el}_Fe'][ind == i]) for i in range(1, len(bins)) if len(ratios_MW[f'{el}_Fe'][ind == i]) > 10]
        plt.plot(xm, ym, '-', lw=0.5, color='tab:blue', label='GALAH DR3')
        plt.plot(xm, ym, '-', lw=1.5, color='black',label=sim_info.simulation_name)
        if el == 'C': C_Fe_median = ym.copy()
        if el == 'Si': Si_Fe_median = ym.copy()
        if el == 'Eu': Eu_Fe_median = ym.copy()
        if el == 'Ba': Ba_Fe_median = ym.copy()

        plt.xlabel("[Fe/H]", labelpad=2)
        plt.ylabel(f"[{el}/Fe]", labelpad=2)
        plt.text(-3.8, 1.2, "MW-type galaxies")
        plt.axis([-4, 1, -2, 1.5])
        plt.legend(loc=[0, 0.02], labelspacing=0.1, handlelength=1.5, handletextpad=0.1, frameon=False, ncol=1,
                   columnspacing=0.02)
        plt.tight_layout()
        plt.savefig(f"{output_path}/{el}_Fe_" + sim_info.simulation_name + ".png", dpi=200)

    ## Here an additional plot for Carbon
    fig = plt.figure(figsize=(3.8, 3))
    ax = plt.subplot(1, 1, 1)
    plt.grid("True")

    plt.plot(Fe_H, ratios_MW['C_Fe'], 'o', ms=0.5, color='grey', alpha=0.5)
    plot_APOGEE_data('C')

    bins = np.arange(-7.2, 1, 0.2)
    ind = np.digitize(Fe_H, bins)
    xm = [np.median(Fe_H[ind == i]) for i in range(1, len(bins)) if len(Fe_H[ind == i]) > 10]
    ym = [np.median(ratios_MW['C_Fe'][ind == i]) for i in range(1, len(bins)) if
          len(ratios_MW[f'{el}_Fe'][ind == i]) > 10]
    plt.plot(xm, ym, '-', lw=0.5, color='blue', label='APOGEE data')
    plt.plot(xm, ym, '-', lw=1.5, color='black', label=sim_info.simulation_name)

    plt.xlabel("[Fe/H]", labelpad=2)
    plt.ylabel("[C/Fe]", labelpad=2)
    plt.text(-3.8, 1.2, "MW-type galaxies")
    plt.axis([-4, 1, -2, 1.5])
    plt.legend(loc=[0, 0.02], labelspacing=0.1, handlelength=1.5, handletextpad=0.1, frameon=False, ncol=1,
               columnspacing=0.02)
    plt.tight_layout()
    plt.savefig(f"{output_path}/C_Fe_" + sim_info.simulation_name + "_2.png", dpi=200)

    # make remaining plots with just GALAH data (Buder+21)
    for el in ['N', 'Ne']:
        fig = plt.figure(figsize=(3.8, 3))
        ax = plt.subplot(1, 1, 1)
        plt.grid("True")

        plt.plot(Fe_H, ratios_MW[f'{el}_Fe'], 'o', ms=0.5, color='grey', alpha=0.5)

        bins = np.arange(-7.2, 1, 0.2)
        ind = np.digitize(Fe_H, bins)
        xm = [np.median(Fe_H[ind == i]) for i in range(1, len(bins)) if len(Fe_H[ind == i]) > 10]
        ym = [np.median(ratios_MW[f'{el}_Fe'][ind == i]) for i in range(1, len(bins)) if
              len(ratios_MW[f'{el}_Fe'][ind == i]) > 10]

        if el == 'N':
            plot_APOGEE_data('N')
            plt.plot(xm, ym, '-', lw=0.5, color='blue', label='APOGEE data')
            N_Fe_median = ym.copy()
        if el == 'Ne': Ne_Fe_median = ym.copy()

        plt.plot(xm, ym, '-', lw=1.5, color='black', label=sim_info.simulation_name)
        plt.xlabel("[Fe/H]", labelpad=2)
        plt.ylabel(f"[{el}/Fe]", labelpad=2)
        plt.text(-3.8, 1.2, "MW-type galaxies")
        plt.axis([-4, 1, -2, 1.5])
        plt.legend(loc=[0, 0.02], labelspacing=0.1, handlelength=1.5, handletextpad=0.1, frameon=False, ncol=1,
                   columnspacing=0.02)
        plt.tight_layout()
        plt.savefig(f"{output_path}/{el}_Fe_" + sim_info.simulation_name + ".png", dpi=200)

    ######
    fig = plt.figure(figsize=(3.8, 3))
    ax = plt.subplot(1, 1, 1)
    plt.grid("True")

    plt.plot(Fe_H, ratios_MW['Sr_Fe'], 'o', ms=0.5, color='grey',alpha=0.2)
    plot_StrontiumObsData()

    bins = np.arange(-7.2, 1, 0.2)
    ind = np.digitize(Fe_H, bins)
    xm = [np.median(Fe_H[ind == i]) for i in range(1, len(bins)) if len(Fe_H[ind == i]) > 10]
    ym = [np.median(ratios_MW['Sr_Fe'][ind == i]) for i in range(1, len(bins)) if
          len(ratios_MW['Sr_Fe'][ind == i]) > 10]
    plt.plot(xm, ym, '-', lw=1.5, color='black', label=sim_info.simulation_name)
    Sr_Fe_median = ym.copy()

    plt.xlabel("[Fe/H]", labelpad=2)
    plt.ylabel("[Sr/Fe]", labelpad=2)
    plt.text(-3.8, 1.2, "MW-type galaxies")
    plt.axis([-4, 1, -2, 1.5])
    plt.legend(loc=[0, 0.02], labelspacing=0.1, handlelength=1.5, handletextpad=0.1, frameon=False, ncol=1,
               columnspacing=0.02)
    plt.tight_layout()
    plt.savefig(f"{output_path}/Sr_Fe_" + sim_info.simulation_name + ".png", dpi=200)

    ########
    fig = plt.figure(figsize=(3.8, 3))
    ax = plt.subplot(1, 1, 1)
    plt.grid("True")

    plt.plot(Fe_H, FeSNIa_Fe, 'o', ms=0.5, color='grey', alpha=0.5)

    bins = np.arange(-7.2, 1, 0.2)
    ind = np.digitize(Fe_H, bins)
    xm = [np.median(Fe_H[ind == i]) for i in range(1, len(bins)) if len(Fe_H[ind == i]) > 10]
    ym = [np.median(FeSNIa_Fe[ind == i]) for i in range(1, len(bins)) if len(FeSNIa_Fe[ind == i]) > 10]
    FeSNIa_Fe_median = ym.copy()
    plt.plot(xm, ym, '-', lw=1.5, color='black', label=sim_info.simulation_name)

    plt.xlabel("[Fe/H]", labelpad=2)
    plt.ylabel("Fe(SNIa)/Fe", labelpad=2)
    plt.yscale('log')
    plt.axis([-4, 1, 1e-2,2])
    ax.tick_params(direction='in', axis='both', which='both')
    plt.legend(loc=[0, 0.02], labelspacing=0.1, handlelength=1.5, handletextpad=0.1, frameon=False, ncol=1,
               columnspacing=0.02)
    plt.tight_layout()
    plt.savefig(f"{output_path}/FeSNIa_Fe_" + sim_info.simulation_name + ".png", dpi=200)

    if abundance_data == None:
        abundance_data = {'Fe_H': Fe_H_median,
                          'O_Fe': O_Fe_median,
                          'Mg_Fe': Mg_Fe_median,
                          'C_Fe': C_Fe_median,
                          'Ba_Fe': Ba_Fe_median,
                          'Sr_Fe': Sr_Fe_median,
                          'Si_Fe': Si_Fe_median,
                          'Eu_Fe': Eu_Fe_median,
                          'N_Fe': N_Fe_median,
                          'Ne_Fe': Ne_Fe_median,
                          'FeSNIa_Fe':FeSNIa_Fe_median,
                          'counter':counter}
    else:
        Fe_H = abundance_data['Fe_H']
        Fe_H = np.append(Fe_H,Fe_H_median)
        O_Fe = abundance_data['O_Fe']
        O_Fe = np.append(O_Fe,O_Fe_median)
        Mg_Fe = abundance_data['Mg_Fe']
        Mg_Fe = np.append(Mg_Fe,Mg_Fe_median)
        C_Fe = abundance_data['C_Fe']
        C_Fe = np.append(C_Fe, C_Fe_median)
        Ba_Fe = abundance_data['Ba_Fe']
        Ba_Fe = np.append(Ba_Fe, Ba_Fe_median)
        Sr_Fe = abundance_data['Sr_Fe']
        Sr_Fe = np.append(Sr_Fe, Sr_Fe_median)
        Si_Fe = abundance_data['Si_Fe']
        Si_Fe = np.append(Si_Fe, Si_Fe_median)
        Eu_Fe = abundance_data['Eu_Fe']
        Eu_Fe = np.append(Eu_Fe, Eu_Fe_median)
        N_Fe = abundance_data['N_Fe']
        N_Fe = np.append(N_Fe, N_Fe_median)
        Ne_Fe = abundance_data['Ne_Fe']
        Ne_Fe = np.append(Ne_Fe, Ne_Fe_median)
        FeSNIa_Fe = abundance_data['FeSNIa_Fe']
        FeSNIa_Fe = np.append(FeSNIa_Fe,FeSNIa_Fe_median)

        counter_sim = abundance_data['counter']
        counter_sim = np.append(counter_sim, counter)
        abundance_data = {'Fe_H': Fe_H,
                          'O_Fe': O_Fe,
                          'Mg_Fe': Mg_Fe,
                          'C_Fe': C_Fe,
                          'Ba_Fe': Ba_Fe,
                          'Sr_Fe': Sr_Fe,
                          'Eu_Fe': Eu_Fe,
                          'Si_Fe': Si_Fe,
                          'N_Fe': N_Fe,
                          'Ne_Fe': Ne_Fe,
                          'FeSNIa_Fe': FeSNIa_Fe,
                          'counter': counter_sim}


    # Additionally , let's try a few things #######

    O_Fe = ratios_MW['O_Fe']
    Mg_Fe = ratios_MW['Mg_Fe']
    Fe_H = ratios_MW['Fe_H']

    # Plot parameters
    params = {
        "font.size": 11,
        "font.family": "STIXGeneral",
        "text.usetex": False,
        "mathtext.fontset": "stix",
        "figure.figsize": (7, 3),
        "figure.subplot.left": 0.1,
        "figure.subplot.right": 0.95,
        "figure.subplot.bottom": 0.15,
        "figure.subplot.top": 0.83,
        "figure.subplot.wspace": 0.3,
        "figure.subplot.hspace": 0.3,
        "lines.markersize": 0.5,
        "lines.linewidth": 0.2,
    }
    rcParams.update(params)

    plt.figure()
    cm = matplotlib.cm.get_cmap('RdYlBu')

    # Box stellar abundance --------------------------------
    ax = plt.subplot(1, 2, 1)
    plt.grid("True")

    sort_age = np.argsort(stars_age)
    sort_age = sort_age[::-1]

    plt.scatter(Fe_H[sort_age], O_Fe[sort_age], c=stars_age[sort_age], s=0.2, vmin=0, vmax=14, cmap=cm)

    # plot_MW_data('O')
    # plot_GALAH_data('O')

    plt.text(-3.8, 1.3, "MW-type galaxies")
    plt.xlabel("[Fe/H]", labelpad=2)
    plt.ylabel("[O/Fe]", labelpad=2)
    plt.axis([-4, 1, -1, 1.5])
    axins = inset_axes(ax, width="100%", height="4%", loc='upper center', borderpad=-2.9)
    cb = plt.colorbar(cax=axins, orientation="horizontal")
    cb.set_ticks([0, 2, 4, 6, 8, 10, 12, 14])
    cb.set_label('Stellar Age [Gyr]', labelpad=-2)

    #####
    ax = plt.subplot(1, 2, 2)
    plt.grid("True")

    inner = np.where((stars_R > 9) & (stars_z < 2))[0]
    plt.plot(Fe_H[inner], O_Fe[inner], 'o', ms=0.5, color='tab:blue', alpha=1,
             label='$R > 9~\mathrm{kpc}, |z|<2~\mathrm{kpc}$')

    inner = np.where((stars_R <= 7) & (np.abs(stars_z) < 2))[0]
    plt.plot(Fe_H[inner], O_Fe[inner], 'o', ms=0.5, color='tab:red', alpha=1,
             label='$R < 7~\mathrm{kpc}, |z|<2~\mathrm{kpc}$')

    inner = np.where((stars_R >= 7) & (stars_R <= 9) & (np.abs(stars_z) < 2))[0]
    plt.plot(Fe_H[inner], O_Fe[inner], 'o', ms=0.5, color='tab:orange', alpha=1,
             label = '$7~\mathrm{kpc} < R < 9~\mathrm{kpc}, |z|<2~\mathrm{kpc}$')

    # plot_MW_data('O')
    # plot_GALAH_data('O')

    plt.text(-3.8, 1.3, "MW-type galaxies")
    plt.xlabel("[Fe/H]", labelpad=2)
    plt.ylabel("[O/Fe]", labelpad=2)
    plt.axis([-4, 1, -1, 1.5])
    plt.legend(loc=[0.0, 0.97], labelspacing=0.1, handlelength=1.5, handletextpad=0.1, frameon=False, ncol=1,
               columnspacing=0.02, fontsize=9)

    plt.savefig(f"{output_path}/O_Fe_" + sim_info.simulation_name + "_test.png", dpi=300)

    #####
    plt.figure()
    cm = matplotlib.cm.get_cmap('RdYlBu')

    # Box stellar abundance --------------------------------
    ax = plt.subplot(1, 2, 1)
    plt.grid("True")

    sort_age = np.argsort(stars_age)
    sort_age = sort_age[::-1]

    plt.scatter(Fe_H[sort_age], Mg_Fe[sort_age], c=stars_age[sort_age], s=0.2, vmin=0, vmax=14, cmap=cm)

    # plot_MW_data('O')
    # plot_GALAH_data('O')

    plt.text(-3.8, 1.2, "MW-type galaxies")
    plt.xlabel("[Fe/H]", labelpad=2)
    plt.ylabel("[Mg/Fe]", labelpad=2)
    plt.axis([-4, 1, -2, 1.5])
    axins = inset_axes(ax, width="100%", height="4%", loc='upper center', borderpad=-2.9)
    cb = plt.colorbar(cax=axins, orientation="horizontal")
    cb.set_ticks([0, 2, 4, 6, 8, 10, 12, 14])
    cb.set_label('Stellar Age [Gyr]', labelpad=-2)

    #####
    ax = plt.subplot(1, 2, 2)
    plt.grid("True")

    inner = np.where((stars_R > 9) & (stars_z < 2))[0]
    plt.plot(Fe_H[inner], Mg_Fe[inner], 'o', ms=0.5, color='tab:blue', alpha=1,
             label='$R > 9~\mathrm{kpc}, |z|<2~\mathrm{kpc}$')

    inner = np.where((stars_R <= 7) & (np.abs(stars_z) < 2))[0]
    plt.plot(Fe_H[inner], Mg_Fe[inner], 'o', ms=0.5, color='tab:red', alpha=1,
             label='$R < 7~\mathrm{kpc}, |z|<2~\mathrm{kpc}$')

    inner = np.where((stars_R >= 7) & (stars_R <= 9) & (np.abs(stars_z) < 2))[0]
    plt.plot(Fe_H[inner], Mg_Fe[inner], 'o', ms=0.5, color='tab:orange', alpha=1,
             label = '$7~\mathrm{kpc} < R < 9~\mathrm{kpc}, |z|<2~\mathrm{kpc}$')

    # plot_MW_data('O')
    # plot_GALAH_data('O')

    plt.text(-3.8, 1.2, "MW-type galaxies")
    plt.xlabel("[Fe/H]", labelpad=2)
    plt.ylabel("[Mg/Fe]", labelpad=2)
    plt.axis([-4, 1, -2, 1.5])
    plt.legend(loc=[0.0, 0.97], labelspacing=0.1, handlelength=1.5, handletextpad=0.1, frameon=False, ncol=1,
               columnspacing=0.02, fontsize=9)

    plt.savefig(f"{output_path}/Mg_Fe_" + sim_info.simulation_name + "_test.png", dpi=300)

    #####
    plt.figure()
    cm = matplotlib.cm.get_cmap('RdYlBu')

    # Box stellar abundance --------------------------------
    ax = plt.subplot(1, 2, 1)
    plt.grid("True")

    select = np.where(stars_age < 8)[0]
    sort_age = np.argsort(stars_age[select])
    sort_age = sort_age[::-1]
    plt.scatter(Fe_H[select[sort_age]], Mg_Fe[select[sort_age]], c=stars_age[select[sort_age]],
                s=0.2, vmin=0, vmax=8, cmap=cm)

    # plot_MW_data('O')
    # plot_GALAH_data('O')

    plt.text(-3.8, 1.2, "MW-type galaxies")
    plt.xlabel("[Fe/H]", labelpad=2)
    plt.ylabel("[Mg/Fe]", labelpad=2)
    plt.axis([-4, 1, -2, 1.5])
    axins = inset_axes(ax, width="100%", height="4%", loc='upper center', borderpad=-2.9)
    cb = plt.colorbar(cax=axins, orientation="horizontal")
    cb.set_ticks([0, 2, 4, 6, 8])
    cb.set_label('Stellar Age [Gyr]', labelpad=-2)

    #####
    ax = plt.subplot(1, 2, 2)
    plt.grid("True")

    select = np.where((stars_age > 4))[0]
    plt.plot(Fe_H[select], Mg_Fe[select], 'o', ms=0.5, color='crimson', alpha=1,
             label='Stellar age $>$ 4 Gyr')

    select = np.where((stars_age > 2) & (stars_age < 4))[0]
    plt.plot(Fe_H[select], Mg_Fe[select], 'o', ms=0.5, color='tab:orange', alpha=1,
             label='2 Gyr $<$ Stellar age $<$ 4 Gyr')

    select = np.where(stars_age < 2)[0]
    plt.plot(Fe_H[select], Mg_Fe[select], 'o', ms=0.5, color='tab:blue', alpha=1,
             label='Stellar age $<$ 2 Gyr')

    # plot_MW_data('O')
    # plot_GALAH_data('O')

    plt.text(-3.8, 1.2, "MW-type galaxies")
    plt.xlabel("[Fe/H]", labelpad=2)
    plt.ylabel("[Mg/Fe]", labelpad=2)
    plt.axis([-4, 1, -2, 1.5])
    plt.legend(loc=[0.0, 0.97], labelspacing=0.1, handlelength=1.5, handletextpad=0.1, frameon=False, ncol=1,
               columnspacing=0.02, fontsize=9)

    plt.savefig(f"{output_path}/Mg_Fe_" + sim_info.simulation_name + "_test_2.png", dpi=300)

    return abundance_data

