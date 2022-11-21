import matplotlib.pylab as plt
from matplotlib.pylab import rcParams
import numpy as np
from tqdm import tqdm
from .morphology import calculate_kappa_co

def compute_abundances(X_O, X_C, X_N, X_H, weight):

    mp_in_cgs = 1.6737236e-24
    mH_in_cgs = 1.00784 * mp_in_cgs
    mO_in_cgs = 15.999 * mp_in_cgs
    mC_in_cgs =  12.0107 * mp_in_cgs
    mN_in_cgs = 14.0067 * mp_in_cgs

    ratio_NO = np.sum( weight * X_N ) / np.sum( weight * X_O )
    ratio_NO = np.log10(ratio_NO) - np.log10(mN_in_cgs / mO_in_cgs)

    ratio_OH = np.sum( weight * X_O ) / np.sum( weight * X_H )
    ratio_OH = 12 + np.log10(ratio_OH) - np.log10(mO_in_cgs / mH_in_cgs)

    ratio_CO = np.sum( weight * X_C ) / np.sum( weight * X_O )
    ratio_CO = np.log10(ratio_CO) - np.log10(mC_in_cgs / mO_in_cgs)

    return {'NO': ratio_NO, 'CO': ratio_CO, 'OH':ratio_OH}

def calculate_galactic_abundances(sim_info):

    select_sample = np.where((sim_info.halo_data.log10_stellar_mass >= 7.0) &
                              (sim_info.halo_data.log10_stellar_mass <= 12))[0]

    select_centrals = np.where(sim_info.halo_data.type[select_sample] == 10)[0]

    sample = select_sample[select_centrals]
    num_sample = len(sample)

    OH = np.zeros(num_sample)
    NO = np.zeros(num_sample)
    CO = np.zeros(num_sample)

    stellar_mass = sim_info.halo_data.log10_stellar_mass[sample]
    SFR = sim_info.halo_data.sfr[sample]
    kappa = np.zeros(num_sample)


    for i in tqdm(range(len(sample))):
        sim_info.make_particle_data(halo_id=sim_info.halo_data.halo_ids[sample[i]])

        # M_gas = sim_info.gas.mass.copy()
        # H_gas = sim_info.gas.hydrogen.copy()
        # O_gas = sim_info.gas.oxygen.copy()
        # N_gas = sim_info.gas.nitrogen.copy()
        # C_gas = sim_info.gas.carbon.copy()
        #
        # ratios = compute_abundances(O_gas, C_gas, N_gas, H_gas, M_gas)

        M_stars = sim_info.stars.mass.copy()
        H_stars = sim_info.stars.hydrogen.copy()
        O_stars = sim_info.stars.oxygen.copy()
        N_stars = sim_info.stars.nitrogen.copy()
        C_stars = sim_info.stars.carbon.copy()

        ratios = compute_abundances(O_stars, C_stars, N_stars, H_stars, M_stars)

        OH[i] = ratios['OH']
        NO[i] = ratios['NO']
        CO[i] = ratios['CO']

        kappa[i] = calculate_kappa_co(sim_info, sample[i])

    return{
        'OH': OH, 'CO': CO, 'NO': NO, 'Mstellar': stellar_mass,
        'StarFormationRate': SFR, 'kappa': kappa
    }

def compute_galactic_abundances(sim_info, galactic_abundance_data):

    # Look for abundance ratios from COLIBRE snaps:
    data = calculate_galactic_abundances(sim_info)

    Mstellar = data['Mstellar']
    SFR = data['StarFormationRate']
    OH = data['OH']
    CO = data['CO']
    NO = data['NO']
    kappa = data['kappa']
    counter = np.array([len(Mstellar)])

    if galactic_abundance_data == None:

        galactic_abundance_data = {
            'Mstellar': Mstellar,
            'StarFormationRate': SFR,
            'OH': OH,
            'CO': CO,
            'NO': NO,
            'kappa': kappa,
            'counter':counter}

    else:

        Mstellar_prev = galactic_abundance_data['Mstellar']
        Mstellar_prev = np.append(Mstellar_prev, Mstellar)
        counter_prev = galactic_abundance_data['counter']
        counter_prev = np.append(counter_prev, counter)

        SFR_prev = galactic_abundance_data['StarFormationRate']
        SFR_prev = np.append(SFR_prev, SFR)

        kappa_prev = galactic_abundance_data['kappa']
        kappa_prev = np.append(kappa_prev, kappa)

        OH_prev = galactic_abundance_data['OH']
        OH_prev = np.append(OH_prev, OH)

        CO_prev = galactic_abundance_data['CO']
        CO_prev = np.append(CO_prev, CO)

        NO_prev = galactic_abundance_data['NO']
        NO_prev = np.append(NO_prev, NO)

        galactic_abundance_data = {
            'StarFormationRate': SFR_prev,
            'OH':OH_prev,
            'CO':CO_prev,
            'NO':NO_prev,
            'Mstellar': Mstellar_prev,
            'kappa': kappa_prev,
            'counter': counter_prev
        }

    return galactic_abundance_data

def plot_berg_2019(ref):

    file = './plotter/obs_data/Berg_2019.txt'
    data = np.loadtxt(file)
    OH = data[:,3]
    if ref == 'CO':
        y = data[:,5]
    if ref == 'NO':
        y = data[:,7]

    plt.plot(OH, y, 'v', ms=2, color='black', label='Berg et al. (2019) (Dwarf galaxies)')

def plot_israelian_2004():

    N_H_Sun = 8.05
    O_H_Sun = 8.93
    N_O_Sun = N_H_Sun - O_H_Sun

    # file = './plotter/obs_data/Israelian_2004.txt'
    # data = np.loadtxt(file)
    # OH = data[:,0] + O_H_Sun
    # NO = data[:,2] + N_O_Sun

    # Define the scatter as offset from the mean value
    # x_scatter = np.array((
    #     data[:,0] - data[:,1]/2 + O_H_Sun,
    #     data[:,0] + data[:,1]/2 + O_H_Sun
    # ))
    # plt.errorbar(OH, NO, xerr= x_scatter, fmt='*',
    #              ls='none',color='lightgreen', ms=2, label='Israelian et al. (2004)')
    # plt.plot(OH, NO, '*', ms=2, color='lightgreen', label='Israelian et al. (2004) (metal-rich stars)')

    file = './plotter/obs_data/Israelian_2004_NH.txt'
    data = np.loadtxt(file)
    OH = data[:,0] + O_H_Sun
    NH = data[:,1]
    NO = NH - data[:,0] + N_O_Sun

    plt.plot(OH, NO, '*', ms=2, color='lightgreen', label='Israelian et al. (2004) (metal-rich stars)')



def plot_CO_OH_relation(CO, OH, Mstellar, kappa, counter, output_name_list, output_file):

    # Plot parameters
    params = {
        "font.size": 12,
        "font.family": "STIXGeneral",
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
    plt.grid("True")

    # Nicholls+ 2017 best-fitting relation
    x = np.arange(7, 9.6, 0.1)
    y = np.log10(10**-0.8 + 10**(x-12+2.72))
    plt.plot(x, y, '--', lw=1.5, color='purple', label='Nicholls et al. (2017)')

    plot_berg_2019('CO')

    count = 0
    color = ['tab:blue','tab:green','tab:orange','crimson','tab:purple']
    for i in range(len(output_name_list)):
        xm = OH[count:count + counter[i]]
        ym = CO[count:count + counter[i]]

        kpm = kappa[count:count + counter[i]]
        mm = Mstellar[count:count + counter[i]]

        count += counter[i]

        # bins = np.arange(7, 9, 0.2)
        # ind = np.digitize(xm, bins)
        # ylo = [np.percentile(10**ym[ind == i],16) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
        # yhi = [np.percentile(10**ym[ind == i],84) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
        # yl = [np.median(10**ym[ind == i]) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
        # xl = [np.median(10**xm[ind == i]) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
        # plt.fill_between(xl, ylo, yhi, color=color[i], alpha=0.2)
        # plt.plot(xl, yl, '-', lw=1.5, color=color[i], label=output_name_list[i])
        # plt.plot(xm, ym, 'o', ms=1.5, color=color[i])

        plt.plot(xm[kpm < 0.4], ym[kpm < 0.4], 'o', ms=2, color=color[i],
                 label=output_name_list[i] + ' (Elliptical galaxies)')
        plt.plot(xm[kpm >= 0.4], ym[kpm >= 0.4], '+', ms=5, color=color[i],
                 label=output_name_list[i] + ' (Disc galaxies)')

    plt.ylabel(r'$\log_{10}$(C/O)', labelpad=2)
    plt.xlabel(r'$12+\log_{10}$(O/H)', labelpad=2)
    ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    plt.axis([7, 9.5, -1.5, 0.5])
    plt.legend(loc='upper left', labelspacing=0.1, handlelength=1.5, handletextpad=0.1, frameon=False, ncol=1,
               fontsize=9, columnspacing=0.02)

    plt.savefig(output_file, dpi=300)

def plot_NO_OH_relation(NO, OH, Mstellar, kappa, counter, output_name_list, output_file):

    # Plot parameters
    params = {
        "font.size": 12,
        "font.family": "STIXGeneral",
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
    plt.grid("True")

    # Nicholls+ 2017 best-fitting relation
    x = np.arange(7, 9.6, 0.1)
    y = np.log10(10**-1.732 + 10**(x-12+2.19))
    plt.plot(x, y, '--', lw=1.5, color='purple', label='Nicholls et al. (2017)')

    plot_berg_2019('NO')
    plot_israelian_2004()

    count = 0
    color = ['tab:blue','tab:green','tab:orange','crimson','tab:purple']
    for i in range(len(output_name_list)):
        xm = OH[count:count + counter[i]]
        ym = NO[count:count + counter[i]]

        kpm = kappa[count:count + counter[i]]
        mm = Mstellar[count:count + counter[i]]

        count += counter[i]

        # bins = np.arange(7, 9, 0.2)
        # ind = np.digitize(xm, bins)
        # ylo = [np.percentile(10**ym[ind == i],16) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
        # yhi = [np.percentile(10**ym[ind == i],84) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
        # yl = [np.median(10**ym[ind == i]) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
        # xl = [np.median(10**xm[ind == i]) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
        # plt.fill_between(xl, ylo, yhi, color=color[i], alpha=0.2)
        # plt.plot(xl, yl, '-', lw=1.5, color=color[i], label=output_name_list[i])

        plt.plot(xm[kpm < 0.4], ym[kpm < 0.4], 'o', ms=2, color=color[i],
                 label=output_name_list[i] + ' (Elliptical galaxies)')
        plt.plot(xm[kpm >= 0.4], ym[kpm >= 0.4], '+', ms=5, color=color[i],
                 label=output_name_list[i] + ' (Disc galaxies)')

    plt.ylabel(r'$\log_{10}$(N/O)', labelpad=2)
    plt.xlabel(r'$12+\log_{10}$(O/H)', labelpad=2)
    ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    plt.axis([7, 9.5, -2, 0])
    plt.legend(loc='upper left', labelspacing=0.1, handlelength=1.5, handletextpad=0.1, frameon=False, ncol=1,
               fontsize=9, columnspacing=0.02)

    plt.savefig(output_file, dpi=300)


def plot_galactic_abundance_relations(sim_data, output_name_list, output_path):

    Mstellar = sim_data['Mstellar']
    counter = sim_data['counter']
    SFR = sim_data['StarFormationRate']
    CO = sim_data['CO']
    OH = sim_data['OH']
    NO = sim_data['NO']
    kappa = sim_data['kappa']

    output_file = f"{output_path}/Galactic_CO_OH_comparison.png"
    plot_CO_OH_relation(CO, OH, Mstellar, kappa, counter, output_name_list, output_file)
    output_file = f"{output_path}/Galactic_NO_OH_comparison.png"
    plot_NO_OH_relation(NO, OH, Mstellar, kappa, counter, output_name_list, output_file)
