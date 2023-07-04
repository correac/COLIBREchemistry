import numpy as np
import matplotlib.pylab as plt
from matplotlib.pylab import rcParams

def test(sim_info):

    sample = np.where(sim_info.halo_data.log10_stellar_mass >= 10.7)[0]

    centrals = np.where(sim_info.halo_data.type[sample] == 10)[0]

    sample = sample[centrals]
    num_sample = len(sample)
    print(sim_info.halo_data.log10_stellar_mass[sample])

    mp_in_cgs = 1.6737236e-24
    mH_in_cgs = 1.00784 * mp_in_cgs
    mFe_in_cgs = 55.845 * mp_in_cgs

    # Asplund et al. (2009)
    Fe_H_Sun = 7.5
    Fe_H_Sun = Fe_H_Sun - 12.0 - np.log10(mH_in_cgs / mFe_in_cgs)
    print(Fe_H_Sun, 10**Fe_H_Sun)

    xbins = np.arange(-7, 1, 0.1)

    # Plot parameters
    params = {
        "font.size": 11,
        "font.family": "STIXGeneral",
        "text.usetex": False,
        "mathtext.fontset": "stix",
        "figure.figsize": (4, 3),
        "figure.subplot.left": 0.18,
        "figure.subplot.right": 0.95,
        "figure.subplot.bottom": 0.15,
        "figure.subplot.top": 0.95,
        "lines.markersize": 0.5,
        "lines.linewidth": 1,
    }
    rcParams.update(params)
    fig = plt.figure()
    ax = plt.subplot(1, 1, 1)
    #plt.grid("True")

    for i in range(1):

        sim_info.make_particle_data(halo_id=sim_info.halo_data.halo_ids[sample[i]])
        H_stars = sim_info.stars.hydrogen.copy()
        Fe_stars = sim_info.stars.iron.copy()
        M_stars = sim_info.stars.mass.copy()

        Fe_H = Fe_stars / H_stars
        Fe_H = np.log10(Fe_H) - Fe_H_Sun

        hist, bin_edges = np.histogram(Fe_H, bins=xbins)
        bin_edges = (bin_edges[1:]+bin_edges[:-1]) * 0.5
        plt.plot(10**bin_edges, hist)

        Fe_H_mass_weighted = np.sum( M_stars * Fe_stars) / np.sum( H_stars * M_stars)
        Fe_H_mass_weighted = np.log10(Fe_H_mass_weighted) - Fe_H_Sun

        Fe_H_ratio_mass_weighted = np.sum(M_stars * (Fe_stars/H_stars)) / np.sum(M_stars)
        Fe_H_ratio_mass_weighted = np.log10(Fe_H_ratio_mass_weighted) - Fe_H_Sun

        # select = np.where(Fe_stars < 1e-4)[0]
        # Fe = Fe_stars.copy()
        # Fe[select] = 1e-4
        # H = H_stars.copy()
        # M = M_stars.copy()
        # Fe_H_log_mass_weighted = np.log10(Fe/H) - Fe_H_Sun
        # Fe_H_log_mass_weighted = np.sum(M * Fe_H_log_mass_weighted) / np.sum(M)

        Fe_H_Sun_2 = 2.82e-5
        Fe = Fe_stars.copy()
        H = H_stars.copy()
        M = M_stars.copy()
        Fe_abundance = Fe / (55.845 * H)
        Fe_H_log_mass_weighted = np.log10( np.clip(Fe_abundance, Fe_H_Sun_2 * 1e-4, np.inf) )
        Fe_H_log_mass_weighted = np.sum(M * Fe_H_log_mass_weighted) / np.sum(M)
        Fe_H_log_mass_weighted -= np.log10(Fe_H_Sun_2)

        # select = np.where(Fe_stars < 1e-3)[0]
        # Fe = Fe_stars.copy()
        # Fe[select] = 1e-3
        # H = H_stars.copy()
        # M = M_stars.copy()
        # Fe_H_log_mass_weighted_lo = np.log10(Fe/H) - Fe_H_Sun
        # Fe_H_log_mass_weighted_lo = np.sum(M * Fe_H_log_mass_weighted_lo) / np.sum(M)
        Fe_H_log_mass_weighted_lo = np.log10( np.clip(Fe_abundance, Fe_H_Sun_2 * 1e-1, np.inf) )
        Fe_H_log_mass_weighted_lo = np.sum(M * Fe_H_log_mass_weighted_lo) / np.sum(M)
        Fe_H_log_mass_weighted_lo -= np.log10(Fe_H_Sun_2)

        print(Fe_H_log_mass_weighted, Fe_H_log_mass_weighted_lo, Fe_H_mass_weighted)

        xrange = np.array([1,1])
        yrange = np.array([1,5e3])
        plt.plot(xrange * 10**Fe_H_mass_weighted, yrange, label='$\log_{10}$(<Fe>/<H>)')
        plt.plot(xrange * 10**Fe_H_ratio_mass_weighted, yrange, label='$\log_{10}$(<Fe/H>)')
        plt.plot(xrange * 10**Fe_H_log_mass_weighted, yrange, label='<[Fe/H]>, [F/H]=-4 floor')
        plt.plot(xrange * 10**Fe_H_log_mass_weighted_lo, yrange, label='<[Fe/H]>, [F/H]=-1 floor')


    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(r'$10^{[\mathrm{Fe/H}]}$', labelpad=2)
    plt.ylabel('Histogram', labelpad=2)
    ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    plt.axis([1e-6, 10, 1, 5e3])
    plt.legend(loc='upper left', labelspacing=0.1, handlelength=1., handletextpad=0.1, frameon=False, ncol=1,
               fontsize=9, columnspacing=0.02)

    plt.savefig('test.png', dpi=300)