import matplotlib.pylab as plt
from matplotlib.pylab import rcParams
import numpy as np
from .loadObservationalData import plot_GALAH_data, plot_MW_data, plot_StrontiumObsData, plot_APOGEE_data
from .plot_mass_metallicity import plot_Kirby_data, plot_Kirby_analysed, \
    plot_gallazzi, plot_gallazzi_2005, plot_Kudritzki_2016, plot_Zahid_2017


def compare_stellar_abundances(sims_data, output_name_list, output_path):

    O_Fe_all = sims_data['O_Fe']
    Fe_H_all = sims_data['Fe_H']
    Mg_Fe_all = sims_data['Mg_Fe']
    FeSNIa_Fe = sims_data['FeSNIa_Fe']
    counter = sims_data['counter']


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
        "lines.linewidth": 0.2,
    }
    rcParams.update(params)
    fig = plt.figure()

    # Box stellar abundance --------------------------------
    ax = plt.subplot(1, 1, 1)
    plt.grid("True")

    count = 0
    color = ['tab:blue', 'tab:green', 'tab:orange', 'crimson', 'tab:purple']
    for i in range(len(output_name_list)):
        xm = Fe_H_all[count:count + counter[i]]
        ym = FeSNIa_Fe[count:count + counter[i]]
        count += counter[i]
        plt.plot(xm, ym, '-', lw=1.5, color=color[i], label=output_name_list[i])

    plt.xlabel("[Fe/H]", labelpad=2)
    plt.ylabel("Fe(SNIa)/Fe", labelpad=2)
    plt.yscale('log')
    plt.axis([-4, 1, 1e-2,2])
    ax.tick_params(direction='in', axis='both', which='both')
    plt.legend(loc=[0.0, 0.02], labelspacing=0.1, handlelength=1.5, handletextpad=0.1, frameon=False, ncol=1,
               columnspacing=0.02)
    plt.savefig(f"{output_path}/FeSNIa_Fe_comparison.png", dpi=200)


    ##################
    fig = plt.figure()

    # Box stellar abundance --------------------------------
    ax = plt.subplot(1, 1, 1)
    plt.grid("True")

    plot_MW_data('O')
    plot_GALAH_data('O')

    count = 0
    color = ['tab:blue','tab:green','tab:orange','crimson','tab:purple']
    for i in range(len(output_name_list)):
        xm = Fe_H_all[count:count+counter[i]]
        ym = O_Fe_all[count:count+counter[i]]
        count += counter[i]

        if i==0 :plt.plot(xm, ym, '-', lw=0.5, color='tab:blue', label='GALAH DR3')
        plt.plot(xm, ym, '-', lw=1.5, color=color[i], label=output_name_list[i])

    plt.text(-3.8, 1.3, "MW-type galaxies")
    plt.xlabel("[Fe/H]", labelpad=2)
    plt.ylabel("[O/Fe]", labelpad=2)
    plt.axis([-4, 1, -1, 1.5])
    plt.legend(loc=[0.0, 0.02], labelspacing=0.1, handlelength=1.5, handletextpad=0.1, frameon=False, ncol=1,
               columnspacing=0.02)
    plt.savefig(f"{output_path}/O_Fe_comparison.png", dpi=200)

    ##################
    fig = plt.figure()

    # Box stellar abundance --------------------------------
    ax = plt.subplot(1, 1, 1)
    plt.grid("True")

    plot_APOGEE_data('O')

    count = 0
    color = ['tab:blue','tab:green','tab:orange','crimson','tab:purple']
    for i in range(len(output_name_list)):
        xm = Fe_H_all[count:count+counter[i]]
        ym = O_Fe_all[count:count+counter[i]]
        count += counter[i]

        if i==0 :plt.plot(xm, ym, '-', lw=0.5, color='blue', label='APOGEE data')
        plt.plot(xm, ym, '-', lw=1.5, color=color[i], label=output_name_list[i], zorder=200)

    plt.text(-3.8, 1.3, "MW-type galaxies")
    plt.xlabel("[Fe/H]", labelpad=2)
    plt.ylabel("[O/Fe]", labelpad=2)
    plt.axis([-4, 1, -1, 1.5])
    plt.legend(loc=[0.0, 0.02], labelspacing=0.1, handlelength=1.5, handletextpad=0.1, frameon=False, ncol=1,
               columnspacing=0.02)
    plt.savefig(f"{output_path}/O_Fe_comparison_2.png", dpi=200)

    ########################
    fig = plt.figure()
    ax = plt.subplot(1, 1, 1)
    plt.grid("True")

    plot_MW_data('Mg')
    plot_GALAH_data('Mg')

    count = 0
    for i in range(len(output_name_list)):
        xm = Fe_H_all[count:count + counter[i]]
        ym = Mg_Fe_all[count:count + counter[i]]
        count += counter[i]

        if i == 0: plt.plot(xm, ym, '-', lw=0.5, color='tab:blue', label='GALAH DR3')
        plt.plot(xm, ym, '-', lw=1.5, color=color[i], label=output_name_list[i])

    plt.xlabel("[Fe/H]", labelpad=2)
    plt.ylabel("[Mg/Fe]", labelpad=2)
    plt.text(-3.8, 1.2, "MW-type galaxies")
    plt.axis([-4, 1, -2, 1.5])
    plt.legend(loc=[0, 0.02], labelspacing=0.1, handlelength=1.5, handletextpad=0.1, frameon=False, ncol=1,
               columnspacing=0.02)
    plt.savefig(f"{output_path}/Mg_Fe_comparison.png", dpi=200)

    ########################
    fig = plt.figure()
    ax = plt.subplot(1, 1, 1)
    plt.grid("True")

    plot_APOGEE_data('MG')

    count = 0
    for i in range(len(output_name_list)):
        xm = Fe_H_all[count:count + counter[i]]
        ym = Mg_Fe_all[count:count + counter[i]]
        count += counter[i]

        if i == 0: plt.plot(xm, ym, '-', lw=0.5, color='blue', label='APOGEE data')
        plt.plot(xm, ym, '-', lw=1.5, color=color[i], label=output_name_list[i], zorder=200)

    plt.xlabel("[Fe/H]", labelpad=2)
    plt.ylabel("[Mg/Fe]", labelpad=2)
    plt.text(-3.8, 1.2, "MW-type galaxies")
    plt.axis([-4, 1, -2, 1.5])
    plt.legend(loc=[0, 0.02], labelspacing=0.1, handlelength=1.5, handletextpad=0.1, frameon=False, ncol=1,
               columnspacing=0.02)
    plt.savefig(f"{output_path}/Mg_Fe_comparison_2.png", dpi=200)

    ########################
    # Load data:
    C_Fe_all = sims_data['C_Fe']
    Ba_Fe_all = sims_data['Ba_Fe']
    Sr_Fe_all = sims_data['Sr_Fe']
    Eu_Fe_all = sims_data['Eu_Fe']
    Si_Fe_all = sims_data['Si_Fe']
    N_Fe_all = sims_data['N_Fe']
    Ne_Fe_all = sims_data['Ne_Fe']

    # make remaining plots with just GALAH data (Buder+21)
    for el in ['C', 'Si', 'Eu', 'Ba']:
        fig = plt.figure(figsize=(3.8, 3))
        ax = plt.subplot(1, 1, 1)
        plt.grid("True")

        count = 0
        for i in range(len(output_name_list)):
            xm = Fe_H_all[count:count + counter[i]]
            if el == 'C': ym = C_Fe_all[count:count + counter[i]]
            if el == 'Ba': ym = Ba_Fe_all[count:count + counter[i]]
            if el == 'Eu': ym = Eu_Fe_all[count:count + counter[i]]
            if el == 'Si': ym = Si_Fe_all[count:count + counter[i]]
            count += counter[i]

            if i == 0: plt.plot(xm, ym, '-', lw=0.5, color='tab:blue', label='GALAH DR3')
            plt.plot(xm, ym, '-', lw=1.5, color=color[i], label=output_name_list[i])

        plot_GALAH_data(el)
        plt.xlabel("[Fe/H]", labelpad=2)
        plt.ylabel(f"[{el}/Fe]", labelpad=2)
        plt.text(-3.8, 1.2, "MW-type galaxies")
        plt.axis([-4, 1, -2, 1.5])
        plt.legend(loc=[0, 0.02], labelspacing=0.1, handlelength=1.5, handletextpad=0.1, frameon=False, ncol=1,
                   columnspacing=0.02)
        plt.tight_layout()
        plt.savefig(f"{output_path}/{el}_Fe_comparison.png", dpi=200)

    #######
    fig = plt.figure(figsize=(3.8, 3))
    ax = plt.subplot(1, 1, 1)
    plt.grid("True")

    count = 0
    for i in range(len(output_name_list)):
        xm = Fe_H_all[count:count + counter[i]]
        ym = C_Fe_all[count:count + counter[i]]
        count += counter[i]

        if i == 0: plt.plot(xm, ym, '-', lw=0.5, color='blue', label='APOGEE data')
        plt.plot(xm, ym, '-', lw=1.5, color=color[i], label=output_name_list[i])

    plot_APOGEE_data('C')
    plt.xlabel("[Fe/H]", labelpad=2)
    plt.ylabel("[C/Fe]", labelpad=2)
    plt.text(-3.8, 1.2, "MW-type galaxies")
    plt.axis([-4, 1, -2, 1.5])
    plt.legend(loc=[0, 0.02], labelspacing=0.1, handlelength=1.5, handletextpad=0.1, frameon=False, ncol=1,
               columnspacing=0.02)
    plt.tight_layout()
    plt.savefig(f"{output_path}/C_Fe_comparison_2.png", dpi=200)

    for el in ['N', 'Ne']:
        fig = plt.figure(figsize=(3.8, 3))
        ax = plt.subplot(1, 1, 1)
        plt.grid("True")

        count = 0
        for i in range(len(output_name_list)):
            xm = Fe_H_all[count:count + counter[i]]
            if el == 'N': ym = N_Fe_all[count:count + counter[i]]
            if el == 'Ne': ym = Ne_Fe_all[count:count + counter[i]]
            count += counter[i]
            if (i == 0) & (el == 'N'): plt.plot(xm, ym, '-', lw=0.5, color='blue', label='APOGEE data')
            plt.plot(xm, ym, '-', lw=1.5, color=color[i], label=output_name_list[i])

        if el == 'N': plot_APOGEE_data('N')
        plt.xlabel("[Fe/H]", labelpad=2)
        plt.ylabel(f"[{el}/Fe]", labelpad=2)
        plt.text(-3.8, 1.2, "MW-type galaxies")
        plt.axis([-4, 1, -2, 1.5])
        plt.legend(loc=[0, 0.02], labelspacing=0.1, handlelength=1.5, handletextpad=0.1, frameon=False, ncol=1,
                   columnspacing=0.02)
        plt.tight_layout()
        plt.savefig(f"{output_path}/{el}_Fe_comparison.png", dpi=200)

    fig = plt.figure()
    ax = plt.subplot(1, 1, 1)
    plt.grid("True")

    count = 0
    for i in range(len(output_name_list)):
        xm = Fe_H_all[count:count + counter[i]]
        ym = Sr_Fe_all[count:count + counter[i]]
        count += counter[i]
        plt.plot(xm, ym, '-', lw=1.5, color=color[i], label=output_name_list[i])

    plot_StrontiumObsData()

    plt.xlabel("[Fe/H]", labelpad=2)
    plt.ylabel("[Sr/Fe]", labelpad=2)
    plt.text(-3.8, 1.2, "MW-type galaxies")
    plt.axis([-4, 1, -2, 1.5])
    plt.legend(loc=[0, 0.02], labelspacing=0.1, handlelength=1.5, handletextpad=0.1, frameon=False, ncol=1,
               columnspacing=0.02)
    plt.savefig(f"{output_path}/Sr_Fe_comparison.png", dpi=200)

def plot_Fe_H_mass_relation(Mstellar, Fe_H, counter, ylabel, output_name_list, output_file):

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

    plot_Kirby_data()
    plot_Kirby_analysed()
    plot_gallazzi_2005()
    plot_Kudritzki_2016()
    plot_Zahid_2017()

    count = 0
    color = ['tab:blue','tab:green','tab:orange','crimson','tab:purple']
    for i in range(len(output_name_list)):
        xm = Mstellar[count:count + counter[i]]
        ym = Fe_H[count:count + counter[i]]
        count += counter[i]

        bins = np.arange(6, 12, 0.2)
        ind = np.digitize(xm, bins)
        ylo = [np.percentile(10**ym[ind == i],16) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
        yhi = [np.percentile(10**ym[ind == i],84) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
        ym = [np.median(10**ym[ind == i]) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
        xm = [np.median(10**xm[ind == i]) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
        plt.fill_between(xm, ylo, yhi, color=color[i], alpha=0.2)
        plt.plot(xm, ym, '-', lw=1.5, color=color[i], label=output_name_list[i])

    plt.ylabel(ylabel, labelpad=2)
    plt.xlabel("Stellar Mass [M$_{\odot}$]", labelpad=2)
    plt.xscale('log')
    plt.yscale('log')
    ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    plt.axis([1e5, 1e12, 1e-3, 1e2])
    handles, labels = plt.gca().get_legend_handles_labels()
    order = np.arange(len(handles)-2)+2
    order = np.append(order,0)
    order = np.append(order,1)

    plt.legend([handles[idx] for idx in order], [labels[idx] for idx in order],
               loc='upper left', labelspacing=0.1, handlelength=1.5, handletextpad=0.1, frameon=False, ncol=1,
               fontsize=9, columnspacing=0.02)

    plt.savefig(output_file, dpi=200)


def plot_Mg_Fe_mass_relation(Mstellar, Mg_Fe, counter, ylabel, output_name_list, output_file):

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

    plot_gallazzi('Mg')

    count = 0
    color = ['tab:blue','tab:green','tab:orange','crimson','tab:purple']
    for i in range(len(output_name_list)):
        xm = Mstellar[count:count + counter[i]]
        ym = Mg_Fe[count:count + counter[i]]
        count += counter[i]

        bins = np.arange(6, 12, 0.2)
        ind = np.digitize(xm, bins)
        ylo = [np.percentile(ym[ind == i],16) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
        yhi = [np.percentile(ym[ind == i],84) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
        ym = [np.median(ym[ind == i]) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
        xm = [np.median(10**xm[ind == i]) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
        plt.fill_between(xm, ylo, yhi, color=color[i], alpha=0.2)
        plt.plot(xm, ym, '-', lw=1.5, color=color[i], label=output_name_list[i])

    plt.ylabel(ylabel, labelpad=2)
    plt.xlabel("Stellar Mass [M$_{\odot}$]", labelpad=2)
    plt.xscale('log')
    plt.axis([1e8, 1e12, -0.2, 0.6])
    ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    plt.legend(loc='upper left', labelspacing=0.1, handlelength=1.5, handletextpad=0.1, frameon=False, ncol=1,
               fontsize=9, columnspacing=0.02)

    plt.savefig(output_file, dpi=200)

def plot_O_Fe_mass_relation(Mstellar, O_Fe, counter, ylabel, output_name_list, output_file):

    plt.figure()

    # Box stellar abundance --------------------------------
    ax = plt.subplot(1, 1, 1)
    plt.grid("True")

    plot_gallazzi('O')

    count = 0
    color = ['tab:blue','tab:green','tab:orange','crimson','tab:purple']
    for i in range(len(output_name_list)):
        xm = Mstellar[count:count + counter[i]]
        ym = O_Fe[count:count + counter[i]]
        count += counter[i]

        bins = np.arange(6, 12, 0.2)
        ind = np.digitize(xm, bins)
        ylo = [np.percentile(ym[ind == i],16) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
        yhi = [np.percentile(ym[ind == i],84) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
        ym = [np.median(ym[ind == i]) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
        xm = [np.median(10**xm[ind == i]) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
        plt.fill_between(xm, ylo, yhi, color=color[i], alpha=0.2)
        plt.plot(xm, ym, '-', lw=1.5, color=color[i], label=output_name_list[i])

    plt.ylabel(ylabel, labelpad=2)
    plt.xlabel("Stellar Mass [M$_{\odot}$]", labelpad=2)
    plt.xscale('log')
    plt.axis([1e8, 1e12, 0.0, 0.6])
    ax.tick_params(direction='in', axis='both', which='both', pad=4.5)

    plt.legend(loc='upper left', labelspacing=0.1, handlelength=1.5, handletextpad=0.1, frameon=False, ncol=1,
               fontsize=9, columnspacing=0.02)

    plt.savefig(output_file, dpi=200)

def plot_metallicity_mass_relation(Mstellar, Z, counter, ylabel, output_name_list, output_file):

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

    # Box stellar abundance --------------------------------
    ax = plt.subplot(1, 1, 1)
    plt.grid("True")

    plot_Kirby_data()
    plot_Kirby_analysed()
    plot_gallazzi_2005()
    plot_Kudritzki_2016()
    plot_Zahid_2017()

    count = 0
    color = ['tab:blue','tab:green','tab:orange','crimson','tab:purple']
    for i in range(len(output_name_list)):
        xm = Mstellar[count:count + counter[i]]
        ym = Z[count:count + counter[i]]
        count += counter[i]

        bins = np.arange(6, 12, 0.2)
        ind = np.digitize(xm, bins)
        ylo = [np.percentile(ym[ind == i],16) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
        yhi = [np.percentile(ym[ind == i],84) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
        ym = [np.median(ym[ind == i]) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
        xm = [np.median(10**xm[ind == i]) for i in range(1, len(bins)) if len(xm[ind == i]) > 2]
        plt.fill_between(xm, ylo, yhi, color=color[i], alpha=0.2)
        plt.plot(xm, ym, '-', lw=1.5, color=color[i], label=output_name_list[i])

    plt.ylabel(ylabel, labelpad=2)
    plt.xlabel("Stellar Mass [M$_{\odot}$]", labelpad=2)
    plt.xscale('log')
    plt.yscale('log')
    ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    plt.axis([1e5, 1e12, 1e-3, 1e2])

    handles, labels = plt.gca().get_legend_handles_labels()
    order = np.arange(len(handles)-2)+2
    order = np.append(order,0)
    order = np.append(order,1)

    plt.legend([handles[idx] for idx in order], [labels[idx] for idx in order],
               loc='upper left', labelspacing=0.1, handlelength=1.5, handletextpad=0.1, frameon=False, ncol=1,
               fontsize=9, columnspacing=0.02)

    plt.savefig(output_file, dpi=200)

def compare_mass_metallicity_relations(sim_data, output_name_list, output_path):

    Mstellar = sim_data['Mstellar']
    counter = sim_data['counter']

    ##############  Fe/H  #################

    FeH_mw = sim_data['FeH_log_mass_weighted']
    FeH_mw_r =sim_data['FeH_log_mass_weighted_ratio']
    FeH_mw_2 = sim_data['FeH_mass_weighted_log']

    ylabel = r"Stellar $10^{[\langle\mathrm{Fe}\rangle_{m}/\langle\mathrm{H}\rangle_{m}]}$"
    output_file = f"{output_path}/Mstellar_FeH_mw_comparison.png"
    plot_Fe_H_mass_relation(Mstellar, FeH_mw, counter, ylabel, output_name_list, output_file)

    ylabel = r"Stellar $10^{[\langle\mathrm{Fe/H}\rangle_{m}]}$"
    output_file = f"{output_path}/Mstellar_FeH_mwr_comparison.png"
    plot_Fe_H_mass_relation(Mstellar, FeH_mw_r, counter, ylabel, output_name_list, output_file)

    ylabel = r"Stellar $10^{\langle[\mathrm{Fe/H}]\rangle_{m}}$"
    output_file = f"{output_path}/Mstellar_FeH_mwr2_comparison.png"
    plot_Fe_H_mass_relation(Mstellar, FeH_mw_2, counter, ylabel, output_name_list, output_file)

    FeH_mw = sim_data['FeH_log_light_weighted']
    FeH_mw_r =sim_data['FeH_log_light_weighted_ratio']
    FeH_mw_2 = sim_data['FeH_light_weighted_log']

    ylabel = r"Stellar $10^{[\langle\mathrm{Fe}\rangle_{l}/\langle\mathrm{H}\rangle_{l}]}$"
    output_file = f"{output_path}/Mstellar_FeH_lw_comparison.png"
    plot_Fe_H_mass_relation(Mstellar, FeH_mw, counter, ylabel, output_name_list, output_file)

    ylabel = r"Stellar $10^{[\langle\mathrm{Fe/H}\rangle_{l}]}$"
    output_file = f"{output_path}/Mstellar_FeH_lwr_comparison.png"
    plot_Fe_H_mass_relation(Mstellar, FeH_mw_r, counter, ylabel, output_name_list, output_file)

    ylabel = r"Stellar $10^{\langle[\mathrm{Fe/H}]\rangle_{l}}$"
    output_file = f"{output_path}/Mstellar_FeH_lwr2_comparison.png"
    plot_Fe_H_mass_relation(Mstellar, FeH_mw_2, counter, ylabel, output_name_list, output_file)

    FeH_mw = sim_data['FeH_log_metallicity_weighted']
    FeH_mw_r = sim_data['FeH_log_metallicity_weighted_ratio']
    FeH_mw_2 = sim_data['FeH_metallicity_weighted_log']

    ylabel = r"Stellar $10^{[\langle\mathrm{Fe}\rangle_{Z}/\langle\mathrm{H}\rangle_{Z}]}$"
    output_file = f"{output_path}/Mstellar_FeH_zw_comparison.png"
    plot_Fe_H_mass_relation(Mstellar, FeH_mw, counter, ylabel, output_name_list, output_file)

    ylabel = r"Stellar $10^{[\langle\mathrm{Fe/H}\rangle_{Z}]}$"
    output_file = f"{output_path}/Mstellar_FeH_zwr_comparison.png"
    plot_Fe_H_mass_relation(Mstellar, FeH_mw_r, counter, ylabel, output_name_list, output_file)

    ylabel = r"Stellar $10^{\langle[\mathrm{Fe/H}]\rangle_{Z}}$"
    output_file = f"{output_path}/Mstellar_FeH_zwr2_comparison.png"
    plot_Fe_H_mass_relation(Mstellar, FeH_mw_2, counter, ylabel, output_name_list, output_file)

    ###################
    ###################

    FeH_mw = sim_data['MgFe_log_mass_weighted']
    FeH_mw_r = sim_data['MgFe_log_mass_weighted_ratio']
    FeH_mw_2 = sim_data['MgFe_mass_weighted_log']

    ylabel = r"$[\langle\mathrm{Mg}\rangle_{m}/\langle\mathrm{Fe}\rangle_{m}]$"
    output_file = f"{output_path}/Mstellar_MgFe_mw_comparison.png"
    plot_Mg_Fe_mass_relation(Mstellar, FeH_mw, counter, ylabel, output_name_list, output_file)

    ylabel = r"$[\langle\mathrm{Mg/Fe}\rangle_{m}]$"
    output_file = f"{output_path}/Mstellar_MgFe_mwr_comparison.png"
    plot_Mg_Fe_mass_relation(Mstellar, FeH_mw_r, counter, ylabel, output_name_list, output_file)

    ylabel = r"$\langle[\mathrm{Mg/Fe}]\rangle_{m}$"
    output_file = f"{output_path}/Mstellar_MgFe_mwr2_comparison.png"
    plot_Mg_Fe_mass_relation(Mstellar, FeH_mw_2, counter, ylabel, output_name_list, output_file)

    FeH_mw = sim_data['MgFe_log_light_weighted']
    FeH_mw_r = sim_data['MgFe_log_light_weighted_ratio']
    FeH_mw_2 = sim_data['MgFe_light_weighted_log']

    ylabel = r"$[\langle\mathrm{Mg}\rangle_{l}/\langle\mathrm{Fe}\rangle_{l}]$"
    output_file = f"{output_path}/Mstellar_MgFe_lw_comparison.png"
    plot_Mg_Fe_mass_relation(Mstellar, FeH_mw, counter, ylabel, output_name_list, output_file)

    ylabel = r"$[\langle\mathrm{Mg/Fe}\rangle_{l}]$"
    output_file = f"{output_path}/Mstellar_MgFe_lwr_comparison.png"
    plot_Mg_Fe_mass_relation(Mstellar, FeH_mw_r, counter, ylabel, output_name_list, output_file)

    ylabel = r"$\langle[\mathrm{Mg/Fe}]\rangle_{l}$"
    output_file = f"{output_path}/Mstellar_MgFe_lwr2_comparison.png"
    plot_Mg_Fe_mass_relation(Mstellar, FeH_mw_2, counter, ylabel, output_name_list, output_file)

    FeH_mw = sim_data['MgFe_log_metallicity_weighted']
    FeH_mw_r = sim_data['MgFe_log_metallicity_weighted_ratio']
    FeH_mw_2 = sim_data['MgFe_metallicity_weighted_log']

    ylabel = r"$[\langle\mathrm{Mg}\rangle_{Z}/\langle\mathrm{Fe}\rangle_{Z}]$"
    output_file = f"{output_path}/Mstellar_MgFe_zw_comparison.png"
    plot_Mg_Fe_mass_relation(Mstellar, FeH_mw, counter, ylabel, output_name_list, output_file)

    ylabel = r"$[\langle\mathrm{Mg/Fe}\rangle_{Z}]$"
    output_file = f"{output_path}/Mstellar_MgFe_zwr_comparison.png"
    plot_Mg_Fe_mass_relation(Mstellar, FeH_mw_r, counter, ylabel, output_name_list, output_file)

    ylabel = r"$\langle[\mathrm{Mg/Fe}]\rangle_{Z}$"
    output_file = f"{output_path}/Mstellar_MgFe_zwr2_comparison.png"
    plot_Mg_Fe_mass_relation(Mstellar, FeH_mw_2, counter, ylabel, output_name_list, output_file)

    ################
    ################

    FeH_mw = sim_data['OFe_log_mass_weighted']
    FeH_mw_r = sim_data['OFe_log_mass_weighted_ratio']
    FeH_mw_2 = sim_data['OFe_mass_weighted_log']

    ylabel = r"$[\langle\mathrm{O}\rangle_{m}/\langle\mathrm{Fe}\rangle_{m}]$"
    output_file = f"{output_path}/Mstellar_OFe_mw_comparison.png"
    plot_O_Fe_mass_relation(Mstellar, FeH_mw, counter, ylabel, output_name_list, output_file)

    ylabel = r"$[\langle\mathrm{O/Fe}\rangle_{m}]$"
    output_file = f"{output_path}/Mstellar_OFe_mwr_comparison.png"
    plot_O_Fe_mass_relation(Mstellar, FeH_mw_r, counter, ylabel, output_name_list, output_file)

    ylabel = r"$\langle[\mathrm{O/Fe}]\rangle_{m}$"
    output_file = f"{output_path}/Mstellar_OFe_mwr2_comparison.png"
    plot_O_Fe_mass_relation(Mstellar, FeH_mw_2, counter, ylabel, output_name_list, output_file)

    FeH_mw = sim_data['OFe_log_light_weighted']
    FeH_mw_r = sim_data['OFe_log_light_weighted_ratio']
    FeH_mw_2 = sim_data['OFe_light_weighted_log']

    ylabel = r"$[\langle\mathrm{O}\rangle_{l}/\langle\mathrm{Fe}\rangle_{l}]$"
    output_file = f"{output_path}/Mstellar_OFe_lw_comparison.png"
    plot_O_Fe_mass_relation(Mstellar, FeH_mw, counter, ylabel, output_name_list, output_file)

    ylabel = r"$[\langle\mathrm{O/Fe}\rangle_{l}]$"
    output_file = f"{output_path}/Mstellar_OFe_lwr_comparison.png"
    plot_O_Fe_mass_relation(Mstellar, FeH_mw_r, counter, ylabel, output_name_list, output_file)

    ylabel = r"$\langle[\mathrm{O/Fe}]\rangle_{l}$"
    output_file = f"{output_path}/Mstellar_OFe_lwr2_comparison.png"
    plot_O_Fe_mass_relation(Mstellar, FeH_mw_2, counter, ylabel, output_name_list, output_file)

    FeH_mw = sim_data['OFe_log_metallicity_weighted']
    FeH_mw_r = sim_data['OFe_log_metallicity_weighted_ratio']
    FeH_mw_2 = sim_data['OFe_metallicity_weighted_log']

    ylabel = r"$[\langle\mathrm{O}\rangle_{Z}/\langle\mathrm{Fe}\rangle_{Z}]$"
    output_file = f"{output_path}/Mstellar_OFe_zw_comparison.png"
    plot_O_Fe_mass_relation(Mstellar, FeH_mw, counter, ylabel, output_name_list, output_file)

    ylabel = r"$[\langle\mathrm{O/Fe}\rangle_{Z}]$"
    output_file = f"{output_path}/Mstellar_OFe_zwr_comparison.png"
    plot_O_Fe_mass_relation(Mstellar, FeH_mw_r, counter, ylabel, output_name_list, output_file)

    ylabel = r"$\langle[\mathrm{O/Fe}]\rangle_{Z}$"
    output_file = f"{output_path}/Mstellar_OFe_zwr2_comparison.png"
    plot_O_Fe_mass_relation(Mstellar, FeH_mw_2, counter, ylabel, output_name_list, output_file)

    ###################
    ###################

    Z_mass_weighted = sim_data['Z_mass_weighted']
    ylabel = "Stellar (mass-weighted) $Z/Z_{\odot}$"
    output_file = f"{output_path}/Mstellar_Z_mw_comparison.png"
    plot_metallicity_mass_relation(Mstellar, Z_mass_weighted, counter, ylabel, output_name_list, output_file)

    Z_light_weighted = sim_data['Z_light_weighted']
    ylabel = "Stellar (light-weighted r-band) $Z/Z_{\odot}$"
    output_file = f"{output_path}/Mstellar_Z_lw_comparison.png"
    plot_metallicity_mass_relation(Mstellar, Z_light_weighted, counter, ylabel, output_name_list, output_file)


