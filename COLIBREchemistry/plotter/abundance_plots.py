import matplotlib
import matplotlib.pylab as plt
from matplotlib.pylab import rcParams
import numpy as np
from tqdm import tqdm
from .loadObservationalData import plot_MW_Satellites, plot_MW_data, \
    plot_GALAH_data, plot_StrontiumObsData, plot_APOGEE_data
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from .stellar_abundances import calculate_abundaces_from_MW_type_galaxies
import h5py


def plot_contour(x, y, xmin, xmax, ymin, ymax, cmap="autumn"):

    ngridx = 30
    ngridy = 30

    # Create grid values first.
    xi = np.linspace(xmin, xmax, ngridx)
    yi = np.linspace(ymin, ymax, ngridy)

    # Create a histogram
    h, xedges, yedges = np.histogram2d(x, y, bins=(xi, yi))
    xbins = xedges[:-1] + (xedges[1] - xedges[0]) / 2
    ybins = yedges[:-1] + (yedges[1] - yedges[0]) / 2

    z = h.T

    levels = np.arange(10, np.ceil(h.max()), 10)

    if len(levels) > 2:
        contour = plt.contourf(xbins, ybins, z, levels=levels, linewidths=0.8, cmap=cmap, alpha=0.5)
    else:
        contour = plt.contour(xbins, ybins, z, levels=levels, linewidths=0.8, cmap=cmap)
    return contour


def make_plot(x, y, GalR, Galz, Galz_min, Galz_max, GalR_min, GalR_max):
    plt.grid(linestyle='-', linewidth=0.3)

    select = np.where((GalR >= GalR_min) & (GalR <= GalR_max)
                      & (np.abs(Galz) <= Galz_max) & (np.abs(Galz) >= Galz_min))[0]
    # plt.plot(x, y, 'o', alpha=0.05, color='grey')
    contour = plot_contour(x[select], y[select], -1, 1, -0.6, 0.5)
    legend_elements, _ = contour.legend_elements()

    plt.axis([-1, 1, -0.6, 0.5])
    plt.tick_params(direction='in', axis='both', which='both', pad=4.5)

    xticks = np.array([-1, -0.5, 0, 0.5])
    labels = ["-1", "-0.5", "0", "0.5"]
    plt.xticks(xticks, labels)

    yticks = np.array([-0.5, -0.25, 0, 0.25, 0.5])
    labels = ["-0.5","-0.25", "0", "0.25", "0.5"]
    plt.yticks(yticks, labels)

    bins = np.arange(-1, 1, 0.2)
    ind = np.digitize(x[select], bins)
    xm = [np.median(x[select[ind == i]]) for i in range(1, len(bins)) if len(x[select[ind == i]]) > 10]
    ym = [np.median(y[select[ind == i]]) for i in range(1, len(bins)) if len(x[select[ind == i]]) > 10]
    plt.plot(xm, ym, '-', lw=1.5, color='black')

    return legend_elements


def add_APOGEE(element, Galz_min, Galz_max, GalR_min, GalR_max):

    file = "./plotter/obs_data/APOGEE_data_V2.hdf5"
    with h5py.File(file, 'r') as APOGEE_dataobj:
        FE_H = APOGEE_dataobj['FE_H'][:]
        el_FE = APOGEE_dataobj[f"{element}_FE"][:]
        Galz = APOGEE_dataobj['Galz'][:]
        GalR = APOGEE_dataobj['GalR'][:]

    Fe_over_H = 7.5
    Fe_over_H_GA07 = 7.45
    if element == "MG":
        Mg_over_H = 7.6
        Mg_over_Fe_AS09 = Mg_over_H - Fe_over_H
        Mg_over_H_GA07 = 7.53
        Mg_over_Fe_GA07 = Mg_over_H_GA07 - Fe_over_H_GA07
        correction = Mg_over_Fe_GA07 - Mg_over_Fe_AS09
    elif element == "O":
        O_over_H = 8.69
        O_over_Fe_AS09 = O_over_H - Fe_over_H
        O_over_H_GA07 = 8.66
        O_over_Fe_GA07 = O_over_H_GA07 - Fe_over_H_GA07
        correction = O_over_Fe_GA07 - O_over_Fe_AS09
    elif element == "N":
        N_over_H = 7.83
        N_over_Fe_AS09 = N_over_H - Fe_over_H
        N_over_H_GA07 = 7.78
        N_over_Fe_GA07 = N_over_H_GA07 - Fe_over_H_GA07
        correction = N_over_Fe_GA07 - N_over_Fe_AS09
    elif element == "C":
        C_over_H = 8.43
        C_over_Fe_AS09 = C_over_H - Fe_over_H
        C_over_H_GA07 = 8.39
        C_over_Fe_GA07 = C_over_H_GA07 - Fe_over_H_GA07
        correction = C_over_Fe_GA07 - C_over_Fe_AS09

    el_FE += correction
    FE_H += Fe_over_H_GA07 - Fe_over_H

    select = np.where((GalR >= GalR_min) & (GalR <= GalR_max)
                      & (np.abs(Galz) <= Galz_max) & (np.abs(Galz) >= Galz_min))[0]

    contour = plot_contour(FE_H[select], el_FE[select], -1, 0.6, -0.5, 0.5, 'winter')
    legend_elements, _ = contour.legend_elements()
    return legend_elements


def make_figure(element, x, y, GalR, Galz, xlabel, ylabel, filename, sim_info, output_path):

    # Plot parameters
    params = {
        "font.size": 15,
        "font.family": "STIXGeneral",
        "text.usetex": False,
        "mathtext.fontset": "stix",
        "figure.figsize": (10, 6.5),
        "figure.subplot.left": 0.08,
        "figure.subplot.right": 0.97,
        "figure.subplot.bottom": 0.08,
        "figure.subplot.top": 0.93,
        "figure.subplot.wspace": 0.0,
        "figure.subplot.hspace": 0.0,
        "lines.markersize": 1,
        "lines.linewidth": 1,
        "figure.max_open_warning": 0,
    }
    plt.rcParams.update(params)

    ### Make plots
    plt.figure()

    ax = plt.subplot(3, 5, 1)
    plt.title('$0<R<3$ kpc')
    legend_elements_sim = make_plot(x, y, GalR, Galz, 1, 2, 0, 3)
    legend_elements_apogee = add_APOGEE(element, 1, 2, 0, 3)
    plt.ylabel(ylabel)

    plt.legend([legend_elements_sim[0],legend_elements_apogee[1]],
               [sim_info.simulation_name,"APOGEE"],
               labelspacing=0.1, handlelength=1.5,
               handletextpad=0.1, frameon=False, ncol=1, columnspacing=0.02)


    ####
    ax = plt.subplot(3, 5, 2)
    plt.title('$3<R<6$ kpc')
    make_plot(x, y, GalR, Galz, 1, 2, 3, 6)
    add_APOGEE(element, 1, 2, 3, 6)
    ax.get_yaxis().set_ticklabels([])

    ####
    ax = plt.subplot(3, 5, 3)
    plt.title('$6<R<9$ kpc')
    make_plot(x, y, GalR, Galz, 1, 2, 6, 9)
    add_APOGEE(element, 1, 2, 6, 9)
    ax.get_yaxis().set_ticklabels([])

    ####
    ax = plt.subplot(3, 5, 4)
    plt.title('$9<R<12$ kpc')
    make_plot(x, y, GalR, Galz, 1, 2, 9, 12)
    add_APOGEE(element, 1, 2, 9, 12)
    ax.get_yaxis().set_ticklabels([])

    ####
    ax = plt.subplot(3, 5, 5)
    plt.title('$12<R<15$ kpc')
    make_plot(x, y, GalR, Galz, 1, 2, 12, 15)
    add_APOGEE(element, 1, 2, 12, 15)
    ax.get_yaxis().set_ticklabels([])

    props = dict(boxstyle='round', fc='white', ec='black', alpha=1)
    ax.text(0.9, 0.9, '$1<|z|<2$ kpc', transform=ax.transAxes,
            fontsize=15, verticalalignment='top', bbox=props, rotation=270)

    ####
    ax = plt.subplot(3, 5, 6)
    make_plot(x, y, GalR, Galz, 0.5, 1.0, 0, 3)
    add_APOGEE(element, 0.5, 1.0, 0, 3)
    ax.get_xaxis().set_ticklabels([])
    plt.ylabel(ylabel)

    ####
    ax = plt.subplot(3, 5, 7)
    make_plot(x, y, GalR, Galz, 0.5, 1.0, 3, 6)
    add_APOGEE(element, 0.5, 1.0, 3, 6)
    ax.get_xaxis().set_ticklabels([])
    ax.get_yaxis().set_ticklabels([])

    ####
    ax = plt.subplot(3, 5, 8)
    make_plot(x, y, GalR, Galz, 0.5, 1.0, 6, 9)
    add_APOGEE(element, 0.5, 1.0, 6, 9)
    ax.get_xaxis().set_ticklabels([])
    ax.get_yaxis().set_ticklabels([])

    ####
    ax = plt.subplot(3, 5, 9)
    make_plot(x, y, GalR, Galz, 0.5, 1.0, 9, 12)
    add_APOGEE(element, 0.5, 1.0, 9, 12)
    ax.get_xaxis().set_ticklabels([])
    ax.get_yaxis().set_ticklabels([])

    ####
    ax = plt.subplot(3, 5, 10)
    make_plot(x, y, GalR, Galz, 0.5, 1.0, 12, 15)
    add_APOGEE(element, 0.5, 1.0, 12, 15)
    ax.get_xaxis().set_ticklabels([])
    ax.get_yaxis().set_ticklabels([])

    props = dict(boxstyle='round', fc='white', ec='black', alpha=1)
    ax.text(0.9, 0.9, '$0.5<|z|<1$ kpc', transform=ax.transAxes,
            fontsize=15, verticalalignment='top', bbox=props, rotation=270)

    ####
    ax = plt.subplot(3, 5, 11)
    make_plot(x, y, GalR, Galz, 0, 0.5, 0, 3)
    add_APOGEE(element, 0, 0.5, 0, 3)
    plt.ylabel(ylabel)
    plt.xlabel(xlabel)

    ####
    ax = plt.subplot(3, 5, 12)
    make_plot(x, y, GalR, Galz, 0, 0.5, 3, 6)
    add_APOGEE(element, 0, 0.5, 3, 6)
    ax.get_yaxis().set_ticklabels([])
    plt.xlabel(xlabel)

    ####
    ax = plt.subplot(3, 5, 13)
    make_plot(x, y, GalR, Galz, 0, 0.5, 6, 9)
    add_APOGEE(element, 0, 0.5, 6, 9)
    ax.get_yaxis().set_ticklabels([])
    plt.xlabel(xlabel)

    ####
    ax = plt.subplot(3, 5, 14)
    make_plot(x, y, GalR, Galz, 0, 0.5, 9, 12)
    add_APOGEE(element, 0, 0.5, 9, 12)
    ax.get_yaxis().set_ticklabels([])
    plt.xlabel(xlabel)

    ####
    ax = plt.subplot(3, 5, 15)
    make_plot(x, y, GalR, Galz, 0, 0.5, 12, 15)
    add_APOGEE(element, 0, 0.5, 12, 15)
    ax.get_yaxis().set_ticklabels([])
    plt.xlabel(xlabel)

    props = dict(boxstyle='round', fc='white', ec='black', alpha=1)
    ax.text(0.9, 0.9, '$0<|z|<0.5$ kpc', transform=ax.transAxes,
            fontsize=15, verticalalignment='top', bbox=props, rotation=270)

    plt.savefig(f"{output_path}/" + filename + "_" + sim_info.simulation_name + ".png", dpi=300)



def plot_stellar_abundances_radial_dependence(sim_info, output_path, abundance_data):

    # Look for abundance ratios from COLIBRE snaps:
    ratios, stars_age, stars_R, stars_z = \
        calculate_abundaces_from_MW_type_galaxies(sim_info)

    make_figure('MG', ratios['Fe_H'], ratios['Mg_Fe'], stars_R, stars_z,
                '[Fe/H]', '[Mg/Fe]', 'MgFe_comparison_radial_cut', sim_info, output_path)
    make_figure('O', ratios['Fe_H'], ratios['O_Fe'], stars_R, stars_z,
                '[Fe/H]', '[O/Fe]', 'OFe_comparison_radial_cut', sim_info, output_path)
    make_figure('N', ratios['Fe_H'], ratios['N_Fe'], stars_R, stars_z,
                '[Fe/H]', '[N/Fe]', 'NFe_comparison_radial_cut', sim_info, output_path)
    make_figure('C', ratios['Fe_H'], ratios['C_Fe'], stars_R, stars_z,
                '[Fe/H]', '[C/Fe]', 'CFe_comparison_radial_cut', sim_info, output_path)
