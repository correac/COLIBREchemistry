import numpy as np
import h5py
import matplotlib.pylab as plt
from matplotlib.pylab import rcParams

def load_satellites_data_Mg_Fe():
    # -------------------------------------------------------------------------------
    # alpha-enhancement (Mg/Fe), extracted manually from Tolstoy, Hill & Tosi (2009)
    # -------------------------------------------------------------------------------
    file = './plotter/obs_data/Fornax.txt'
    data = np.loadtxt(file)
    FeH_fornax = data[:, 0]
    MgFe_fornax = data[:, 1]

    file = './plotter/obs_data/Sculptor.txt'
    data = np.loadtxt(file)
    FeH_sculptor = data[:, 0]
    MgFe_sculptor = data[:, 1]

    file = './plotter/obs_data/Sagittarius.txt'
    data = np.loadtxt(file)
    FeH_sagittarius = data[:, 0]
    MgFe_sagittarius = data[:, 1]

    file = './plotter/obs_data/Carina.txt'
    data = np.loadtxt(file)
    FeH_carina = data[:, 0]
    MgFe_carina = data[:, 1]

    return FeH_fornax, MgFe_fornax, FeH_sculptor, MgFe_sculptor, \
           FeH_sagittarius, MgFe_sagittarius, FeH_carina, MgFe_carina

def load_satellites_data():
    # compute COLIBRE standard ratios
    Fe_over_H = 12. - 4.5
    Mg_over_H = 12. - 4.4
    O_over_H = 12. - 3.31
    Mg_over_Fe = Mg_over_H - Fe_over_H
    O_over_Fe = O_over_H - Fe_over_H

    # tabulate/compute the same ratios from Anders & Grevesse (1989)
    Fe_over_H_AG89 = 7.67
    Mg_over_H_AG89 = 7.58
    O_over_H_AG89 = 8.93

    # --
    Mg_over_Fe_AG89 = Mg_over_H_AG89 - Fe_over_H_AG89
    O_over_Fe_AG89 = O_over_H_AG89 - Fe_over_H_AG89

    ## I assume these works use Grevesser & Anders solar metallicity

    file = './plotter/obs_data/Letarte_2007.txt'
    data = np.loadtxt(file, skiprows=1)
    FeH_fornax = data[:, 0] + Fe_over_H_AG89 - Fe_over_H
    OFe_fornax = data[:, 4] + O_over_Fe_AG89 - O_over_Fe

    file = './plotter/obs_data/Sbordone_2007.txt'
    data = np.loadtxt(file, skiprows=1)
    FeH_sg = data[:, 0] + Fe_over_H_AG89 - Fe_over_H
    OFe_sg = data[:, 4] + O_over_Fe_AG89 - O_over_Fe

    file = './plotter/obs_data/Koch_2008.txt'
    data = np.loadtxt(file, skiprows=1)
    FeH_ca = data[:, 0] + Fe_over_H_AG89 - Fe_over_H
    OFe_ca = data[:, 4] + O_over_Fe_AG89 - O_over_Fe

    file = './plotter/obs_data/Geisler_2005.txt'
    data = np.loadtxt(file, skiprows=3)
    FeH_scu = data[:, 0] + Fe_over_H_AG89 - Fe_over_H
    OFe_scu = data[:, 4] - data[:, 0] + O_over_Fe_AG89 - O_over_Fe

    return FeH_fornax, OFe_fornax, FeH_sg, OFe_sg, FeH_ca, OFe_ca, FeH_scu, OFe_scu

def plot_MW_Satellites(element):

    if element == 'O':
        # Load Satellite data:
        FeH_fornax, OFe_fornax, FeH_sg, OFe_sg, FeH_ca, OFe_ca, FeH_scu, OFe_scu = load_satellites_data()

        plt.plot(FeH_ca, OFe_ca, '>', color='tab:purple', ms=4, label='Carina')
        plt.plot(FeH_scu, OFe_scu, '*', ms=4, color='tab:green', label='Sculptor')
        plt.plot(FeH_fornax, OFe_fornax, 'o', color='tab:orange', ms=4, label='Fornax')
        plt.plot(FeH_sg, OFe_sg, 'v', ms=4, color='crimson', label='Sagittarius')

    if element == 'Mg':

        FeH_fornax, MgFe_fornax, FeH_sculptor, MgFe_sculptor, \
        FeH_sagittarius, MgFe_sagittarius, FeH_carina, MgFe_carina = load_satellites_data_Mg_Fe()
        plt.plot(FeH_carina, MgFe_carina, 'o', color='crimson', ms=4, label='Carina')
        plt.plot(FeH_sculptor, MgFe_sculptor, '>', color='khaki', ms=4, label='Sculptor')
        plt.plot(FeH_fornax, MgFe_fornax, '<', color='royalblue', ms=4, label='Fornax')
        plt.plot(FeH_sagittarius, MgFe_sagittarius, '*', ms=4, color='lightblue', label='Sagittarius')


def load_MW_data():
    #compute COLIBRE standard ratios
    Fe_over_H = 12. - 4.5
    Mg_over_H = 12. - 4.4
    O_over_H = 12. - 3.31
    Mg_over_Fe = Mg_over_H - Fe_over_H
    O_over_Fe = O_over_H - Fe_over_H

    # tabulate/compute the same ratios from Anders & Grevesse (1989)
    Fe_over_H_AG89 = 7.67
    Mg_over_H_AG89 = 7.58
    O_over_H_AG89 = 8.93

    # --
    Mg_over_Fe_AG89 = Mg_over_H_AG89 - Fe_over_H_AG89
    O_over_Fe_AG89 = O_over_H_AG89 - Fe_over_H_AG89

    # MW data
    FeH_MW = []
    OFe_MW = []

    file = './plotter/obs_data/Koch_2008.txt'
    data = np.loadtxt(file, skiprows=3)
    FeH_koch = data[:, 1] + Fe_over_H_AG89 - Fe_over_H
    OH_koch = data[:, 2]
    OFe_koch = OH_koch - FeH_koch + O_over_Fe_AG89 - O_over_Fe

    FeH_MW = np.append(FeH_MW, FeH_koch)
    OFe_MW = np.append(OFe_MW, OFe_koch)

    file = './plotter/obs_data/Bai_2004.txt'
    data = np.loadtxt(file, skiprows=3, usecols=[1, 2])
    FeH_bai = data[:, 0] + Fe_over_H_AG89 - Fe_over_H
    OFe_bai = data[:, 1] + O_over_Fe_AG89 - O_over_Fe

    FeH_MW = np.append(FeH_MW, FeH_bai)
    OFe_MW = np.append(OFe_MW, OFe_bai)

    file = './plotter/obs_data/Cayrel_2004.txt'
    data = np.loadtxt(file, skiprows=18, usecols=[2, 6])
    FeH_cayrel = data[:, 0] + Fe_over_H_AG89 - Fe_over_H
    OFe_cayrel = data[:, 1] + O_over_Fe_AG89 - O_over_Fe

    FeH_MW = np.append(FeH_MW, FeH_cayrel)
    OFe_MW = np.append(OFe_MW, OFe_cayrel)

    file = './plotter/obs_data/Israelian_1998.txt'
    data = np.loadtxt(file, skiprows=3, usecols=[1, 3])
    FeH_isra = data[:, 0] + Fe_over_H_AG89 - Fe_over_H
    OFe_isra = data[:, 1] + O_over_Fe_AG89 - O_over_Fe

    FeH_MW = np.append(FeH_MW, FeH_isra)
    OFe_MW = np.append(OFe_MW, OFe_isra)

    file = './plotter/obs_data/Mishenina_1999.txt'
    data = np.loadtxt(file, skiprows=3, usecols=[1, 3])
    FeH_mish = data[:, 0] + Fe_over_H_AG89 - Fe_over_H
    OFe_mish = data[:, 1] + O_over_Fe_AG89 - O_over_Fe

    FeH_MW = np.append(FeH_MW, FeH_mish)
    OFe_MW = np.append(OFe_MW, OFe_mish)

    file = './plotter/obs_data/Zhang_Zhao_2005.txt'
    data = np.loadtxt(file, skiprows=3)
    FeH_zhang = data[:, 0] + Fe_over_H_AG89 - Fe_over_H
    OFe_zhang = data[:, 1] + O_over_Fe_AG89 - O_over_Fe

    FeH_MW = np.append(FeH_MW, FeH_zhang)
    OFe_MW = np.append(OFe_MW, OFe_zhang)

    return FeH_MW, OFe_MW

def plot_MW_data(element):

    if element == 'O':
        # Load MW data:
        FeH_MW, OFe_MW = load_MW_data()

        plt.plot(FeH_MW, OFe_MW, '+', color='tab:blue', ms=4, label='MW')

    if element == 'Mg':

        # Load MW data:
        FeH_MW, MgFe_MW = load_MW_data_with_Mg_Fe()
        plt.plot(FeH_MW, MgFe_MW, '+', color='orange', ms=4, label='MW')


def load_GALAH_data():
    file = "./plotter/obs_data/Buder21_data.hdf5"
    GALAH_dataobj = h5py.File(file, 'r')
    return GALAH_dataobj


def plot_GALAH_data(element):

    galah_data = load_GALAH_data()
    galah_edges = np.array(galah_data["abundance_bin_edges"])

    obs_plane = np.array(galah_data[f"{element}_enrichment_vs_Fe_abundance"]).T
    obs_plane[obs_plane < 10] = None

    contour = plt.contour(np.log10(obs_plane), origin='lower',
                          extent=[galah_edges[0], galah_edges[-1],
                                  galah_edges[0], galah_edges[-1]],
                          zorder=100,
                          cmap='winter',
                          linewidths=0.5)
    # contour.collections[0].set_label(['GALAH DR3'])


def load_strontium_data_Zhao():
    file = './plotter/obs_data/Zhao_2016.txt'
    # Roeder are calculated wrt to Lodder+ solar metallicity
    Fe_H_Sun = 7.50
    Eu_H_Sun = 0.53
    Ba_H_Sun = 2.18
    Sr_H_Sun = 2.9
    Sr_Fe_Sun = Sr_H_Sun - Fe_H_Sun

    # Asplund et al. (2009)
    Fe_H_Sun_As09 = 7.50
    Eu_H_Sun_As09 = 0.52
    Ba_H_Sun_As09 = 2.18
    Sr_H_Sun_As09 = 2.87

    Eu_Fe_Sun_As09 = Eu_H_Sun_As09 - Fe_H_Sun_As09
    Ba_Fe_Sun_As09 = Ba_H_Sun_As09 - Fe_H_Sun_As09
    Sr_Fe_Sun_As09 = Sr_H_Sun_As09 - Fe_H_Sun_As09

    data = np.loadtxt(file)
    FeH= data[:, 0]
    SrFe = data[:, 1] + Sr_Fe_Sun - Sr_Fe_Sun_As09
    return FeH, SrFe

def load_strontium_data_Roeder():
    file = './plotter/obs_data/Roeder_2014.txt'
    # Roeder are calculated wrt to Asplund+ solar metallicity
    data = np.loadtxt(file)
    FeH= data[:, 0]
    SrFe = data[:, 1]
    return FeH, SrFe

def load_strontium_data_Spite():
    file = './plotter/obs_data/Spite_2018.txt'
    # Roeder are calculated wrt to Lodder+ solar metallicity
    Fe_H_Sun = 7.50
    Eu_H_Sun = 0.53
    Ba_H_Sun = 2.18
    Sr_H_Sun = 2.9
    Sr_Fe_Sun = Sr_H_Sun - Fe_H_Sun

    # Asplund et al. (2009)
    Fe_H_Sun_As09 = 7.50
    Eu_H_Sun_As09 = 0.52
    Ba_H_Sun_As09 = 2.18
    Sr_H_Sun_As09 = 2.87

    Eu_Fe_Sun_As09 = Eu_H_Sun_As09 - Fe_H_Sun_As09
    Ba_Fe_Sun_As09 = Ba_H_Sun_As09 - Fe_H_Sun_As09
    Sr_Fe_Sun_As09 = Sr_H_Sun_As09 - Fe_H_Sun_As09

    data = np.loadtxt(file)
    FeH = data[:, 0]
    SrFe = data[:, 1] + Sr_Fe_Sun - Sr_Fe_Sun_As09

    return FeH, SrFe

def load_MW_data_with_Mg_Fe():
    file = './plotter/obs_data/MW.txt'
    data = np.loadtxt(file)
    FeH_mw = data[:, 0]
    MgFe_mw = data[:, 1]
    return FeH_mw, MgFe_mw

def plot_StrontiumObsData():

    FeH_Ro, SrFe_Ro = load_strontium_data_Roeder()
    FeH_Sp, SrFe_Sp = load_strontium_data_Spite()
    FeH_Zh, SrFe_Zh = load_strontium_data_Zhao()
    plt.plot(FeH_Ro, SrFe_Ro, '+', color='crimson', ms=4, label='Roederer et al. (2014)')
    plt.plot(FeH_Sp, SrFe_Sp, 'x', color='tab:blue', ms=4, label='Spite et al. (2018)')
    plt.plot(FeH_Zh, SrFe_Zh, 'o', color='tab:orange', ms=3, label='Zhao et al. (2016)')

def load_APOGEE_data():
    file = "./plotter/obs_data/APOGEE_data.hdf5"
    APOGEE_dataobj = h5py.File(file, 'r')
    return APOGEE_dataobj

def plot_APOGEE_data(element):

    apogee_data = load_APOGEE_data()
    FE_H = apogee_data['FE_H'][:]
    el_FE = apogee_data[f"{element}_FE"][:]

    # compute COLIBRE assumed abundances ( Asplund et al. 2009 )
    Fe_over_H = 7.5
    Mg_over_H = 7.6
    O_over_H = 8.69
    C_over_H = 8.43
    N_over_H = 7.83

    Mg_over_Fe_AS09 = Mg_over_H - Fe_over_H
    O_over_Fe_AS09 = O_over_H - Fe_over_H
    C_over_Fe_AS09 = C_over_H - Fe_over_H
    N_over_Fe_AS09 = N_over_H - Fe_over_H

    # tabulate/compute the same ratios from Grevesse, Asplund & Sauval (2007)
    Fe_over_H_GA07 = 7.45
    Mg_over_H_GA07 = 7.53
    O_over_H_GA07 = 8.66
    C_over_H_GA07 = 8.39
    N_over_H_GA07 = 7.78

    # --
    Mg_over_Fe_GA07 = Mg_over_H_GA07 - Fe_over_H_GA07
    O_over_Fe_GA07 = O_over_H_GA07 - Fe_over_H_GA07
    C_over_Fe_GA07 = C_over_H_GA07 - Fe_over_H_GA07
    N_over_Fe_GA07 = N_over_H_GA07 - Fe_over_H_GA07

    FE_H += Fe_over_H_GA07 - Fe_over_H
    if element == 'O': el_FE += O_over_Fe_GA07 - O_over_Fe_AS09
    if element == 'MG': el_FE += Mg_over_Fe_GA07 - Mg_over_Fe_AS09
    if element == 'N': el_FE += N_over_Fe_GA07 - N_over_Fe_AS09
    if element == 'C': el_FE += C_over_Fe_GA07 - C_over_Fe_AS09

    x = FE_H
    y = el_FE

    xmin = -3
    xmax = 1
    ymin = -1
    ymax = 1

    ngridx = 100
    ngridy = 50

    # Create grid values first.
    xi = np.linspace(xmin, xmax, ngridx)
    yi = np.linspace(ymin, ymax, ngridy)

    # Create a histogram
    h, xedges, yedges = np.histogram2d(x, y, bins=(xi, yi))
    xbins = xedges[:-1] + (xedges[1] - xedges[0]) / 2
    ybins = yedges[:-1] + (yedges[1] - yedges[0]) / 2

    z = h.T

    binsize = 0.25
    grid_min = np.log10(10)
    grid_max = np.log10(np.ceil(h.max()))
    levels = np.arange(grid_min, grid_max, binsize)
    levels = 10 ** levels

    #levels = np.arange(50, np.ceil(h.max()), 100)

    contour = plt.contour(xbins, ybins, z,
                          levels=levels, linewidths=0.5,
                          cmap='winter', zorder=100)
