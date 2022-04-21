import numpy as np
from .html import add_web_section, PlotsInPipeline
from typing import List

def loadAbundancePlots(
    web, output_path: str, name_list: List[int]
):
    """
    @TODO Create separate .yaml config containing all necessary information about the plots
    """

    PlotsInWeb = PlotsInPipeline()

    num_sims = len(name_list)

    for i in range(num_sims):
        title = name_list[i]
        caption = "Carbon abundance [C/Fe] as a function of Iron abundance [Fe/H]. "
        caption += "Grey dots correspond to individual stars in MW-type haloes, solid line shows the median relation. "
        caption += "Contours corresponds to GALAH DR3 data (Buder+21) and are constructed from histogram of star counts "
        caption += "using a log scale with a minimum star count of 10."
        filename = "C_Fe_"+name_list[i]+".png"
        id = abs(hash("Carbon %i" %i))
        PlotsInWeb.load_plots(title, caption, filename, id)

    for i in range(num_sims):
        title = name_list[i]
        caption = "Carbon abundance [C/Fe] as a function of Iron abundance [Fe/H]. "
        caption += "Grey dots correspond to individual stars in MW-type haloes, solid line shows the median relation. "
        caption += "Contours corresponds to APOGEE data (Majewski et al. 2017, García Pérez et al. 2016) and are constructed from histogram of star counts "
        caption += "using a log scale with a minimum star count of 10."
        filename = "C_Fe_"+name_list[i]+"_2.png"
        id = abs(hash("Carbon 2 %i" %i))
        PlotsInWeb.load_plots(title, caption, filename, id)

    if num_sims > 1:
        title = "Comparison"
        caption = "Comparison between the [C/Fe]-[Fe/H] median relations from each simulation listed in this catalogue."
        filename = "C_Fe_comparison.png"
        id = abs(hash("Carbon comparison %i" % i))
        PlotsInWeb.load_plots(title, caption, filename, id)


    title = "Carbon"
    id = abs(hash("Carbon section"))
    plots = PlotsInWeb.plots_details
    caption = " "
    add_web_section(web, title, caption, id, plots)
    PlotsInWeb.reset_plots_list()

    for i in range(num_sims):
        title = name_list[i]
        caption = "Nitrogen abundance [N/Fe] as a function of Iron abundance [Fe/H]. "
        caption += "Grey dots correspond to individual stars within MW-type haloes, the solid line shows the median relation. "
        caption += "Contours corresponds to APOGEE data (Majewski et al. 2017, García Pérez et al. 2016) and are constructed from histogram of star counts "
        caption += "using a log scale with a minimum star count of 10."
        filename = "N_Fe_"+name_list[i]+".png"
        id = abs(hash("Nitrogen %i" %i))
        PlotsInWeb.load_plots(title, caption, filename, id)

    if num_sims > 1:
        title = "Comparison"
        caption = "Comparison between the [N/Fe]-[Fe/H] median relations from each simulation listed in this catalogue."
        filename = "N_Fe_comparison.png"
        id = abs(hash("Nitrogen_comparison %i" % i))
        PlotsInWeb.load_plots(title, caption, filename, id)

    title = "Nitrogen"
    id = abs(hash("Nitrogen section"))
    plots = PlotsInWeb.plots_details
    caption = " "
    add_web_section(web, title, caption, id, plots)
    PlotsInWeb.reset_plots_list()

    for i in range(num_sims):
        title = name_list[i]
        caption = "Oxygen abundance [O/Fe] as a function of Iron abundance [Fe/H]. "
        caption += "Left panel shows the stellar abundance from MW-type haloes, whereas right panel shows the abundance from satellites. "
        caption += "Grey dots correspond to individual stars, solid line shows the median relation. "
        caption += "Contours corresponds to GALAH DR3 data (Buder+21) and are constructed from histogram of star counts "
        caption += "using a log scale with a minimum star count of 10. "
        caption += "The observational data for MW, Carina, Fornax, Sculptor and Sagittarious, corresponds to a data compilation presented by Tolstoy, Hill & Tosi (2009)"
        filename = "O_Fe_"+name_list[i]+".png"
        id = abs(hash("Oxygen %i" %i))
        PlotsInWeb.load_plots(title, caption, filename, id)

    for i in range(num_sims):
        title = name_list[i]
        caption = "Oxygen abundance [O/Fe] as a function of Iron abundance [Fe/H]. "
        caption += "Both panels show the stellar abundance from MW-type haloes, with grey dots corresponding to individual stars and solid line showing the median relations. "
        caption += "The left panel shows contours corresponding to GALAH DR3 data (Buder+21), these are constructed from histogram of star counts "
        caption += "using a log scale with a minimum star count of 10. "
        caption += "The right panel shows APOGEE data (Majewski et al. 2017, García Pérez et al. 2016), where the contours are also constructed from histogram of star counts "
        caption += "using a log scale with a minimum star count of 10."
        caption += "The observational data for MW, Carina, Fornax, Sculptor and Sagittarious, corresponds to a data compilation presented by Tolstoy, Hill & Tosi (2009)"
        filename = "O_Fe_"+name_list[i]+"_2.png"
        id = abs(hash("Oxygen 2 %i" %i))
        PlotsInWeb.load_plots(title, caption, filename, id)

    if num_sims > 1:
        title = "Comparison"
        caption = "Comparison between the [O/Fe]-[Fe/H] median relations from each simulation listed in this catalogue."
        filename = "O_Fe_comparison.png"
        id = abs(hash("Oxygen comparison %i" % i))
        PlotsInWeb.load_plots(title, caption, filename, id)

    title = "Oxygen"
    id = abs(hash("Oxygen section"))
    plots = PlotsInWeb.plots_details
    caption = " "
    add_web_section(web, title, caption, id, plots)
    PlotsInWeb.reset_plots_list()

    for i in range(num_sims):
        title = name_list[i]
        caption = "Neon abundance [Ne/Fe] as a function of Iron abundance [Fe/H]. "
        caption += "Grey dots correspond to individual stars within MW-type haloes, the solid line shows the median relation. "
        filename = "Ne_Fe_"+name_list[i]+".png"
        id = abs(hash("Neon %i" %i))
        PlotsInWeb.load_plots(title, caption, filename, id)

    if num_sims > 1:
        title = "Comparison"
        caption = "Comparison between the [Ne/Fe]-[Fe/H] median relations from each simulation listed in this catalogue."
        filename = "Ne_Fe_comparison.png"
        id = abs(hash("Neon_comparison %i" % i))
        PlotsInWeb.load_plots(title, caption, filename, id)

    title = "Neon"
    id = abs(hash("Neon section"))
    plots = PlotsInWeb.plots_details
    caption = " "
    add_web_section(web, title, caption, id, plots)
    PlotsInWeb.reset_plots_list()

    for i in range(num_sims):
        title = name_list[i]
        caption = "Magensium abundance [Mg/Fe] as a function of Iron abundance [Fe/H]. "
        caption += "Left panel shows the stellar abundance from MW-type haloes, whereas right panel shows the abundance from satellites. "
        caption += "Grey dots correspond to individual stars, solid line shows the median relation. "
        caption += "Contours corresponds to GALAH DR3 data (Buder+21) and are constructed from histogram of star counts "
        caption += "using a log scale with a minimum star count of 10. "
        caption += "The observational data for MW, Carina, Fornax, Sculptor and Sagittarious, corresponds to a data compilation presented by Tolstoy, Hill & Tosi (2009)"
        filename = "Mg_Fe_"+name_list[i]+".png"
        id = abs(hash("Magnesium %i" %i))
        PlotsInWeb.load_plots(title, caption, filename, id)

    for i in range(num_sims):
        title = name_list[i]
        caption = "Magensium abundance [Mg/Fe] as a function of Iron abundance [Fe/H]. "
        caption += "Both panels show the stellar abundance from MW-type haloes, with grey dots corresponding to individual stars and solid line showing the median relations. "
        caption += "The left panel includes contours showing the GALAH DR3 data (Buder+21), these are constructed from histogram of star counts "
        caption += "using a log scale with a minimum star count of 10. "
        caption += "The right panel shows the APOGEE data (Majewski et al. 2017, García Pérez et al. 2016), where the contours are also constructed from histogram of star counts "
        caption += "using a log scale with a minimum star count of 10."
        caption += "The observational data for MW, Carina, Fornax, Sculptor and Sagittarious, corresponds to a data compilation presented by Tolstoy, Hill & Tosi (2009)"
        filename = "Mg_Fe_"+name_list[i]+"_2.png"
        id = abs(hash("Magnesium 2 %i" %i))
        PlotsInWeb.load_plots(title, caption, filename, id)

    if num_sims > 1:
        title = "Comparison"
        caption = "Comparison between the [Mg/Fe]-[Fe/H] median relations from each simulation listed in this catalogue."
        filename = "Mg_Fe_comparison.png"
        id = abs(hash("Magnesium comparison %i" % i))
        PlotsInWeb.load_plots(title, caption, filename, id)

    title = "Magnesium"
    id = abs(hash("Magnesium section"))
    plots = PlotsInWeb.plots_details
    caption = " "
    add_web_section(web, title, caption, id, plots)
    PlotsInWeb.reset_plots_list()

    for i in range(num_sims):
        title = name_list[i]
        caption = "Silicon abundance [Si/Fe] as a function of Iron abundance [Fe/H]. "
        caption += "Grey dots correspond to individual stars within MW-type haloes, the solid line shows the median relation. "
        caption += "Contours corresponds to GALAH DR3 data (Buder+21) and are constructed from histogram of star counts "
        caption += "using a log scale with a minimum star count of 10. "
        filename = "Si_Fe_"+name_list[i]+".png"
        id = abs(hash("Silicon %i" %i))
        PlotsInWeb.load_plots(title, caption, filename, id)

    if num_sims > 1:
        title = "Comparison"
        caption = "Comparison between the [Si/Fe]-[Fe/H] median relations from each simulation listed in this catalogue."
        filename = "Si_Fe_comparison.png"
        id = abs(hash("Silicon comparison %i" % i))
        PlotsInWeb.load_plots(title, caption, filename, id)

    title = "Silicon"
    id = abs(hash("Silicon section"))
    plots = PlotsInWeb.plots_details
    caption = " "
    add_web_section(web, title, caption, id, plots)
    PlotsInWeb.reset_plots_list()

    for i in range(num_sims):
        title = name_list[i]
        caption = "Strontium abundance [Sr/Fe] as a function of Iron abundance [Fe/H]. "
        caption += "Grey dots correspond to individual stars within MW-type haloes, the solid line shows the median relation. "
        caption += "Coloured symbols correspond to observational data from Roederer et al. (2014) and Spite et al. (2018)."
        filename = "Sr_Fe_"+name_list[i]+".png"
        id = abs(hash("Strontium %i" %i))
        PlotsInWeb.load_plots(title, caption, filename, id)

    if num_sims > 1:
        title = "Comparison"
        caption = "Comparison between the [Sr/Fe]-[Fe/H] median relations from each simulation listed in this catalogue."
        filename = "Sr_Fe_comparison.png"
        id = abs(hash("Strontium_comparison %i" % i))
        PlotsInWeb.load_plots(title, caption, filename, id)

    title = "Strontium"
    id = abs(hash("Strontium section"))
    plots = PlotsInWeb.plots_details
    caption = " "
    add_web_section(web, title, caption, id, plots)
    PlotsInWeb.reset_plots_list()

    for i in range(num_sims):
        title = name_list[i]
        caption = "Barium abundance [Ba/Fe] as a function of Iron abundance [Fe/H]. "
        caption += "Grey dots correspond to individual stars within MW-type haloes, the solid line shows the median relation. "
        caption += "Contours corresponds to GALAH DR3 data (Buder+21) and are constructed from histogram of star counts "
        caption += "using a log scale with a minimum star count of 10. "
        filename = "Ba_Fe_"+name_list[i]+".png"
        id = abs(hash("Barium %i" %i))
        PlotsInWeb.load_plots(title, caption, filename, id)

    if num_sims > 1:
        title = "Comparison"
        caption = "Comparison between the [Ba/Fe]-[Fe/H] median relations from each simulation listed in this catalogue."
        filename = "Ba_Fe_comparison.png"
        id = abs(hash("Barium comparison %i" % i))
        PlotsInWeb.load_plots(title, caption, filename, id)

    title = "Barium"
    id = abs(hash("Barium section"))
    plots = PlotsInWeb.plots_details
    caption = " "
    add_web_section(web, title, caption, id, plots)
    PlotsInWeb.reset_plots_list()

    for i in range(num_sims):
        title = name_list[i]
        caption = "Europium abundance [Eu/Fe] as a function of Iron abundance [Fe/H]. "
        caption += "Grey dots correspond to individual stars within MW-type haloes, the solid line shows the median relation. "
        caption += "Contours corresponds to GALAH DR3 data (Buder+21) and are constructed from histogram of star counts "
        caption += "using a log scale with a minimum star count of 10. "
        filename = "Eu_Fe_"+name_list[i]+".png"
        id = abs(hash("Europium %i" %i))
        PlotsInWeb.load_plots(title, caption, filename, id)

    if num_sims > 1:
        title = "Comparison"
        caption = "Comparison between the [Eu/Fe]-[Fe/H] median relations from each simulation listed in this catalogue."
        filename = "Eu_Fe_comparison.png"
        id = abs(hash("Europium comparison %i" % i))
        PlotsInWeb.load_plots(title, caption, filename, id)

    title = "Europium"
    id = abs(hash("Europium section"))
    plots = PlotsInWeb.plots_details
    caption = " "
    add_web_section(web, title, caption, id, plots)
    PlotsInWeb.reset_plots_list()

    title = "Stellar Mass - Z/Zsun relation (light-weighted r-band, 100 kpc aperture)"
    caption = "Comparison between the Stellar mass-Z/Zsun median relations from each simulation listed in this catalogue. "
    caption += "The values of Z/Zsun are obtained by calculated the light-weighted r-band average of the stars metallicity."
    filename = "Mstellar_Z_lw_comparison.png"
    id = abs(hash("Mstellar Z light weighted r band comparison %i" % i))
    PlotsInWeb.load_plots(title, caption, filename, id)

    title = "Stellar Mass - Z/Zsun relation (mass-weighted, 100 kpc aperture)"
    caption = "Comparison between the Stellar mass-Z/Zsun median relations from each simulation listed in this catalogue. "
    caption += "The values of Z/Zsun are obtained by calculated the mass-weighted average of the stars total metallicity."
    filename = "Mstellar_Z_mw_comparison.png"
    id = abs(hash("Mstellar Z mass weighted comparison %i" % i))
    PlotsInWeb.load_plots(title, caption, filename, id)

    title = "Stellar Metallicity"
    id = abs(hash("Stellar Metallicity section"))
    plots = PlotsInWeb.plots_details
    caption = " "
    add_web_section(web, title, caption, id, plots)
    PlotsInWeb.reset_plots_list()

    title = "Stellar Mass - [Fe/H] relation (mass-weight, 100 kpc aperture)"
    caption = "Comparison between the Stellar mass-[Fe/H] relations from each simulation listed in this catalogue. "
    caption += "The values of [Fe/H] are obtained by calculating &#x3c\Fe&#x3e_m = sum(Fe x mi)/sum(mi) (same for H). "
    caption += "Then doing log10(&#x3c\Fe&#x3e_m/&#x3cH&#x3e_m) and normalising by the corresponding solar abundances."
    filename = "Mstellar_FeH_mw_comparison.png"
    id = abs(hash("Mstellar FeH mass weighted comparison %i" % i))
    PlotsInWeb.load_plots(title, caption, filename, id)

    title = "Stellar Mass - [Fe/H] relation (mass-weight of ratio, 100 kpc aperture)"
    caption = "Comparison between the Stellar mass-[Fe/H] relations from each simulation listed in this catalogue. "
    caption += "The values of [Fe/H] are obtained by calculating &#x3c\Fe/H&#x3e_m = sum((Fe/H) x mi)/sum(mi). "
    caption += "Then doing log10(&#x3c\Fe/H&#x3e_m) and normalising by the corresponding solar abundances."
    filename = "Mstellar_FeH_mwr_comparison.png"
    id = abs(hash("Mstellar FeH ratio weighted comparison %i" % i))
    PlotsInWeb.load_plots(title, caption, filename, id)

    title = "Stellar Mass - [Fe/H] relation (mass-weight of log, 100 kpc aperture)"
    caption = "Comparison between the Stellar mass-[Fe/H] relations from each simulation listed in this catalogue. "
    caption += "The values of [Fe/H] are obtained by calculating log10(Fe/H), then weighting it by mass, "
    caption += "&#x3c log10(Fe/H)&#x3e_m, and normalising by the corresponding solar abundances."
    filename = "Mstellar_FeH_mwr2_comparison.png"
    id = abs(hash("Mstellar FeH ratio2 weighted comparison %i" % i))
    PlotsInWeb.load_plots(title, caption, filename, id)

    title = "Stellar Mass - [Fe/H] relation (light-weight r-band, 100 kpc aperture)"
    caption = "Comparison between the Stellar mass-[Fe/H] relations from each simulation listed in this catalogue. "
    caption += "The values of [Fe/H] are obtained by calculating &#x3c\Fe&#x3e_l = sum(Fe x li)/sum(li) (same for H, note "
    caption += "here li stands for r-band luminosity), "
    caption += "Then doing log10(&#x3c\Fe&#x3e_m/&#x3cH&#x3e_l) and normalising by the corresponding solar abundances."
    filename = "Mstellar_FeH_lw_comparison.png"
    id = abs(hash("Mstellar FeH light weighted comparison %i" % i))
    PlotsInWeb.load_plots(title, caption, filename, id)

    title = "Stellar Mass - [Fe/H] relation (light-weight r-band of ratio, 100 kpc aperture)"
    caption = "Comparison between the Stellar mass-[Fe/H] relations from each simulation listed in this catalogue. "
    caption += "The values of [Fe/H] are obtained by calculating &#x3c\Fe/H&#x3e_l = sum((Fe/H) x li)/sum(li). "
    caption += "Then doing log10(&#x3c\Fe/H&#x3e_l) and normalising by the corresponding solar abundances."
    filename = "Mstellar_FeH_lwr_comparison.png"
    id = abs(hash("Mstellar FeH ratio light weighted comparison %i" % i))
    PlotsInWeb.load_plots(title, caption, filename, id)

    title = "Stellar Mass - [Fe/H] relation (light-weight r-band of log, 100 kpc aperture)"
    caption = "Comparison between the Stellar mass-[Fe/H] relations from each simulation listed in this catalogue. "
    caption += "The values of [Fe/H] are obtained by calculating log10(Fe/H), then weighting it by r-band luminosity, "
    caption += "&#x3c log10(Fe/H)&#x3e_l, and normalising by the corresponding solar abundances."
    filename = "Mstellar_FeH_lwr2_comparison.png"
    id = abs(hash("Mstellar FeH ratio2 light weighted comparison %i" % i))
    PlotsInWeb.load_plots(title, caption, filename, id)

    title = "Stellar Mass - [Fe/H] relation (metallicity-weight, 100 kpc aperture)"
    caption = "Comparison between the Stellar mass-[Fe/H] relations from each simulation listed in this catalogue. "
    caption += "The values of [Fe/H] are obtained by calculating &#x3c\Fe&#x3e_Z = sum(Fe x Zi)/sum(Zi) (same for H, note "
    caption += "here Zi stands for metallicity). "
    caption += "Then doing log10(&#x3c\Fe&#x3e_m/&#x3cH&#x3e_Z) and normalising by the corresponding solar abundances."
    filename = "Mstellar_FeH_zw_comparison.png"
    id = abs(hash("Mstellar FeH metallicity weighted comparison %i" % i))
    PlotsInWeb.load_plots(title, caption, filename, id)

    title = "Stellar Mass - [Fe/H] relation (metallicity-weight of ratio, 100 kpc aperture)"
    caption = "Comparison between the Stellar mass-[Fe/H] relations from each simulation listed in this catalogue. "
    caption += "The values of [Fe/H] are obtained by calculating &#x3c\Fe/H&#x3e_Z = sum((Fe/H) x Zi)/sum(Zi). "
    caption += "Then doing log10(&#x3c\Fe/H&#x3e_Z) and normalising by the corresponding solar abundances."
    filename = "Mstellar_FeH_zwr_comparison.png"
    id = abs(hash("Mstellar FeH ratio metallicity weighted comparison %i" % i))
    PlotsInWeb.load_plots(title, caption, filename, id)

    title = "Stellar Mass - [Fe/H] relation (metallicity-weight of log, 100 kpc aperture)"
    caption = "Comparison between the Stellar mass-[Fe/H] relations from each simulation listed in this catalogue. "
    caption += "The values of [Fe/H] are obtained by calculating log10(Fe/H), then weighting it by metallicity, "
    caption += "&#x3c log10(Fe/H)&#x3e_Z, and normalising by the corresponding solar abundances."
    filename = "Mstellar_FeH_zwr2_comparison.png"
    id = abs(hash("Mstellar FeH ratio2 metallicity weighted comparison %i" % i))
    PlotsInWeb.load_plots(title, caption, filename, id)

    title = "Stellar Mass-[Fe/H]"
    id = abs(hash("Stellar Fe/H section"))
    plots = PlotsInWeb.plots_details
    caption = "In this section we explore three different ways of calculating weighted [Fe/H]: "
    caption += "(1) log10(&#x3c\Fe&#x3e_w / &#x3cH&#x3e_w), "
    caption += "(2) log10(&#x3c\Fe/H&#x3e_w), (3) &#x3clog10(Fe/H)&#x3e_w, where _w indicates weighted (meaning x_w = sum(wi x xi)/sum(wi)). "
    caption += "We weight [Fe/H] by mass, r-band luminosity and metallicity."
    add_web_section(web, title, caption, id, plots)
    PlotsInWeb.reset_plots_list()

    ##############################

    ## Oxygen ##

    title = "Stellar Mass - [O/Fe] relation (mass-weight, 100 kpc aperture)"
    caption = "Comparison between the Stellar mass-[O/Fe] relations from each simulation listed in this catalogue. "
    caption += "The values of [O/Fe] are obtained by calculating &#x3c\O&#x3e_m = sum(O x mi)/sum(mi) (same for Fe). "
    caption += "Then doing log10(&#x3c\O &#x3e_m/&#x3c Fe &#x3e_m) and normalising by the corresponding solar abundances."
    filename = "Mstellar_OFe_mw_comparison.png"
    id = abs(hash("Mstellar OFe mass weighted comparison %i" % i))
    PlotsInWeb.load_plots(title, caption, filename, id)

    title = "Stellar Mass - [O/Fe] relation (mass-weight of ratio, 100 kpc aperture)"
    caption = "Comparison between the Stellar mass-[O/Fe] relations from each simulation listed in this catalogue. "
    caption += "The values of [O/Fe] are obtained by calculating &#x3c\O/Fe &#x3e_m = sum((O/Fe) x mi)/sum(mi). "
    caption += "Then doing log10(&#x3c\O/Fe &#x3e_m) and normalising by the corresponding solar abundances."
    filename = "Mstellar_OFe_mwr_comparison.png"
    id = abs(hash("Mstellar OFe ratio weighted comparison %i" % i))
    PlotsInWeb.load_plots(title, caption, filename, id)

    title = "Stellar Mass - [O/Fe] relation (mass-weight of log, 100 kpc aperture)"
    caption = "Comparison between the Stellar mass-[O/Fe] relations from each simulation listed in this catalogue. "
    caption += "The values of [O/Fe] are obtained by calculating log10(O/Fe), then weighting it by mass, "
    caption += "&#x3c log10(O/Fe)&#x3e_m, and normalising by the corresponding solar abundances."
    filename = "Mstellar_OFe_mwr2_comparison.png"
    id = abs(hash("Mstellar OFe ratio2 weighted comparison %i" % i))
    PlotsInWeb.load_plots(title, caption, filename, id)

    title = "Stellar Mass - [O/Fe] relation (light-weight r-band, 100 kpc aperture)"
    caption = "Comparison between the Stellar mass-[O/Fe] relations from each simulation listed in this catalogue. "
    caption += "The values of [O/Fe] are obtained by calculating &#x3c\O&#x3e_l = sum(O x li)/sum(li) (same for Fe, note here li stands for r-band luminosity). "
    caption += "Then doing log10(&#x3c\O &#x3e_l/&#x3c Fe &#x3e_l) and normalising by the corresponding solar abundances."
    filename = "Mstellar_OFe_lw_comparison.png"
    id = abs(hash("Mstellar OFe light weighted comparison %i" % i))
    PlotsInWeb.load_plots(title, caption, filename, id)

    title = "Stellar Mass - [O/Fe] relation (light-weight r-band of ratio, 100 kpc aperture)"
    caption = "Comparison between the Stellar mass-[O/Fe] relations from each simulation listed in this catalogue. "
    caption += "The values of [O/Fe] are obtained by calculating &#x3c\O/Fe &#x3e_l = sum((O/Fe) x mi)/sum(mi). "
    caption += "Then doing log10(&#x3c\O/Fe &#x3e_l) and normalising by the corresponding solar abundances."
    filename = "Mstellar_OFe_lwr_comparison.png"
    id = abs(hash("Mstellar OFe ratio light weighted comparison %i" % i))
    PlotsInWeb.load_plots(title, caption, filename, id)

    title = "Stellar Mass - [O/Fe] relation (light-weight r-band of log, 100 kpc aperture)"
    caption = "Comparison between the Stellar mass-[O/Fe] relations from each simulation listed in this catalogue. "
    caption += "The values of [O/Fe] are obtained by calculating log10(O/Fe), then weighting it by r-band luminosity, "
    caption += "&#x3c log10(O/Fe)&#x3e_l, and normalising by the corresponding solar abundances."
    filename = "Mstellar_OFe_lwr2_comparison.png"
    id = abs(hash("Mstellar OFe ratio2 light weighted comparison %i" % i))
    PlotsInWeb.load_plots(title, caption, filename, id)

    title = "Stellar Mass - [O/Fe] relation (metallicity-weight, 100 kpc aperture)"
    caption = "Comparison between the Stellar mass-[O/Fe] relations from each simulation listed in this catalogue. "
    caption += "The values of [O/Fe] are obtained by calculating &#x3c\O&#x3e_Z = sum(O x Zi)/sum(Zi) (same for Fe, note here Zi stands for metallicity). "
    caption += "Then doing log10(&#x3c\O &#x3e_Z/&#x3c Fe &#x3e_Z) and normalising by the corresponding solar abundances."
    filename = "Mstellar_OFe_zw_comparison.png"
    id = abs(hash("Mstellar OFe metallicity weighted comparison %i" % i))
    PlotsInWeb.load_plots(title, caption, filename, id)

    title = "Stellar Mass - [O/Fe] relation (metallicity-weight of ratio, 100 kpc aperture)"
    caption = "Comparison between the Stellar mass-[O/Fe] relations from each simulation listed in this catalogue. "
    caption += "The values of [O/Fe] are obtained by calculating &#x3c\O/Fe &#x3e_Z = sum((O/Fe) x Zi)/sum(Zi). "
    caption += "Then doing log10(&#x3c\O/Fe &#x3e_Z) and normalising by the corresponding solar abundances."
    filename = "Mstellar_OFe_zwr_comparison.png"
    id = abs(hash("Mstellar OFe ratio metallicity weighted comparison %i" % i))
    PlotsInWeb.load_plots(title, caption, filename, id)

    title = "Stellar Mass - [O/Fe] relation (metallicity-weight of log, 100 kpc aperture)"
    caption = "Comparison between the Stellar mass-[O/Fe] relations from each simulation listed in this catalogue. "
    caption += "The values of [O/Fe] are obtained by calculating log10(O/Fe), then weighting it by metallicity, "
    caption += "&#x3c log10(O/Fe)&#x3e_Z, and normalising by the corresponding solar abundances."
    filename = "Mstellar_OFe_zwr2_comparison.png"
    id = abs(hash("Mstellar OFe ratio2 metallicity weighted comparison %i" % i))
    PlotsInWeb.load_plots(title, caption, filename, id)

    title = "Stellar Mass-[O/Fe]"
    id = abs(hash("Stellar Mass-[O/Fe] section"))
    plots = PlotsInWeb.plots_details
    caption = "In this section we explore three different ways of calculating weighted [O/Fe]: "
    caption += "(1) log10(&#x3c\O&#x3e_w / &#x3c\Fe&#x3e_w), "
    caption += "(2) log10(&#x3cO/Fe&#x3e_w), (3) &#x3clog10(O/Fe)&#x3e_w, where _w indicates weighted (meaning x_w = sum(wi x xi)/sum(wi)). "
    caption += "We weight [O/Fe] by mass, r-band luminosity and metallicity."
    add_web_section(web, title, caption, id, plots)
    PlotsInWeb.reset_plots_list()

    ## Magnesium ##

    title = "Stellar Mass - [Mg/Fe] relation (mass-weight, 100 kpc aperture)"
    caption = "Comparison between the Stellar mass-[Mg/Fe] relations from each simulation listed in this catalogue. "
    caption += "The values of [Mg/Fe] are obtained by calculating &#x3c\Mg &#x3e_m = sum(Mg x mi)/sum(mi) (same for Fe). "
    caption += "Then doing log10(&#x3c\Mg &#x3e_m/&#x3c Fe &#x3e_m) and normalising by the corresponding solar abundances."
    filename = "Mstellar_MgFe_mw_comparison.png"
    id = abs(hash("Mstellar MgFe mass weighted comparison %i" % i))
    PlotsInWeb.load_plots(title, caption, filename, id)

    title = "Stellar Mass - [Mg/Fe] relation (mass-weight of ratio, 100 kpc aperture)"
    caption = "Comparison between the Stellar mass-[Mg/Fe] relations from each simulation listed in this catalogue. "
    caption += "The values of [Mg/Fe] are obtained by calculating &#x3c\Mg/Fe &#x3e_m = sum((Mg/Fe) x mi)/sum(mi). "
    caption += "Then doing log10(&#x3c\Mg/Fe &#x3e_m) and normalising by the corresponding solar abundances."
    filename = "Mstellar_MgFe_mwr_comparison.png"
    id = abs(hash("Mstellar MgFe ratio weighted comparison %i" % i))
    PlotsInWeb.load_plots(title, caption, filename, id)

    title = "Stellar Mass - [Mg/Fe] relation (mass-weight of log, 100 kpc aperture)"
    caption = "Comparison between the Stellar mass-[Mg/Fe] relations from each simulation listed in this catalogue. "
    caption += "The values of [Mg/Fe] are obtained by calculating log10(Mg/Fe), then weighting it by mass, "
    caption += "&#x3c log10(Mg/Fe)&#x3e_m, and normalising by the corresponding solar abundances."
    filename = "Mstellar_MgFe_mwr2_comparison.png"
    id = abs(hash("Mstellar MgFe ratio2 weighted comparison %i" % i))
    PlotsInWeb.load_plots(title, caption, filename, id)

    title = "Stellar Mass - [Mg/Fe] relation (light-weight r-band, 100 kpc aperture)"
    caption = "Comparison between the Stellar mass-[Mg/Fe] relations from each simulation listed in this catalogue. "
    caption += "The values of [Mg/Fe] are obtained by calculating &#x3c\Mg &#x3e_l = sum(Mg x li)/sum(li) (same for Fe, note here li stands for r-band luminosity). "
    caption += "Then doing log10(&#x3c\Mg &#x3e_l/&#x3c Fe &#x3e_l) and normalising by the corresponding solar abundances."
    filename = "Mstellar_MgFe_lw_comparison.png"
    id = abs(hash("Mstellar MgFe light weighted comparison %i" % i))
    PlotsInWeb.load_plots(title, caption, filename, id)

    title = "Stellar Mass - [Mg/Fe] relation (light-weight r-band of ratio, 100 kpc aperture)"
    caption = "Comparison between the Stellar mass-[Mg/Fe] relations from each simulation listed in this catalogue. "
    caption += "The values of [Mg/Fe] are obtained by calculating &#x3c\Mg/Fe &#x3e_l = sum((Mg/Fe) x mi)/sum(mi). "
    caption += "Then doing log10(&#x3c\Mg/Fe &#x3e_l) and normalising by the corresponding solar abundances."
    filename = "Mstellar_MgFe_lwr_comparison.png"
    id = abs(hash("Mstellar MgFe ratio light weighted comparison %i" % i))
    PlotsInWeb.load_plots(title, caption, filename, id)

    title = "Stellar Mass - [Mg/Fe] relation (light-weight r-band of log, 100 kpc aperture)"
    caption = "Comparison between the Stellar mass-[Mg/Fe] relations from each simulation listed in this catalogue. "
    caption += "The values of [Mg/Fe] are obtained by calculating log10(Mg/Fe), then weighting it by r-band luminosity, "
    caption += "&#x3c log10(Mg/Fe)&#x3e_l, and normalising by the corresponding solar abundances."
    filename = "Mstellar_MgFe_lwr2_comparison.png"
    id = abs(hash("Mstellar MgFe ratio2 light weighted comparison %i" % i))
    PlotsInWeb.load_plots(title, caption, filename, id)

    title = "Stellar Mass - [Mg/Fe] relation (metallicity-weight, 100 kpc aperture)"
    caption = "Comparison between the Stellar mass-[Mg/Fe] relations from each simulation listed in this catalogue. "
    caption += "The values of [Mg/Fe] are obtained by calculating &#x3c\Mg &#x3e_Z = sum(Mg x Zi)/sum(Zi) (same for Fe, note here Zi stands for metallicity). "
    caption += "Then doing log10(&#x3c\Mg &#x3e_Z/&#x3c Fe &#x3e_Z) and normalising by the corresponding solar abundances."
    filename = "Mstellar_MgFe_zw_comparison.png"
    id = abs(hash("Mstellar MgFe metallicity weighted comparison %i" % i))
    PlotsInWeb.load_plots(title, caption, filename, id)

    title = "Stellar Mass - [Mg/Fe] relation (metallicity-weight of ratio, 100 kpc aperture)"
    caption = "Comparison between the Stellar mass-[Mg/Fe] relations from each simulation listed in this catalogue. "
    caption += "The values of [Mg/Fe] are obtained by calculating &#x3c Mg/Fe &#x3e_Z = sum((Mg/Fe) x Zi)/sum(Zi). "
    caption += "Then doing log10(&#x3c\Mg/Fe &#x3e_Z) and normalising by the corresponding solar abundances."
    filename = "Mstellar_MgFe_zwr_comparison.png"
    id = abs(hash("Mstellar MgFe ratio metallicity weighted comparison %i" % i))
    PlotsInWeb.load_plots(title, caption, filename, id)

    title = "Stellar Mass - [Mg/Fe] relation (metallicity-weight of log, 100 kpc aperture)"
    caption = "Comparison between the Stellar mass-[Mg/Fe] relations from each simulation listed in this catalogue. "
    caption += "The values of [Mg/Fe] are obtained by calculating log10(Mg/Fe), then weighting it by metallicity, "
    caption += "&#x3c log10(Mg/Fe)&#x3e_Z, and normalising by the corresponding solar abundances."
    filename = "Mstellar_MgFe_zwr2_comparison.png"
    id = abs(hash("Mstellar MgFe ratio2 metallicity weighted comparison %i" % i))
    PlotsInWeb.load_plots(title, caption, filename, id)

    title = "Stellar Mass-[Mg/Fe]"
    id = abs(hash("Stellar Mass-[Mg/Fe] section"))
    plots = PlotsInWeb.plots_details
    caption = "In this section we explore three different ways of calculating weighted [Mg/Fe]: "
    caption += "(1) log10(&#x3c\Mg&#x3e_w / &#x3c\Fe&#x3e_w), "
    caption += "(2) log10(&#x3c Mg/Fe&#x3e_w), (3) &#x3clog10(Mg/Fe)&#x3e_w, where _w indicates weighted (meaning x_w = sum(wi x xi)/sum(wi)). "
    caption += "We weight [Mg/Fe] by mass, r-band luminosity and metallicity."
    add_web_section(web, title, caption, id, plots)
    PlotsInWeb.reset_plots_list()

    ##############################

    for i in range(num_sims):
        title = name_list[i]
        caption = "Fe SNIa Fractions as a function of Iron abundance [Fe/H]. "
        caption += "Grey dots correspond to individual stars within MW-type haloes, the solid line shows the median relation. "
        filename = "FeSNIa_Fe_"+name_list[i]+".png"
        id = abs(hash("FeSNIa Fe %i" %i))
        PlotsInWeb.load_plots(title, caption, filename, id)

    if num_sims > 1:
        title = "Fe SNIa Fractions"
        caption = "Comparison between the fraction of Fe from SNIa from each simulation listed in this catalogue. "
        caption += "The curves show the median relations for the mass fraction of Fe from SNIa vs [Fe/H]."
        filename = "FeSNIa_Fe_comparison.png"
        id = abs(hash("FeSNIa Fe comparison %i" % i))
        PlotsInWeb.load_plots(title, caption, filename, id)

    title = "Fe(SNIa)/Fe"
    id = abs(hash("Fe(SNIa)/Fe section"))
    plots = PlotsInWeb.plots_details
    caption = " "
    add_web_section(web, title, caption, id, plots)
    PlotsInWeb.reset_plots_list()

    title = "SNIa Rates"
    caption = "Comparison between the SNIa rates from each simulation listed in this catalogue. "
    filename = "SNIa_rates_comparison.png"
    id = abs(hash("SNIa rates comparison %i" % i))
    PlotsInWeb.load_plots(title, caption, filename, id)

    title = "SFH"
    caption = "Comparison between the SFH from each simulation listed in this catalogue. "
    filename = "SFH_comparison.png"
    id = abs(hash("SFH comparison %i" % i))
    PlotsInWeb.load_plots(title, caption, filename, id)

    title = "SNIa Rates"
    id = abs(hash("SNIa rates section"))
    plots = PlotsInWeb.plots_details
    caption = " "
    add_web_section(web, title, caption, id, plots)
    PlotsInWeb.reset_plots_list()

    ##############################

    # title = "Stellar Mass - Z/Zsun relation (light-weighted i-band, 100 kpc aperture)"
    # caption = "Comparison between the Stellar mass-Z/Zsun median relations from each simulation listed in this catalogue. "
    # caption += "The values of Z/Zsun are obtained by calculated the light-weighted i-band average of the stars metallicity."
    # filename = "Mstellar_Z_light_weighted_i_band_comparison.png"
    # id = abs(hash("Mstellar Z light weighted i band comparison %i" % i))
    # PlotsInWeb.load_plots(title, caption, filename, id)
    #
    # title = "Stellar Mass - Z/Zsun relation (light-weighted z-band, 100 kpc aperture)"
    # caption = "Comparison between the Stellar mass-Z/Zsun median relations from each simulation listed in this catalogue. "
    # caption += "The values of Z/Zsun are obtained by calculated the light-weighted z-band average of the stars metallicity."
    # filename = "Mstellar_Z_light_weighted_z_band_comparison.png"
    # id = abs(hash("Mstellar Z light weighted z band comparison %i" % i))
    # PlotsInWeb.load_plots(title, caption, filename, id)
    #
    # title = "Stellar Mass - [Fe/H] relation (light-weighted i-band, 100 kpc aperture)"
    # caption = "Comparison between the Stellar mass-[Fe/H] median relations from each simulation listed in this catalogue. "
    # caption += "The values of [Fe/H] are obtained by calculated the light-weighted i-band average of [Fe/H]."
    # filename = "Mstellar_Fe_H_light_weighted_i_band_comparison.png"
    # id = abs(hash("Mstellar FeH light weighted i band comparison %i" % i))
    # PlotsInWeb.load_plots(title, caption, filename, id)
    #
    # title = "Stellar Mass - [Fe/H] relation (light-weighted z-band, 100 kpc aperture)"
    # caption = "Comparison between the Stellar mass-[Fe/H] median relations from each simulation listed in this catalogue. "
    # caption += "The values of [Fe/H] are obtained by calculated the light-weighted z-band average of [Fe/H]."
    # filename = "Mstellar_Fe_H_light_weighted_z_band_comparison.png"
    # id = abs(hash("Mstellar FeH light weighted z band comparison %i" % i))
    # PlotsInWeb.load_plots(title, caption, filename, id)
    #
    # title = "Extra plots"
    # id = abs(hash("Extra plots section"))
    # plots = PlotsInWeb.plots_details
    # caption = " "
    # add_web_section(web, title, caption, id, plots)
    # PlotsInWeb.reset_plots_list()
