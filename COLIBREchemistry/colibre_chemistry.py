"""
Welcome to COLIBRE chemistry, yet another plotting pipeline for the COLIBRE simulations
"""

from argumentparser import ArgumentParser
from plotter.stellar_abundances import plot_stellar_abundances
from plotter.compare_abundances import compare_stellar_abundances, compare_mass_metallicity_relations
from plotter.plot_mass_metallicity import compute_metallicity_relation
from simulation import simulation_data
from plotter.loadplots import loadAbundancePlots
from plotter.plot_SNIa_rates import read_SNIa_rates, plot_SNIa_rates
from plotter.plot_sfh import read_SFH, plot_SFH
from plotter import html
from time import time


def main(config: ArgumentParser):

    time_start = time()
    output_name_list = []
    web = None
    abundance_data = None
    metallicity_data = None
    SNIa_data = None
    SFH_data = None

    # Loop over simulation list
    for sim in range(config.number_of_inputs):

        # Fetch relevant input parameters from lists
        directory = config.directory_list[sim]
        snapshot = config.snapshot_list[sim]
        catalogue = config.catalogue_list[sim]
        sim_name = config.name_list[sim]

        # Load all data and save it in SimInfo class
        sim_info = simulation_data.SimInfo(
            directory=directory,
            snapshot=snapshot,
            catalogue=catalogue,
            name=sim_name,
            galaxy_min_stellar_mass=config.min_stellar_mass,
        )

        output_name_list.append(sim_info.simulation_name)

        # Make initial part of the webpage
        if sim == 0:
            web = html.make_web(sim_info.snapshot)
        elif web is not None:
            html.add_metadata_to_web(web, sim_info.snapshot)

        # Load luminosity tables
        simulation_data.SimInfo.load_photometry_grid()

        metallicity_data = compute_metallicity_relation(sim_info, metallicity_data)

        abundance_data = plot_stellar_abundances(sim_info, config.output_directory, abundance_data)

        SNIa_data = read_SNIa_rates(sim_info, SNIa_data)

        SFH_data = read_SFH(sim_info, SFH_data)


    if len(output_name_list) > 1: compare_stellar_abundances(abundance_data, output_name_list, config.output_directory)

    compare_mass_metallicity_relations(metallicity_data, output_name_list, config.output_directory)

    plot_SNIa_rates(SNIa_data, output_name_list, config.output_directory)
    plot_SFH(SFH_data, output_name_list, config.output_directory)

    loadAbundancePlots(web, config.output_directory, output_name_list)

    # Finish and output html file
    html.render_web(web, config.output_directory)

    # Compute how much time it took to run the script
    time_end = time()
    script_runtime = time_end - time_start
    m, s = divmod(script_runtime, 60)
    h, m = divmod(m, 60)
    print(f"The script was run in {h:.0f} hours, {m:.0f} minutes, and {s:.0f} seconds")

    return


if __name__ == "__main__":

    config_parameters = ArgumentParser()
    main(config_parameters)
