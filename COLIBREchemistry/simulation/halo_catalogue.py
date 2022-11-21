import numpy as np
import unyt
from velociraptor import load


class HaloCatalogue:
    """
    General class containing halo properties
    """

    def __init__(
        self, path_to_catalogue: str, 
        galaxy_min_stellar_mass: unyt.array.unyt_quantity,
        galaxy_min_gas_mass: unyt.array.unyt_quantity,
    ):
        """
        Parameters
        ----------
        path_to_catalogue: str
        Path to the catalogue with halo properties

        galaxy_min_stellar_mass: unyt.array.unyt_quantity
        Minimum stellar mass in units of Msun. Objects whose stellar mass is lower than this
        threshold are disregarded. Same for gas mass.
        """

        self.path_to_catalogue = path_to_catalogue

        # Load catalogue using velociraptor python library
        catalogue = load(self.path_to_catalogue)

        # Selecting central galaxies whose stellar mass is larger than
        # 'galaxy_min_stellar_mass'
        # mask = np.logical_and(
        #     catalogue.apertures.mass_star_30_kpc >= galaxy_min_stellar_mass,
        #     catalogue.structure_type.structuretype == 10,
        # )
        # They also need to contain at least one gas particle
        # mask = np.logical_and(
        #     mask, catalogue.apertures.mass_gas_30_kpc > unyt.unyt_quantity(0.0, "Msun")
        # )

        #mask = np.logical_and(
        #    catalogue.apertures.mass_star_30_kpc >= galaxy_min_stellar_mass,
        #    catalogue.apertures.mass_gas_30_kpc > unyt.unyt_quantity(1e7, "Msun")
        #)
  
        mask = np.logical_and(
            catalogue.apertures.mass_star_30_kpc >= galaxy_min_stellar_mass,
            catalogue.apertures.mass_gas_30_kpc > galaxy_min_gas_mass
        )

        # Compute the number of haloes following the selection mask
        self.number_of_haloes = mask.sum()

        # Log10 stellar mass in units of Msun
        self.log10_stellar_mass = np.log10(
            catalogue.apertures.mass_star_30_kpc.to("Msun").value[mask]
        )
        # Log10 gas mass in units of Msun
        self.log10_gas_mass = np.log10(
            catalogue.apertures.mass_gas_30_kpc.to("Msun").value[mask]
        )
        # Log10 halo mass in units of Msun
        self.log10_halo_mass = np.log10(
            catalogue.masses.mass_200crit.to("Msun").value[mask]
        )

        # Galaxy type, either central (=10) or satellite (>10)
        self.type = catalogue.structure_type.structuretype.value[mask]

        # Half mass radius in units of kpc (stars)
        self.half_mass_radius_star = catalogue.radii.r_halfmass_star.to("kpc").value[
            mask
        ]
        # Half mass radius in units of kpc (gas)
        self.half_mass_radius_gas = catalogue.radii.r_halfmass_gas.to("kpc").value[mask]

        # Star formation rate in units of Msun/yr
        self.sfr = (
            catalogue.apertures.sfr_gas_30_kpc.value[mask] * 10227144.8879616 / 1e9
        )

        # Metallicity of star-forming gas
        self.metallicity_gas_sfr = catalogue.apertures.zmet_gas_sf_30_kpc.value[mask]

        # Metallicity of all gas
        self.metallicity_gas = catalogue.apertures.zmet_gas_30_kpc.value[mask]

        self.metallicity = catalogue.apertures.zmet_star_30_kpc.value[mask]

        # Ids of haloes satisfying the selection criterion
        self.halo_ids = np.array([i for i in range(len(mask)) if mask[i] == True])

        self.xminpot = catalogue.positions.xcminpot.to("kpc").value[mask]
        self.yminpot = catalogue.positions.ycminpot.to("kpc").value[mask]
        self.zminpot = catalogue.positions.zcminpot.to("kpc").value[mask]

        self.vxminpot = catalogue.velocities.vxcminpot.to("km/s").value[mask]
        self.vyminpot = catalogue.velocities.vycminpot.to("km/s").value[mask]
        self.vzminpot = catalogue.velocities.vzcminpot.to("km/s").value[mask]

        self.kappa = None
