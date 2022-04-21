#!/bin/bash

python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_04_COLIBRE/L025N188_NewYields_SNIILin8hi40_DTDExp2nu1p6_Fe05O07Mg07 \
                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_04_COLIBRE/L025N188_NewYields_SNIILin8hi40_DTDExp2nu1p6_Fe1O07Mg07 \
                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_04_COLIBRE/L025N188_NewYields_SNIILin8hi40_DTDExp2nu2_Fe05O07Mg07 \
                     -c halo_0036.properties halo_0036.properties halo_0036.properties \
                     -s colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 \
                     -n DTDnu1p6Fe05 DTDnu1p6Fe1 DTDnu2Fe05 \
                     -o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/AbundancePlots/2022_04_comparison/ComparisonDTDnuFeO07Mg07

#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_04_COLIBRE/L050N376_NewYields_SNIILin8hi40_DTDExp2nu1p6_ZdepB \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_04_COLIBRE/L025N188_NewYields_SNIILin8hi40_DTDExp2nu1p6_Fe05_ZdepBoost_04Zsun \
#                     -o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/AbundancePlots/2022_04_comparison/ConvergenceTest3 \
#                     -c halo_0036.properties halo_0036.properties \
#                     -s colibre_0036.hdf5 colibre_0036.hdf5 \
#                     -n L50DTDtau2nu1p6Zdep L25DTDtau2nu1p6Zdep \
#                     -m 5e9

#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_04_COLIBRE/L025N188_NewYields_SNIILin8hi40_DTDExp2nu1p6_Fe05_vBF \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_04_COLIBRE/L025N188_NewYields_SNIILin8hi40_DTDExp2nu1p6_Fe05_ZdepBoost_008Zsun \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_04_COLIBRE/L025N188_NewYields_SNIILin8hi40_DTDExp2nu1p6_Fe05_ZdepBoost_02Zsun \
#                     -o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/AbundancePlots/2022_04_comparison/L025SNII8to40Tau2nu1p6Fe05ZdepBoost \
#                     -c halo_0036.properties halo_0036.properties halo_0036.properties \
#                     -s colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 \
#                     -n SNII8to40Tau2nu1p6Fe05vB SNII8to40Tau2nu1p6Fe05Zdep008Zsun SNII8to40Tau2nu1p6Fe05vBZdep02Zsun \
#                     -m 5e8

#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_04_COLIBRE/L025N188_NewYields_SNIILin8hi40_DTDExp2nu1p6_vBF \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_04_COLIBRE/L025N188_NewYields_SNIILin8hi40_DTDExp2nu1p6_Fe1_ZdepBoost_008Zsun \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_04_COLIBRE/L025N188_NewYields_SNIILin8hi40_DTDExp2nu1p6_Fe1_ZdepBoost_02Zsun \
#                     -o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/AbundancePlots/2022_04_comparison/L025SNII8to40Tau2nu1p6Fe1ZdepBoost \
#                     -c halo_0036.properties halo_0036.properties halo_0036.properties \
#                     -s colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 \
#                     -n SNII8to40Tau2nu1p6Fe1vB SNII8to40Tau2nu1p6Fe1Zdep008Zsun SNII8to40Tau2nu1p6Fe1vBZdep02Zsun \
#                     -m 5e8

#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_03_COLIBRE/Test_Constant_Boost_factors/L025N188_NewYields_SNIILinlo8hi40_DTDExptau2nu2_Fe025 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_03_COLIBRE/L025N188_NewYields_SNIILinlo8hi40_DTDExptau2nu2_Fe05 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_03_COLIBRE/Test_Constant_Boost_factors/L025N188_NewYields_SNIILinlo8hi40_DTDExptau2nu2_Fe075 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_03_COLIBRE/Test_Constant_Boost_factors/L025N188_NewYields_SNIILinlo8hi40_DTDExptau2nu2_Fe1 \
#                     -o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/AbundancePlots/2022_03_comparisons/TestConstBoostFactorTau2nu2 \
#                     -c halo_0036.properties halo_0036.properties halo_0036.properties halo_0036.properties \
#                     -s colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 \
#                     -n DTDtau2nu2Fe025 DTDtau2nu2Fe05 DTDtau2nu2Fe075 DTDtau2nu2Fe1 \
#                     -m 5e8

#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_03_COLIBRE/Test_Metallicity_Dependent_Boost_Factor/L025N188_NewYields_SNIILinlo8hi40_DTDExptau2nu1p6_Fe05Zhigh1 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_03_COLIBRE/Test_Metallicity_Dependent_Boost_Factor/L025N188_NewYields_SNIILinlo8hi40_DTDExptau2nu1p6_Fe05Zhigh2 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_03_COLIBRE/Test_Metallicity_Dependent_Boost_Factor/L025N188_NewYields_SNIILinlo8hi40_DTDExptau2nu2_Fe05Zhigh1 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_03_COLIBRE/Test_Metallicity_Dependent_Boost_Factor/L025N188_NewYields_SNIILinlo8hi40_DTDExptau2nu2_Fe05Zhigh2 \
#                     -o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/AbundancePlots/2022_03_new_plots/L025NYSampleB \
#                     -c halo_0036.properties halo_0036.properties halo_0036.properties halo_0036.properties \
#                     -s colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 \
#                     -n DTDtau2nu1p6Fe05Z0p4Zsun DTDtau2nu1p6Fe05Z0p8Zsun DTDtau2nu2Fe05Z0p4Zsun DTDExptau2nu2Fe05Z0p8Zsun \
#                     -m 5e8


#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_03_COLIBRE/Test_Metallicity_Dependent_Boost_Factor/L025N188_NewYields_SNIILinlo8hi40_DTDExptau2nu1p6_BoostZhigh1 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_03_COLIBRE/Test_Metallicity_Dependent_Boost_Factor/L025N188_NewYields_SNIILinlo8hi40_DTDExptau2nu1p6_BoostZhigh2 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_03_COLIBRE/Test_Metallicity_Dependent_Boost_Factor/L025N188_NewYields_SNIILinlo8hi40_DTDExptau2nu1p6_BoostZhigh3 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo8hi40_DTDtau2nu1p6 \
#                     -o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/AbundancePlots/2022_03_comparisons/TestMetalOMgDepBoostFactorTau2nu1p6 \
#                     -c halo_0036.properties halo_0036.properties halo_0036.properties halo_0036.properties \
#                     -s colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 \
#                     -n DTDtau2nu1p6OMg05Z0p4Zsun DTDtau2nu1p6OMg05Z0p8Zsun DTDtau2nu1p6OMg05Z1p5Zsun DTDExptau2nu1p6 \
#                     -m 5e8

#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_03_COLIBRE/Test_Metallicity_Dependent_Boost_Factor/L025N188_NewYields_SNIILinlo8hi40_DTDExptau2nu1p6_Fe05Zhigh1 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_03_COLIBRE/Test_Metallicity_Dependent_Boost_Factor/L025N188_NewYields_SNIILinlo8hi40_DTDExptau2nu1p6_Fe05Zhigh2 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_03_COLIBRE/Test_Metallicity_Dependent_Boost_Factor/L025N188_NewYields_SNIILinlo8hi40_DTDExptau2nu1p6_Fe05Zhigh3 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_03_COLIBRE/L025N188_NewYields_SNIILinlo8hi40_DTDExptau2nu1p6_Fe05 \
#                     -o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/AbundancePlots/2022_03_comparisons/TestMetalDepBoostFactorTau2nu1p6 \
#                     -c halo_0036.properties halo_0036.properties halo_0036.properties halo_0036.properties \
#                     -s colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 \
#                     -n DTDtau2nu1p6Fe05Z0p4Zsun DTDtau2nu1p6Fe05Z0p8Zsun DTDtau2nu1p6Fe05Z1p5Zsun DTDExptau2nu1p6Fe05 \
#                     -m 5e8

#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_03_COLIBRE/Test_Metallicity_Dependent_Boost_Factor/L025N188_NewYields_SNIILinlo8hi40_DTDExptau2nu1p6_Fe02Zlow1 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_03_COLIBRE/Test_Metallicity_Dependent_Boost_Factor/L025N188_NewYields_SNIILinlo8hi40_DTDExptau2nu1p6_Fe02Zlow2 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_03_COLIBRE/Test_Metallicity_Dependent_Boost_Factor/L025N188_NewYields_SNIILinlo8hi40_DTDExptau2nu1p6_Fe02Zlow3 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_03_COLIBRE/L025N188_NewYields_SNIILinlo8hi40_DTDExptau2nu1p6_Fe05 \
#                     -o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/AbundancePlots/2022_03_comparisons/TestMetallicityDepBoostFactor1 \
#                     -c halo_0036.properties halo_0036.properties halo_0036.properties halo_0036.properties \
#                     -s colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 \
#                     -n DTDExptau2nu1p6Fe02Zlow1 DTDExptau2nu1p6Fe02Zlow2 DTDExptau2nu1p6Fe02Zlow3 DTDExptau2nu1p6Fe05 \
#                     -m 5e8

#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_03_COLIBRE/Test_Metallicity_Dependent_Boost_Factor/L025N188_NewYields_SNIILinlo8hi40_DTDExptau2nu1p6_Fe07FeZhigh1 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_03_COLIBRE/Test_Metallicity_Dependent_Boost_Factor/L025N188_NewYields_SNIILinlo8hi40_DTDExptau2nu1p6_Fe07FeZhigh2 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_03_COLIBRE/Test_Metallicity_Dependent_Boost_Factor/L025N188_NewYields_SNIILinlo8hi40_DTDExptau2nu1p6_Fe07FeZhigh3 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo8hi40_DTDtau2nu1p6 \
#                     -o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/AbundancePlots/2022_03_comparisons/TestMetalDepBoostFactorZhigh2 \
#                     -c halo_0036.properties halo_0036.properties halo_0036.properties halo_0036.properties \
#                     -s colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 \
#                     -n DTDtau2nu1p6Fe07Fe4Z0p4Zsun DTDtau2nu1p6Fe07Fe4Z0p6Zsun DTDtau2nu1p6Fe07Fe4Z0p8Zsun DTDtau2nu1p6Fe05 \
#                     -m 5e8

#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_03_COLIBRE/Test_Metallicity_Dependent_Boost_Factor/L025N188_NewYields_SNIILinlo8hi40_DTDExptau2nu1p6_Fe4Zhigh1 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_03_COLIBRE/Test_Metallicity_Dependent_Boost_Factor/L025N188_NewYields_SNIILinlo8hi40_DTDExptau2nu1p6_Fe4Zhigh2 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_03_COLIBRE/L025N188_NewYields_SNIILinlo8hi40_DTDExptau2nu1p6_Fe05 \
#                     -o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/AbundancePlots/2022_03_comparisons/TestMetallicityDepBoostFactor2 \
#                     -c halo_0036.properties halo_0036.properties halo_0036.properties \
#                     -s colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 \
#                     -n DTDExptau2nu1p6Fe4Zhigh1 DTDExptau2nu1p6Fe4Zhigh2 DTDExptau2nu1p6Fe05 \
#                     -m 5e8

#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_03_COLIBRE/L025N188_NewYields_SNIILinlo10hi40_DTDExptau2nu4 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_03_COLIBRE/L025N188_NewYields_SNIILinlo10hi40_DTDExptau4nu1_O05 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_03_COLIBRE/L025N188_NewYields_SNIILinlo10hi40_DTDExptau4nu1p5 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_03_COLIBRE/L025N188_NewYields_SNIILinlo10hi40_DTDExptau4nu1p5_O05Mg15 \
#                     -o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/AbundancePlots/2022_03_comparisons/DTDExptau4NewComparison \
#                     -c halo_0036.properties halo_0036.properties halo_0036.properties halo_0036.properties \
#                     -s colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 \
#                     -n DTDExptau2nu4 DTDExptau4nu1O05 DTDExptau4nu1p5 DTDtau4nu1p5O05Mg15 \
#                     -m 5e8


#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_03_COLIBRE/L025N188_NewYields_SNIILinlo8hi40_DTDExptau2nu2 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_03_COLIBRE/L025N188_NewYields_SNIILinlo8hi40_DTDExptau2nu2_Fe05 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_03_COLIBRE/L025N188_NewYields_SNIILinlo8hi40_DTDExptau2nu2_O15Mg15 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_03_COLIBRE/L025N188_NewYields_SNIILinlo8hi40_DTDExptau2nu2_Fe07 \
#                     -o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/AbundancePlots/2022_03_new_plots/L025NYSampleC \
#                     -c halo_0036.properties halo_0036.properties halo_0036.properties halo_0036.properties \
#                     -s colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 \
#                     -n DTDExptau2nu2 DTDExptau2nu2Fe05 DTDExptau2nu2O15Mg15 DTDtau2nu2Fe07 \
#                     -m 5e8

#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_03_COLIBRE/L025N188_NewYields_SNIILinlo8hi40_DTDExptau2nu4 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_03_COLIBRE/L025N188_NewYields_SNIILinlo10hi40_DTDExptau2nu4 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_03_COLIBRE/L025N188_NewYields_SNIILinlo8hi40_DTDExptau2nu3 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_03_COLIBRE/L025N188_NewYields_SNIILinlo8hi40_DTDExptau2nu3_Fe05 \
#                     -o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/AbundancePlots/2022_03_comparisons/DTDExpNewComparison \
#                     -c halo_0036.properties halo_0036.properties halo_0036.properties halo_0036.properties \
#                     -s colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 \
#                     -n 8to40DTDExptau2nu4 10to40DTDExptau2nu4 8to40DTDExptau2nu3 8to40DTDExptau2nu3Fe05 \
#                     -m 5e8

#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo10hi40_DTDtau1p5nu1p6_Fe05O05Mg05 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo10hi40_DTDtau2nu1p6 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo10hi40_DTDtau2nu1_O05Fe05Si05 \
#                     -o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/AbundancePlots/2022_03_comparisons/L025NYSNIIDTDExpTau2Comparison \
#                     -c halo_0036.properties halo_0036.properties halo_0036.properties \
#                     -s colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 \
#                     -n DTDtau1p5nu1p6Fe05O05 DTDtau2nu1p6 DTDtau2nu1Fe05O05 \
#                     -m 5e8

#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo8hi40_DTDPowerLaw \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo8hi40_DTDPowerLaw_Fe05 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo8hi40_DTDtau1p5nu1p6_Fe05O05Mg05 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo8hi40_DTDtau2nu1p6 \
#                     -o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/AbundancePlots/2022_02_comparisons/L025NYSNIIlo8hi40DTDComparison \
#                     -c halo_0036.properties halo_0036.properties halo_0036.properties halo_0036.properties \
#                     -s colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 \
#                     -n DTDPowerLawBeta05nu1 DTDtauPowerLawBeta05nu1Fe05 DTDtau1p5nu1p6Fe05O05Mg05 DTDtau2nu1p6 \
#                     -m 5e8


#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_StandardYields_SNIIlo8hi100_DTDtau2nu1p6 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo8hi100_DTDtau2nu1p6 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo8hi40_DTDtau2nu1p6 \
#                     -o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/AbundancePlots/2022_02_comparisons/L025CompareYields \
#                     -c halo_0036.properties halo_0036.properties halo_0036.properties \
#                     -s colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 \
#                     -n StandardYields NewYieldsSNIIlo8hi100 NewYieldsSNIIlo8hi40 \
#                     -m 5e8

#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo8hi40_DTDPowerLaw \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo8hi40_DTDPowerLaw_Fe05 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo10hi40_DTDPowerLaw \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo10hi40_DTDPowerLaw_Fe05 \
#                     -o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/AbundancePlots/2022_02_comparisons/L025NYSNIIMassPowerLawComparison \
#                     -c halo_0036.properties halo_0036.properties halo_0036.properties halo_0036.properties \
#                     -s colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 \
#                     -n lo8hi40DTDPowerLawBeta05nu1 lo8hi40DTDtauPowerLawBeta05nu1Fe05 lo10hi40DTDPowerLawBeta05nu1 lo10hi40DTDPowerLawBeta05nu1Fe05 \
#                     -m 5e8

#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo10hi40_DTDPowerLaw \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo10hi40_DTDPowerLaw_Fe05 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo10hi40_DTDPowerLawbeta08 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo10hi40_DTDPowerLawbeta08_Fe05 \
#                     -o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/AbundancePlots/2022_02_comparisons/L025NYSNIIlo10hi40DTDPowerLawComparison \
#                     -c halo_0036.properties halo_0036.properties halo_0036.properties halo_0036.properties \
#                     -s colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 \
#                     -n DTDPowerLawBeta05nu1 DTDtauPowerLawBeta05nu1Fe05 DTDPowerLawBeta08nu1 DTDPowerLawBeta08nu1Fe05 \
#                     -m 5e8

#python morpholopy.py -d /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo10hi40_DTDPowerLaw \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo10hi40_DTDPowerLaw_Fe05 \
#                        /cosma7/data/dp004/dc-corr1/SIMULATION_RUNS/2022_02_COLIBRE/L025N188_NewYields_SNIILinlo10hi40_DTDPowerLawbeta08 \
#                     -o /cosma/home/www/swift.dur.ac.uk/public_html/COLIBREPlots/ccorrea/AbundancePlots/2022_02_comparisons/L025NYSNIIlo10hi40DTDPowerLawComparison \
#                     -c halo_0036.properties halo_0036.properties halo_0036.properties \
#                     -s colibre_0036.hdf5 colibre_0036.hdf5 colibre_0036.hdf5 \
#                     -n DTDPowerLaw DTDtauPowerLawFe05 DTDPowerLawBeta08 \
#                     -m 5e8


