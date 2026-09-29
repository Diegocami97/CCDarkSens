# ============================================================================
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  CMakeLists.txt -- Build description of the framework: the ccdarksens_core
#  library (io, model, experiment, detector, response, backgrounds, stats and
#  scan helpers) and every app executable (scans, band, plots, validation and
#  example programs), linked against ROOT, Minuit2 and nlohmann_json.
# ============================================================================

cmake_minimum_required(VERSION 3.18)
project(ccdarksens LANGUAGES CXX)

# Match the ROOT toolchain in use
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ---- Dependencies ----
# Minuit2: for profile_minimizer "minuit". Omit Minuit2 from the list if the ROOT install has no Minuit2.
find_package(ROOT REQUIRED COMPONENTS Core RIO Hist Graf Gpad Minuit2)
find_package(nlohmann_json CONFIG REQUIRED)
find_package(Threads REQUIRED)

# ---- Headers in repo ----
# (target_include_directories after the targets would also work)
include_directories(${PROJECT_SOURCE_DIR}/include)

# ---- Core library ----
add_library(ccdarksens_core
    # IO / config
    src/io/ConfigManager.cc
    src/io/RateTable.cc
    src/model/DMElectronModel.cc
    src/model/DarkPhotonModel.cc
    src/model/MigdalModel.cc
    src/model/WimpNucleonModel.cc
    src/model/ModelFactory.cc

    # Experiment / detector
    src/experiment/ExperimentSetup.cc
    src/detector/Detector.cc

    # Response
    src/response/ChargeIonization.cc
    src/response/Diffusion.cc
    src/response/PatternEfficiency.cc
    src/response/DetectorResponsePipeline.cc
    src/response/PCDCalculator.cc
    src/response/PCDBasedResponse.cc
    src/response/PixelSimulator.cc
    src/response/PatternClassifier.cc
    src/response/EfficiencyMC.cc
    src/response/PatternImageGenerator.cc
    src/response/PatternRates.cc
    src/response/ChargeTransport.cc
    src/response/ClusterFitModel.cc
    src/response/ClusterFitEngine.cc
    src/response/NoiseTailCalibrator.cc
    src/response/ClusterFitMC.cc
    src/response/ClusterEnergyRates.cc
    src/response/PatternResponseFold.cc
    src/response/NeSpaceResponseFold.cc
    src/response/ClusterEnergyResponseFold.cc
    src/response/FlatBackgroundSpectrum.cc
    src/response/ResponseFactory.cc
    src/response/BackgroundFactory.cc
    src/response/BackgroundEfficiencyTable.cc


    # Backgrounds
    src/backgrounds/PoissonDarkCurrent.cc
    src/backgrounds/BackgroundBuilder.cc

    # Rates
    src/rates/RateSource.cc

    # Stats
    src/stats/PoissonAsimovPLR.cc
    src/stats/ProfileLikelihood.cc
    src/stats/TestStatisticFactory.cc

    # Scan utilities
    src/scan/ScanUtils.cc

)

# Ensure nlohmann json headers are seen BEFORE ROOT's vendored headers
get_target_property(NLOHMANN_INC_DIRS nlohmann_json::nlohmann_json INTERFACE_INCLUDE_DIRECTORIES)
target_include_directories(ccdarksens_core
    PUBLIC  ${NLOHMANN_INC_DIRS}                  # searched first by our targets
    SYSTEM PRIVATE ${ROOT_INCLUDE_DIRS}           # searched later; also silences warnings
)

target_link_libraries(ccdarksens_core
    PUBLIC ROOT::Hist ROOT::RIO ROOT::Core
           nlohmann_json::nlohmann_json
)
target_compile_definitions(ccdarksens_core PRIVATE CCDARKSENS_USE_MINUIT2)
target_link_libraries(ccdarksens_core PUBLIC ROOT::Minuit2)
message(STATUS "Profile likelihood: Minuit2 enabled (profile_minimizer: minuit)")

# ---- App ----
# NOTE: ccdarksens_sensitivity, ccdarksens_sensitivity_fast removed — synthetic flat-rate
# tutorial apps superseded by ccdarksens_example_one_point_pattern. Sources kept in apps/.

add_custom_target(gen_qedark_rates
  COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/utils/qedark_generate_grid.py configs/qedark_generate_Si_heavy.json
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

add_executable(ccdarksens_check_dmelectron apps/ccdarksens_check_dmelectron.cc)
target_link_libraries(ccdarksens_check_dmelectron
  PRIVATE ccdarksens_core
          ROOT::Core ROOT::RIO ROOT::Hist ROOT::Graf ROOT::Gpad
          nlohmann_json::nlohmann_json)


# NOTE: ccdarksens_scan_dmelectron_grid_dualspace removed — superseded by
# ccdarksens_scan_dmelectron_pattern which handles both ne and pattern spaces via config.

add_executable(ccdarksens_xcheck_ne_pcd
  apps/ccdarksens_xcheck_ne_pcd.cc)
target_link_libraries(ccdarksens_xcheck_ne_pcd
  PRIVATE ccdarksens_core
          ROOT::Core ROOT::RIO ROOT::Hist ROOT::Graf ROOT::Gpad
          nlohmann_json::nlohmann_json)

add_executable(ccdarksens_scan_dmelectron_pattern
  apps/ccdarksens_scan_dmelectron_pattern.cc)
target_link_libraries(ccdarksens_scan_dmelectron_pattern
  PRIVATE ccdarksens_core
          ROOT::Core ROOT::RIO ROOT::Hist ROOT::Graf ROOT::Gpad
          nlohmann_json::nlohmann_json)

add_executable(ccdarksens_scan_generic
  apps/ccdarksens_scan_generic.cc)
target_link_libraries(ccdarksens_scan_generic
  PRIVATE ccdarksens_core
          ROOT::Core ROOT::RIO ROOT::Hist ROOT::Graf ROOT::Gpad
          nlohmann_json::nlohmann_json)
set_target_properties(ccdarksens_scan_generic PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/build)

add_executable(ccdarksens_scan_srdm_pattern_csv
  apps/ccdarksens_scan_srdm_pattern_csv.cc)
target_link_libraries(ccdarksens_scan_srdm_pattern_csv
  PRIVATE ccdarksens_core
          ROOT::Core ROOT::RIO ROOT::Hist ROOT::Graf ROOT::Gpad ROOT::Tree
          nlohmann_json::nlohmann_json)

# Generic sensitivity-band orchestrator (subprocess driver of any scan binary
# satisfying the contract in docs/Plan_sensitivity_band_app.md §2).
add_executable(ccdarksens_band
  apps/ccdarksens_band.cc)
target_link_libraries(ccdarksens_band
  PRIVATE ROOT::Core ROOT::RIO ROOT::Hist ROOT::Graf ROOT::Tree
          nlohmann_json::nlohmann_json
          Threads::Threads)
set_target_properties(ccdarksens_band PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/build)

add_executable(ccdarksens_example_one_point_pattern
  apps/ccdarksens_example_one_point_pattern.cc)
target_link_libraries(ccdarksens_example_one_point_pattern
  PRIVATE ccdarksens_core
          ROOT::Core ROOT::RIO ROOT::Hist ROOT::Graf ROOT::Gpad
          nlohmann_json::nlohmann_json)
set_target_properties(ccdarksens_example_one_point_pattern PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/build)


add_executable(ccdarksens_scan_dmelectron_pattern_Edep
  apps/ccdarksens_scan_dmelectron_pattern_Edep.cc)
target_link_libraries(ccdarksens_scan_dmelectron_pattern_Edep
  PRIVATE ccdarksens_core
          ROOT::Core ROOT::RIO ROOT::Hist ROOT::Graf ROOT::Gpad
          nlohmann_json::nlohmann_json)

add_executable(ccdarksens_plot_dmelectron_limit
  apps/ccdarksens_plot_dmelectron_limit.cc
)
target_link_libraries(ccdarksens_plot_dmelectron_limit
  PRIVATE ccdarksens_core
    ROOT::Core ROOT::RIO ROOT::Hist ROOT::Graf ROOT::Gpad ROOT::Tree
)

# Fork of ccdarksens_plot_dmelectron_limit (NOT the same file -- the original
# stays untouched as the trusted reference) adding --wimp mode (WIMP-nucleon
# SI literature overlay). See docs/GenericScanApp_Design.md.
add_executable(ccdarksens_plot_limit
  apps/ccdarksens_plot_limit.cc
)
target_link_libraries(ccdarksens_plot_limit
  PRIVATE ccdarksens_core
    ROOT::Core ROOT::RIO ROOT::Hist ROOT::Graf ROOT::Gpad ROOT::Tree
)

# NOTE: ccdarksens_scan_dmelectron_toy removed — hardcoded 3-event tutorial app,
# no config file, no output file. Source kept in apps/.

add_executable(ccdarksens_pattern_background_eff
  apps/ccdarksens_pattern_background_eff.cc)
target_link_libraries(ccdarksens_pattern_background_eff
  PRIVATE ccdarksens_core
          ROOT::Core ROOT::RIO ROOT::Hist
          nlohmann_json::nlohmann_json)
set_target_properties(ccdarksens_pattern_background_eff PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/build)

# CSV → ROOT data converter (pattern counts, same style as S_pat)
add_executable(ccdarksens_csv_to_root_data
  apps/ccdarksens_csv_to_root_data.cc)
target_link_libraries(ccdarksens_csv_to_root_data
  PRIVATE ccdarksens_core
          ROOT::Core ROOT::RIO ROOT::Hist
          nlohmann_json::nlohmann_json)
set_target_properties(ccdarksens_csv_to_root_data PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/build)




 

add_executable(ccdarksens_validate_pattern_efficiency
  apps/ccdarksens_validate_pattern_efficiency.cc)
target_link_libraries(ccdarksens_validate_pattern_efficiency
  PRIVATE ccdarksens_core
          ROOT::Core ROOT::Hist
          nlohmann_json::nlohmann_json)
set_target_properties(ccdarksens_validate_pattern_efficiency PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/build)

add_executable(ccdarksens_validate_cluster_fit_engine
  apps/ccdarksens_validate_cluster_fit_engine.cc)
target_link_libraries(ccdarksens_validate_cluster_fit_engine
  PRIVATE ccdarksens_core
          ROOT::Core ROOT::Hist)
set_target_properties(ccdarksens_validate_cluster_fit_engine PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/build)

add_executable(ccdarksens_calibrate_noise_tail
  apps/ccdarksens_calibrate_noise_tail.cc)
target_link_libraries(ccdarksens_calibrate_noise_tail
  PRIVATE ccdarksens_core
          ROOT::Core ROOT::Hist)
set_target_properties(ccdarksens_calibrate_noise_tail PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/build)

add_executable(ccdarksens_build_cluster_fit_kernel
  apps/ccdarksens_build_cluster_fit_kernel.cc)
target_link_libraries(ccdarksens_build_cluster_fit_kernel
  PRIVATE ccdarksens_core
          ROOT::Core ROOT::Hist)
set_target_properties(ccdarksens_build_cluster_fit_kernel PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/build)

add_executable(ccdarksens_example_one_point_cluster
  apps/ccdarksens_example_one_point_cluster.cc)
target_link_libraries(ccdarksens_example_one_point_cluster
  PRIVATE ccdarksens_core
          ROOT::Core ROOT::Hist
          nlohmann_json::nlohmann_json)
set_target_properties(ccdarksens_example_one_point_cluster PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/build)

add_executable(ccdarksens_validate_wimp_nucleon_paper_repro
  apps/ccdarksens_validate_wimp_nucleon_paper_repro.cc)
target_link_libraries(ccdarksens_validate_wimp_nucleon_paper_repro
  PRIVATE ccdarksens_core
          ROOT::Core ROOT::Hist
          nlohmann_json::nlohmann_json)
set_target_properties(ccdarksens_validate_wimp_nucleon_paper_repro PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/build)

add_executable(ccdarksens_validate_response_factory
  apps/ccdarksens_validate_response_factory.cc)
target_link_libraries(ccdarksens_validate_response_factory
  PRIVATE ccdarksens_core
          ROOT::Core ROOT::Hist
          nlohmann_json::nlohmann_json)
set_target_properties(ccdarksens_validate_response_factory PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/build)

add_executable(ccdarksens_validate_background_efficiency
  apps/ccdarksens_validate_background_efficiency.cc)
target_link_libraries(ccdarksens_validate_background_efficiency
  PRIVATE ccdarksens_core
          ROOT::Core ROOT::Hist
          nlohmann_json::nlohmann_json)
set_target_properties(ccdarksens_validate_background_efficiency PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/build)

add_executable(ccdarksens_compute_brc
  apps/ccdarksens_compute_brc.cc)
target_link_libraries(ccdarksens_compute_brc
  PRIVATE ccdarksens_core
          nlohmann_json::nlohmann_json)
set_target_properties(ccdarksens_compute_brc PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/build)

message(STATUS "CXX standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "ROOT include dirs: ${ROOT_INCLUDE_DIRS}")
