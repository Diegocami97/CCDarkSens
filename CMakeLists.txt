cmake_minimum_required(VERSION 3.18)
project(ccdarksens LANGUAGES CXX)

# Match your ROOT toolchain
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ---- Dependencies ----
# Minuit2: for profile_minimizer "minuit". Omit Minuit2 from the list if your ROOT has no Minuit2.
find_package(ROOT REQUIRED COMPONENTS Core RIO Hist Graf Gpad Minuit2)
find_package(nlohmann_json CONFIG REQUIRED)
find_package(Threads REQUIRED)

# ---- Headers in repo ----
# (You can also use target_include_directories after targets; this is fine too)
include_directories(${PROJECT_SOURCE_DIR}/include)

# ---- Core library ----
add_library(ccdarksens_core
    # IO / config
    src/io/ConfigManager.cc
    src/io/RateTable.cc
    src/model/DMElectronModel.cc

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


    # Backgrounds
    src/backgrounds/PoissonDarkCurrent.cc
    src/backgrounds/BackgroundBuilder.cc

    # Rates
    src/rates/RateSource.cc

    # Stats
    src/stats/PoissonAsimovPLR.cc
    src/stats/ProfileLikelihood.cc
    src/stats/TestStatisticFactory.cc

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
add_executable(ccdarksens_sensitivity
    apps/ccdarksens_sensitivity.cc
    
)

# Make sure the app also sees nlohmann first (inherits via core, but explicit is fine)
target_include_directories(ccdarksens_sensitivity
    PUBLIC  ${NLOHMANN_INC_DIRS}
    SYSTEM PRIVATE ${ROOT_INCLUDE_DIRS}
)

target_link_libraries(ccdarksens_sensitivity
    PRIVATE ccdarksens_core
)


# Optional: put binary in project/build/
set_target_properties(ccdarksens_sensitivity PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/build
)

add_executable(ccdarksens_sensitivity_fast
  apps/ccdarksens_sensitivity_fast.cc
)
target_link_libraries(ccdarksens_sensitivity_fast PRIVATE ccdarksens_core)

# Optional: put binaries under build/
set_target_properties(ccdarksens_sensitivity_fast PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/build
)

add_custom_target(gen_qedark_rates
  COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/utils/qedark_generate_grid.py configs/qedark_generate_Si_heavy.json
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

add_executable(ccdarksens_check_dmelectron apps/ccdarksens_check_dmelectron.cc)
target_link_libraries(ccdarksens_check_dmelectron
  PRIVATE ccdarksens_core
          ROOT::Core ROOT::RIO ROOT::Hist ROOT::Graf ROOT::Gpad
          nlohmann_json::nlohmann_json)


add_executable(ccdarksens_scan_dmelectron_grid_dualspace 
  apps/ccdarksens_scan_dmelectron_grid_dualspace.cc)
target_link_libraries(ccdarksens_scan_dmelectron_grid_dualspace
  PRIVATE ccdarksens_core
          ROOT::Core ROOT::RIO ROOT::Hist ROOT::Graf ROOT::Gpad
          nlohmann_json::nlohmann_json)

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

add_executable(ccdarksens_scan_dmelectron_toy
  apps/ccdarksens_scan_dmelectron_toy.cc)
target_link_libraries(ccdarksens_scan_dmelectron_toy
  PRIVATE ccdarksens_core
          ROOT::Core ROOT::RIO ROOT::Hist ROOT::Graf ROOT::Gpad
          nlohmann_json::nlohmann_json)

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
