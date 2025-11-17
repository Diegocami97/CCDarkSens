cmake_minimum_required(VERSION 3.18)
project(ccdarksens LANGUAGES CXX)

# Match your ROOT toolchain
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ---- Dependencies ----
find_package(ROOT REQUIRED COMPONENTS Core RIO Hist Graf Gpad)
find_package(nlohmann_json CONFIG REQUIRED)

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
    src/response/ClusterMC.cc

    # Backgrounds
    src/backgrounds/PoissonDarkCurrent.cc
    src/backgrounds/BackgroundBuilder.cc

    # Rates
    src/rates/RateSource.cc
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

add_executable(ccdarksens_sensitivity_cluster_mc
  apps/ccdarksens_sensitivity_cluster_mc.cc
)
target_link_libraries(ccdarksens_sensitivity_cluster_mc PRIVATE ccdarksens_core)

# Optional: put binaries under build/
set_target_properties(ccdarksens_sensitivity_fast PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/build
)
set_target_properties(ccdarksens_sensitivity_cluster_mc PROPERTIES
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

add_executable(ccdarksens_check_dmelectron_chain_cluster_mc
  apps/ccdarksens_check_dmelectron_chain_cluster_mc.cc)
target_link_libraries(ccdarksens_check_dmelectron_chain_cluster_mc
  PRIVATE ccdarksens_core
          ROOT::Core ROOT::RIO ROOT::Hist ROOT::Graf ROOT::Gpad
          nlohmann_json::nlohmann_json)



message(STATUS "CXX standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "ROOT include dirs: ${ROOT_INCLUDE_DIRS}")
