cmake_minimum_required(VERSION 3.20)
project(ccdarksens LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# -------- JSON (system or FetchContent fallback) --------
include(FetchContent)
find_package(nlohmann_json QUIET)
if(NOT nlohmann_json_FOUND)
  message(STATUS "nlohmann_json not found; fetching...")
  FetchContent_Declare(
    nlohmann_json
    GIT_REPOSITORY https://github.com/nlohmann/json.git
    GIT_TAG        v3.11.3
  )
  FetchContent_MakeAvailable(nlohmann_json)
endif()

# -------- ROOT --------
find_package(ROOT REQUIRED COMPONENTS Hist)
message(STATUS "Found ROOT: ${ROOT_USE_FILE}")

add_library(ccdarksens_core
  src/detector/Detector.cc
  src/experiment/ExperimentSetup.cc
  src/io/ConfigManager.cc
  src/response/ChargeIonization.cc
  src/response/PatternEfficiency.cc
  src/response/Diffusion.cc
  src/response/DetectorResponsePipeline.cc
  src/backgrounds/PoissonDarkCurrent.cc
  src/rates/RateSource.cc

)
target_include_directories(ccdarksens_core PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(ccdarksens_core PUBLIC nlohmann_json::nlohmann_json ROOT::Hist)

add_executable(ccdarksens_sensitivity apps/ccdarksens_sensitivity.cc)
target_link_libraries(ccdarksens_sensitivity PRIVATE ccdarksens_core)

# Optional: put executables in build/bin/
# set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
