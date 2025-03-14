# ==================================================================
#  codac / basics example - cmake configuration file
# ==================================================================

cmake_minimum_required(VERSION 3.0.2)
project(peibos_example LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Adding Codac

# In case you installed Codac in a local directory, you need 
# to specify its path with the CMAKE_PREFIX_PATH option.
# set(CMAKE_PREFIX_PATH "~/codac/build_install")

  find_package(CODAC REQUIRED)
  message(STATUS "Found Codac version ${CODAC_VERSION}")

# Initializing IBEX

  ibex_init_common()

# Compilation

if(FAST_RELEASE)
  add_compile_definitions(FAST_RELEASE)
  message(STATUS "You are running Codac in fast release mode. (option -DCMAKE_BUILD_TYPE=Release is required)")
endif()

# Adding PEIBOS

find_package(peibos REQUIRED)

add_executable(${PROJECT_NAME} henon.cpp)
target_compile_options(${PROJECT_NAME} PUBLIC ${CODAC_CXX_FLAGS})
target_include_directories(${PROJECT_NAME} SYSTEM PUBLIC ${CODAC_INCLUDE_DIRS} )
target_link_libraries(${PROJECT_NAME} PUBLIC ${CODAC_LIBRARIES} peibos)