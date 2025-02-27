cmake_minimum_required(VERSION 3.10)
project(peibos)

# Create an INTERFACE library for the header-only library
add_library(peibos INTERFACE)

# Correct include paths for build and install
target_include_directories(peibos INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>  # During build
    $<INSTALL_INTERFACE:include>  # After installation
)

# Install headers
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_PREFIX}/include)

# Install the CMake package
install(TARGETS peibos EXPORT peibosConfig)
install(EXPORT peibosConfig DESTINATION ${CMAKE_INSTALL_PREFIX}/share/peibos/cmake)

