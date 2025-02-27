cmake_minimum_required(VERSION 3.10)
project(peibos)

add_library(peibos INTERFACE)

target_include_directories(peibos INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

# Install headers
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_PREFIX}/include)

# Install the CMake package
install(TARGETS peibos EXPORT peibosConfig)
install(EXPORT peibosConfig DESTINATION ${CMAKE_INSTALL_PREFIX}/share/peibos/cmake)

