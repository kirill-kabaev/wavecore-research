cmake_minimum_required(VERSION 3.20)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
    message(FATAL_ERROR "In-source builds are not supported. Use: cmake -S . -B build")
endif()

project(
    WaveCore
    VERSION 0.1.0
    DESCRIPTION "Reproducible computational electromagnetics research"
    LANGUAGES CXX
)

function(wavecore_enable_warnings target_name)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(
            ${target_name}
            PRIVATE
                -Wall
                -Wextra
                -Wpedantic
        )
    endif()
endfunction()

add_library(
    wavecore
    STATIC
        src/helmholtz/assembly.cpp
        src/helmholtz/csr_matrix.cpp
        src/helmholtz/grid.cpp
        src/helmholtz/gmres.cpp
        src/helmholtz/version.cpp
)
add_library(WaveCore::wavecore ALIAS wavecore)

target_compile_features(wavecore PUBLIC cxx_std_20)
set_target_properties(wavecore PROPERTIES CXX_EXTENSIONS OFF)

target_include_directories(
    wavecore
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_compile_definitions(
    wavecore
    PRIVATE
        WAVECORE_VERSION="${PROJECT_VERSION}"
)

wavecore_enable_warnings(wavecore)

add_executable(
    helmholtz_2d
    src/apps/helmholtz_2d.cpp
)

target_link_libraries(helmholtz_2d PRIVATE WaveCore::wavecore)
set_target_properties(helmholtz_2d PROPERTIES CXX_EXTENSIONS OFF)
wavecore_enable_warnings(helmholtz_2d)

include(CTest)

if(BUILD_TESTING)
    add_executable(
        wavecore_basic_test
        tests/helmholtz/test_library.cpp
    )

    target_link_libraries(wavecore_basic_test PRIVATE WaveCore::wavecore)
    target_compile_definitions(
        wavecore_basic_test
        PRIVATE
            WAVECORE_EXPECTED_VERSION="${PROJECT_VERSION}"
    )
    set_target_properties(wavecore_basic_test PROPERTIES CXX_EXTENSIONS OFF)
    wavecore_enable_warnings(wavecore_basic_test)

    add_test(NAME wavecore.library_metadata COMMAND wavecore_basic_test)

    add_executable(
        wavecore_grid_test
        tests/helmholtz/test_grid.cpp
    )

    target_link_libraries(wavecore_grid_test PRIVATE WaveCore::wavecore)
    set_target_properties(wavecore_grid_test PROPERTIES CXX_EXTENSIONS OFF)
    wavecore_enable_warnings(wavecore_grid_test)

    add_test(NAME wavecore.grid COMMAND wavecore_grid_test)

    add_executable(
        wavecore_csr_matrix_test
        tests/helmholtz/test_csr_matrix.cpp
    )

    target_link_libraries(wavecore_csr_matrix_test PRIVATE WaveCore::wavecore)
    set_target_properties(wavecore_csr_matrix_test PROPERTIES CXX_EXTENSIONS OFF)
    wavecore_enable_warnings(wavecore_csr_matrix_test)

    add_test(NAME wavecore.csr_matrix COMMAND wavecore_csr_matrix_test)

    add_executable(
        wavecore_assembly_test
        tests/helmholtz/test_assembly.cpp
    )

    target_link_libraries(wavecore_assembly_test PRIVATE WaveCore::wavecore)
    set_target_properties(wavecore_assembly_test PROPERTIES CXX_EXTENSIONS OFF)
    wavecore_enable_warnings(wavecore_assembly_test)

    add_test(NAME wavecore.assembly COMMAND wavecore_assembly_test)

    add_executable(wavecore_gmres_test tests/helmholtz/test_gmres.cpp)
    target_link_libraries(wavecore_gmres_test PRIVATE WaveCore::wavecore)
    set_target_properties(wavecore_gmres_test PROPERTIES CXX_EXTENSIONS OFF)
    wavecore_enable_warnings(wavecore_gmres_test)
    add_test(NAME wavecore.gmres COMMAND wavecore_gmres_test)

    add_executable(wavecore_integration_test tests/helmholtz/test_integration.cpp)
    target_link_libraries(wavecore_integration_test PRIVATE WaveCore::wavecore)
    set_target_properties(wavecore_integration_test PROPERTIES CXX_EXTENSIONS OFF)
    wavecore_enable_warnings(wavecore_integration_test)
    foreach(integration_case IN ITEMS small medium altered restarted)
        add_test(NAME wavecore.integration.${integration_case}
                 COMMAND wavecore_integration_test ${integration_case})
    endforeach()

    add_executable(wavecore_convergence_test tests/helmholtz/test_convergence.cpp)
    target_link_libraries(wavecore_convergence_test PRIVATE WaveCore::wavecore)
    set_target_properties(wavecore_convergence_test PROPERTIES CXX_EXTENSIONS OFF)
    wavecore_enable_warnings(wavecore_convergence_test)
    add_test(NAME wavecore.convergence COMMAND wavecore_convergence_test)

    add_executable(wavecore_multimode_convergence_test
                   tests/helmholtz/test_multimode_convergence.cpp)
    target_link_libraries(wavecore_multimode_convergence_test PRIVATE WaveCore::wavecore)
    set_target_properties(wavecore_multimode_convergence_test PROPERTIES CXX_EXTENSIONS OFF)
    wavecore_enable_warnings(wavecore_multimode_convergence_test)
    add_test(NAME wavecore.convergence.multimode COMMAND wavecore_multimode_convergence_test)

    add_executable(wavecore_verification_checks_test tests/helmholtz/test_verification_checks.cpp)
    target_link_libraries(wavecore_verification_checks_test PRIVATE WaveCore::wavecore)
    set_target_properties(wavecore_verification_checks_test PROPERTIES CXX_EXTENSIONS OFF)
    wavecore_enable_warnings(wavecore_verification_checks_test)
    add_test(NAME wavecore.verification.acceptance COMMAND wavecore_verification_checks_test)
endif()
