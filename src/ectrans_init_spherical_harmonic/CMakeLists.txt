enable_language(CXX)
ecbuild_add_library(TARGET ectrans_init_spherical_harmonic
    SOURCES
        ectrans_init_spherical_harmonic.h
        ectrans_init_spherical_harmonic.cc
        ectrans_init_spherical_harmonic_mod.F90
    PUBLIC_INCLUDES
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
)
target_compile_features( ectrans_init_spherical_harmonic PRIVATE cxx_std_11 )
