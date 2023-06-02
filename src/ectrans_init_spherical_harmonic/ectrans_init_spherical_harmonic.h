/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */


#pragma once

#ifdef __cplusplus
extern "C" {
#endif

/// \brief An analytic function that provides a spherical harmonic on a 2D sphere
///
/// \param n Total wave number
/// \param m Zonal wave number
/// \param lon Longitude in degrees
/// \param lat Latitude in degrees
/// \return spherical harmonic
double ectrans_init_spherical_harmonic(int n, int m, double lon, double lat, bool imag);
double ectrans_init_spherical_harmonic_hardcoded(int n, int m, double lon, double lat, bool imag);
double ectrans_init_spherical_harmonic_eastwest_derivative(int n, int m, double lon, double lat, bool imag);
double ectrans_init_spherical_harmonic_northsouth_derivative(int n, int m, double lon, double lat, bool imag);

#ifdef __cplusplus
}
#endif

