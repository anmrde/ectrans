/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <cmath>
#include <stdexcept>
#include<iostream>

#include "ectrans_init_spherical_harmonic.h"

namespace ectrans {
namespace init {
namespace {

constexpr double earth_radius = 6371229.;

static double factorial(double v) {
    if (v == 0) {
        return 1;
    }
    double result = v;
    while (--v > 0) {
        result *= v;
    }
    return result;
}

static double double_factorial(double x) {
    if (x == 0 || x == -1) {
        return 1;
    }

    double result = x;
    while ((x -= 2) > 0) {
        result *= x;
    }
    return result;
}

// Associated Legendre Polynomial
static double P(const int n, const int m, const double x) {
    // No recursive calculation needed
    if (n == m) {
        return (double_factorial(2 * m - 1) * std::pow(std::sqrt(1. - x * x), m));
    }

    if (n == m + 1) {
        return x * (2 * m + 1) * P(m, m, x);
    }

    // Formula 1
    //std::cout << "P: n: " << n << std::endl;
    return (x * (2 * n - 1) * P(n - 1, m, x) - (n + m - 1) * P(n - 2, m, x)) / (n - m);
}

static double K(const int n, const int m) {
    //When m is less than 0, multiply - 1 to pass in
    return std::sqrt(((2 * n + 1) * factorial(n - m)) / (factorial(n + m)));
}
}  // namespace

extern "C" {
double ectrans_init_spherical_harmonic(int n, int m, double lon, double lat, bool imag) {
    const int abs_m = std::abs(m);

    if (n < abs_m) {
	throw std::runtime_error("Error in ectrans_init_spherical_harmonic: assertion (n >= abs_m) failed");
    }

    double colat = M_PI_2 - lat;

    if (m == 0) {
        if (imag) {
            return 0.0;
        } else {
            return (K(n, m) * P(n, m, std::sin(lat)));
        }
    }

    if (m > 0) {
        if (imag) {
            return (-2 * K(n, m) * std::sin(m * lon) * P(n, m, std::sin(lat)));
        } else {
            return (2 * K(n, m) * std::cos(m * lon) * P(n, m, std::sin(lat)));
        }
    }

    // When m is less than 0, multiply - 1 in advance and send it to K
    if (imag) {
        return (-2 * K(n, abs_m) * std::cos(abs_m * lon) * P(n, abs_m, std::sin(lat)));
    } else {
        return (2 * K(n, abs_m) * std::sin(abs_m * lon) * P(n, abs_m, std::sin(lat)));
    }
}
}

extern "C" {
double ectrans_init_spherical_harmonic_eastwest_derivative(int n, int m, double lon, double lat, bool imag) {
    if (imag) {
        return (- m * ectrans_init_spherical_harmonic(n, m, lon, lat, !imag) / (earth_radius * cos(lat)));
    } else {
        return (m * ectrans_init_spherical_harmonic(n, m, lon, lat, !imag) / (earth_radius * cos(lat)));
    }
}
}

extern "C" {
double ectrans_init_spherical_harmonic_northsouth_derivative(int n, int m, double lon, double lat, bool imag) {

    double sinlat = std::sin(lat), coslat = std::cos(lat);

    const int abs_m = std::abs(m);

    if (n < abs_m) {
	throw std::runtime_error("Error in ectrans_init_spherical_harmonic: assertion (n >= abs_m) failed");
    }

    double colat = M_PI_2 - lat;

    double coeff_a = (n+1)*std::sqrt((n*n-m*m)/(4.0*n*n-1));
    //if(m==1 && n==1) std::cout << "coeff_a: " << coeff_a << std::endl;
    double coeff_b = -n*std::sqrt(((n+1)*(n+1)-m*m)/(4.0*(n+1)*(n+1)-1));
    //if(m==1 && n==1) std::cout << "coeff_b: " << coeff_b << std::endl;
    std::cout << "m=" << m << " n=" << n << std::endl;

    if (m == 0) {
        if (n > m) {
            if (imag) {
                return 0.0;
            } else {
                return (K(n-1, m) * coeff_a * P(n-1, m, std::sin(lat)) + K(n+1, m) * coeff_b * P(n+1, m, std::sin(lat))) / (earth_radius * coslat);
            }
        } else {
            if (imag) {
                return 0.0;
            } else {
                return (K(n+1, m) * coeff_b * P(n+1, m, std::sin(lat))) / (earth_radius * coslat);
            }
        }
    }

    if (m > 0) {
        if (n > m) {
            if (imag) {
                return -2 * std::sin(m * lon) * (K(n-1, m) * coeff_a * P(n-1, m, std::sin(lat)) + K(n+1, m) * coeff_b * P(n+1, m, std::sin(lat))) / (earth_radius * coslat);
            } else {
                return 2 * std::cos(m * lon) * (K(n-1, m) * coeff_a * P(n-1, m, std::sin(lat)) + K(n+1, m) * coeff_b * P(n+1, m, std::sin(lat))) / (earth_radius * coslat);
            }
        } else {
            if (imag) {
                return -2 * std::sin(m * lon) * (K(n+1, m) * coeff_b * P(n+1, m, std::sin(lat))) / (earth_radius * coslat);
            } else {
                return 2 * std::cos(m * lon) * (K(n+1, m) * coeff_b * P(n+1, m, std::sin(lat))) / (earth_radius * coslat);
            }
        }
    }

}
}

//-----------------------------------------------------------------------------
// Routine to compute the spherical harmonics analytically at one point
// (up to wave number 3)
//
// Author:
// Andreas Mueller *ECMWF*
//
double sphericalharmonics_analytic_point(
    const int n,         // total wave number (implemented so far for n<4
    const int m,         // zonal wave number (implemented so far for m<4, m<n
    const int imag,      // 0: test real part, 1: test imaginary part
    const double lon,    // longitude in radians
    const double lat,    // latitude in radians
    const int ivar_in,   // variable that is set to 1 for wave number n,m. 0:
                         // vorticity, 1: divergence, 2: scalar
    const int ivar_out)  // variable returned by this function. 0: u, 1: v, 2: scalar
{
    double latsin = std::sin(lat), latcos = std::cos(lat);
    double lonsin = std::sin(m * lon), loncos = std::cos(m * lon);
    double a = 6371229.;
    // Fourier part of the spherical harmonics:
    double rft = 1.;
    if (m > 0) {
        rft *= 2.;  // the famous factor 2 that noone really understands
    }
    if (imag == 0) {
        rft *= loncos;
    }
    else {
        rft *= -lonsin;
    }
    // Legendre part of the spherical harmonics (following
    // http://mathworld.wolfram.com/SphericalHarmonic.html
    // multiplied with -2*sqrt(pi) due to different normalization and different
    // coordinates):
    // (can also be computed on http://www.wolframalpha.com with:
    // LegendreP[n, m, x]/Sqrt[1/2*Integrate[LegendreP[n, m, y]^2, {y, -1, 1}]])
    // n, m need to be replaced by hand with the correct values
    // (otherwise the command will be too long for the free version of
    // wolframalpha)

    // scalar:
    if (ivar_in == 2) {
        if (ivar_out == 2) {
            if (m == 0 && n == 0) {
                return rft;
            }
            if (m == 0 && n == 1) {
                return std::sqrt(3.) * latsin * rft;
            }
            if (m == 0 && n == 2) {
                return std::sqrt(5.) / 2. * (3. * latsin * latsin - 1.) * rft;  // sign?
            }
            if (m == 0 && n == 3) {
                return std::sqrt(7.) / 2. * (5. * latsin * latsin - 3.) * latsin * rft;  // sign?
            }
            if (m == 1 && n == 1) {
                return std::sqrt(3. / 2.) * latcos * rft;  // sign?
            }
            if (m == 1 && n == 2) {
                return std::sqrt(15. / 2.) * latsin * latcos * rft;  // sign?
            }
            if (m == 1 && n == 3) {
                return std::sqrt(21.) / 4. * latcos * (5. * latsin * latsin - 1.) * rft;  // sign?
            }
            if (m == 2 && n == 2) {
                return std::sqrt(15. / 2.) / 2. * latcos * latcos * rft;
            }
            if (m == 2 && n == 3) {
                return std::sqrt(105. / 2.) / 2. * latcos * latcos * latsin * rft;
            }
            if (m == 3 && n == 3) {
                return std::sqrt(35.) / 4. * latcos * latcos * latcos * rft;  // sign?
            }
            if (m == 4 && n == 4) {
                return (3 * std::sqrt(17.5) * std::pow(latcos, 4)) / 8. * rft;
            }
            if (m == 5 && n == 5) {
                return (3 * std::sqrt(77) * std::pow(latcos, 5)) / 16. * rft;
            }
            if (m == 6 && n == 6) {
                return (std::sqrt(3003) * std::pow(latcos, 6)) / 32. * rft;
            }
            if (m == 7 && n == 7) {
                return (3 * std::sqrt(357.5) * std::pow(latcos, 7)) / 32. * rft;
            }
            if (m == 8 && n == 8) {
                return (3 * std::sqrt(6077.5) * std::pow(latcos, 8)) / 128. * rft;
            }
            if (m == 9 && n == 9) {
                return (std::sqrt(230945) * std::pow(latcos, 9)) / 256. * rft;
            }
            if (m == 10 && n == 10) {
                return (std::sqrt(969969) * std::pow(latcos, 10)) / 512. * rft;
            }
            if (m == 11 && n == 11) {
                return (std::sqrt(1.0140585e6) * std::pow(latcos, 11)) / 512. * rft;
            }
            if (m == 12 && n == 12) {
                return (5 * std::sqrt(676039) * std::pow(latcos, 12)) / 2048. * rft;
            }
            if (m == 13 && n == 13) {
                return (15 * std::sqrt(78004.5) * std::pow(latcos, 13)) / 2048. * rft;
            }
            if (m == 14 && n == 14) {
                return (15 * std::sqrt(323161.5) * std::pow(latcos, 14)) / 4096. * rft;
            }
            if (m == 15 && n == 15) {
                return (3 * std::sqrt(33393355) * std::pow(latcos, 15)) / 8192. * rft;
            }
            if (m == 16 && n == 16) {
                return (3 * std::sqrt(5.509903575e8) * std::pow(latcos, 16)) / 32768. * rft;
            }
            if (m == 17 && n == 17) {
                return (15 * std::sqrt(90751353) * std::pow(latcos, 17)) / 65536. * rft;
            }
            if (m == 18 && n == 18) {
                return (5 * std::sqrt(3357800061) * std::pow(latcos, 18)) / 131072. * rft;
            }
            if (m == 19 && n == 19) {
                return (15 * std::sqrt(3.829070245e8) * std::pow(latcos, 19)) / 131072. * rft;
            }
            if (m == 20 && n == 20) {
                return (3 * std::sqrt(156991880045) * std::pow(latcos, 20)) / 524288. * rft;
            }
            if (m == 21 && n == 21) {
                return (std::sqrt(1.4465680375575e12) * std::pow(latcos, 21)) / 524288. * rft;
            }
            if (m == 22 && n == 22) {
                return (15 * std::sqrt(2.63012370465e10) * std::pow(latcos, 22)) / 1.048576e6 * rft;
            }
            if (m == 23 && n == 23) {
                return (15 * std::sqrt(107492012277) * std::pow(latcos, 23)) / 2.097152e6 * rft;
            }
            if (m == 24 && n == 24) {
                return (105 * std::sqrt(35830670759) * std::pow(latcos, 24)) / 8.388608e6 * rft;
            }
            if (m == 25 && n == 25) {
                return (21 * std::sqrt(9.136821043545e11) * std::pow(latcos, 25)) / 8.388608e6 * rft;
            }
            if (m == 26 && n == 26) {
                return (21 * std::sqrt(3.7250116562145e12) * std::pow(latcos, 26)) / 1.6777216e7 * rft;
            }
            if (m == 27 && n == 27) {
                return (7 * std::sqrt(136583760727865.) * std::pow(latcos, 27)) / 3.3554432e7 * rft;
            }
            if (m == 28 && n == 28) {
                return (std::sqrt(2.7248460265209068e16) * std::pow(latcos, 28)) / 6.7108864e7 * rft;
            }
            if (m == 29 && n == 29) {
                return (std::sqrt(110873045217057585.) * std::pow(latcos, 29)) / 1.34217728e8 * rft;
            }
            if (m == 30 && n == 30) {
                return (std::sqrt(450883717216034179.) * std::pow(latcos, 30)) / 2.68435456e8 * rft;
            }
            if (m == 31 && n == 31) {
                return (21 * std::sqrt(1.0389025742304935e15) * std::pow(latcos, 31)) / 2.68435456e8 * rft;
            }
            if (m == 32 && n == 32) {
                return (21 * std::sqrt(6.752866732498208e16) * std::pow(latcos, 32)) / 2.147483648e9 * rft;
            }
            if (m == 33 && n == 33) {
                return (7 * std::sqrt(2467865842240254105.) * std::pow(latcos, 33)) / 4.294967296e9 * rft;
            }
            if (m == 34 && n == 34) {
                return (21 * std::sqrt(1112959105324036165.) * std::pow(latcos, 34)) / 8.589934592e9 * rft;
            }
            if (m == 35 && n == 35) {
                return (3 * std::sqrt(5.53140675346046e19) * std::pow(latcos, 35)) / 8.589934592e9 * rft;
            }
            if (m == 36 && n == 36) {
                return (std::sqrt(8075853860052271220473.) * std::pow(latcos, 36)) / 3.4359738368e10 * rft;
            }
            if (m == 37 && n == 37) {
                return (5 * std::sqrt(3.2739948081292994e20) * std::pow(latcos, 37)) / 3.4359738368e10 * rft;
            }
            if (m == 38 && n == 38) {
                return (35 * std::sqrt(2.707815254843781e19) * std::pow(latcos, 38)) / 6.8719476736e10 * rft;
            }
            if (m == 39 && n == 39) {
                return (35 * std::sqrt(109701233401363445369.) * std::pow(latcos, 39)) / 1.37438953472e11 * rft;
            }
            if (m == 40 && n == 40) {
                return (63 * std::sqrt(548506167006817226845.) * std::pow(latcos, 40)) / 5.49755813888e11 * rft;
            }
            if (m == 41 && n == 41) {
                return (63 * std::sqrt(5.551952666044613e20) * std::pow(latcos, 41)) / 5.49755813888e11 * rft;
            }
            if (m == 42 && n == 42) {
                return (15 * std::sqrt(3.964094203555854e22) * std::pow(latcos, 42)) / 1.099511627776e12 * rft;
            }
            if (m == 43 && n == 43) {
                return (45 * std::sqrt(17823059209786010066579.) * std::pow(latcos, 43)) / 2.199023255552e12 * rft;
            }
            if (m == 44 && n == 44) {
                return (45 * std::sqrt(7.210237589413431e22) * std::pow(latcos, 44)) / 4.398046511104e12 * rft;
            }
            if (m == 45 && n == 45) {
                return (21 * std::sqrt(1339044123748208678378695.) * std::pow(latcos, 45)) / 8.796093022208e12 * rft;
            }
            //    return
            //    std::pow(latcos,45)*rft*21.*std::sqrt(1339044123748208678378695.)/8796093022208.;
            //    // sign?
        }
        else {
            return 0.;
        }
    }

    // for the remainder the factor 2 from rft is already included in the
    // formulas:

    // vorticity:
    if (ivar_in == 0) {
        if (ivar_out == 0) {  // u:
            if (m == 0 && n == 0) {
                return 0.;
            }
            if (m == 0 && n == 1) {
                if (imag == 0) {
                    return std::sqrt(3.) * a / 2. * latcos;
                }
                else {
                    return 0.;
                }
            }
            if (m == 1 && n == 1) {
                if (imag == 0) {
                    return -a * std::sqrt(3. / 2.) * loncos * latsin;
                }
                else {
                    return a * std::sqrt(3. / 2.) * lonsin * latsin;
                }
            }
        }
        else if (ivar_out == 1) {  // v:
            if (m == 0 && n == 0) {
                return 0.;
            }
            if (m == 0 && n == 1) {
                return 0.;
            }
            if (m == 1 && n == 1) {
                if (imag == 0) {
                    return a * std::sqrt(3. / 2.) * lonsin;
                }
                else {
                    return a * std::sqrt(3. / 2.) * loncos;
                }
            }
        }
        else {
            return 0.;
        }
    }

    // divergence:
    if (ivar_in == 1) {
        if (ivar_out == 0) {  // u:
            if (m == 0 && n == 0) {
                return 0.;
            }
            if (m == 0 && n == 1) {
                return 0.;
            }
            if (m == 1 && n == 1) {
                if (imag == 0) {
                    return a * std::sqrt(3. / 2.) * lonsin;
                }
                else {
                    return a * std::sqrt(3. / 2.) * loncos;
                }
            }
        }
        else if (ivar_out == 1) {  // v:
            if (m == 0 && n == 0) {
                return 0.;
            }
            if (m == 0 && n == 1) {
                if (imag == 0) {
                    return -std::sqrt(3.) * a / 2. * latcos;
                }
                else {
                    return 0.;
                }
            }
            if (m == 1 && n == 1) {
                if (imag == 0) {
                    return a * std::sqrt(3. / 2.) * loncos * latsin;
                }
                else {
                    return -a * std::sqrt(3. / 2.) * lonsin * latsin;  // sign?
                }
            }
        }
        else {
            return 0.;
        }
    }

    return -1.;
}

extern "C" {
double ectrans_init_spherical_harmonic_hardcoded(int n, int m, double lon, double lat, bool imag) {
    const int abs_m = std::abs(m);

    if (n < abs_m) {
	throw std::runtime_error("Error in ectrans_init_spherical_harmonic: assertion (n >= abs_m) failed");
    }

    return sphericalharmonics_analytic_point(n, m, imag, lon, lat, 2, 2);
}
}

extern "C" {
double ectrans_init_spherical_harmonic_northsouth_derivative_hardcoded(int n, int m, double lon, double lat, bool imag) {

    double latsin = std::sin(lat), latcos = std::cos(lat);
    double lonsin = std::sin(m * lon), loncos = std::cos(m * lon);

    double rcosinv = 1. / (earth_radius);
    double rft = 1.;
    if (m > 0) {
        rft *= 2.;  // the famous factor 2 that noone really understands
    }
    if (imag == 0) {
        rft *= lonsin; // DEBUGGING: should be loncos!!!
    }
    else {
        rft *= -lonsin;
    }

    if (m == 0 && n == 1) {
        return std::sqrt(3.) * latcos * rcosinv;
    } else if (m == 1 && n == 1) {
        return -std::sqrt(3. / 2.) * rft * rcosinv;
    }
}
}

}  // namespace init 
}  // namespace ectrans
