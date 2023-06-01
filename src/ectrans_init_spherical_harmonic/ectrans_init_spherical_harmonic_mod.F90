module ectrans_init_spherical_harmonic_mod

public

interface

  ! \brief An analytic function that provides a spherical harmonic on a 2D sphere
  !
  ! \param n Total wave number
  ! \param m Zonal wave number
  ! \param lon Longitude in degrees
  ! \param lat Latitude in degrees
  ! \return spherical harmonic
  !
  !double ectrans_init_spherical_harmonic(int n, int m, double lon, double lat)
  function ectrans_init_spherical_harmonic(n, m, lon, lat) bind(c,name="ectrans_init_spherical_harmonic")
    use iso_c_binding, only: c_double, c_int
    implicit none
    real(c_double)                    :: ectrans_init_spherical_harmonic
    integer(c_int), intent(in), value :: n
    integer(c_int), intent(in), value :: m
    real(c_double), intent(in), value :: lon
    real(c_double), intent(in), value :: lat
  end function
end interface

contains

end module
