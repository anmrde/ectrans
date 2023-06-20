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
  function ectrans_init_spherical_harmonic(n, m, lon, lat, imag) bind(c,name="ectrans_init_spherical_harmonic")
    use iso_c_binding, only: c_double, c_int, c_bool
    implicit none
    real(c_double)                     :: ectrans_init_spherical_harmonic
    integer(c_int), intent(in), value  :: n
    integer(c_int), intent(in), value  :: m
    real(c_double), intent(in), value  :: lon
    real(c_double), intent(in), value  :: lat
    logical, intent(in), value :: imag
  end function

  function ectrans_init_spherical_harmonic_eastwest_derivative(n, m, lon, lat, imag) bind(c,name="ectrans_init_spherical_harmonic_eastwest_derivative")
    use iso_c_binding, only: c_double, c_int, c_bool
    implicit none
    real(c_double)                     :: ectrans_init_spherical_harmonic_eastwest_derivative
    integer(c_int), intent(in), value  :: n
    integer(c_int), intent(in), value  :: m
    real(c_double), intent(in), value  :: lon
    real(c_double), intent(in), value  :: lat
    logical, intent(in), value :: imag
  end function

  function ectrans_init_spherical_harmonic_northsouth_derivative(n, m, lon, lat, imag) bind(c,name="ectrans_init_spherical_harmonic_northsouth_derivative")
    use iso_c_binding, only: c_double, c_int, c_bool
    implicit none
    real(c_double)                     :: ectrans_init_spherical_harmonic_northsouth_derivative
    integer(c_int), intent(in), value  :: n
    integer(c_int), intent(in), value  :: m
    real(c_double), intent(in), value  :: lon
    real(c_double), intent(in), value  :: lat
    logical, intent(in), value :: imag
  end function

  function ectrans_init_spherical_harmonic_northsouth_derivative_hardcoded(n, m, lon, lat, imag) bind(c,name="ectrans_init_spherical_harmonic_northsouth_derivative_hardcoded")
    use iso_c_binding, only: c_double, c_int, c_bool
    implicit none
    real(c_double)                     :: ectrans_init_spherical_harmonic_northsouth_derivative_hardcoded
    integer(c_int), intent(in), value  :: n
    integer(c_int), intent(in), value  :: m
    real(c_double), intent(in), value  :: lon
    real(c_double), intent(in), value  :: lat
    logical, intent(in), value :: imag
  end function

  function ectrans_init_spherical_harmonic_hardcoded(n, m, lon, lat, imag) bind(c,name="ectrans_init_spherical_harmonic_hardcoded")
    use iso_c_binding, only: c_double, c_int, c_bool
    implicit none
    real(c_double)                     :: ectrans_init_spherical_harmonic_hardcoded
    integer(c_int), intent(in), value  :: n
    integer(c_int), intent(in), value  :: m
    real(c_double), intent(in), value  :: lon
    real(c_double), intent(in), value  :: lat
    logical, intent(in), value :: imag
  end function
end interface

contains

end module
