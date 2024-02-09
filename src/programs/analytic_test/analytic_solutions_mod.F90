module analytic_solutions_mod

  use parkind1, only: jpim, jprd

  real(kind=jprd) :: z_pi = 4.d0*datan(1.d0)
  real(kind=jprd), allocatable :: zmu(:), legpolys(:,:,:), legpolys_ectrans(:,:,:), gelam(:,:), gelat(:,:)
  integer(kind=jpim), allocatable :: nlatidxs(:,:), nmeng(:)
  integer(kind=jpim) :: nfirstlat, nlastlat

#include "trans_inq.h"

  contains

  !===================================================================================================
  ! Subroutine analytic_init:
  ! Compute with the help of trans_inq the geographic longitude gelam and latitude gelat.
  ! Also create a helper array nlatidxs(nproma,ngpblks) which contains for each blocked point the
  ! global latitude index. This is used later to retrieve the corresponding Legendre polynomial.
  !===================================================================================================

  subroutine analytic_init(nproma, ngpblks, ndgl, n_regions_ns, n_regions_ew, nloen)

    implicit none
  
    integer(kind=jpim), intent(in) :: nproma, ngpblks, ndgl, n_regions_ns, n_regions_ew
    integer(kind=jpim), dimension(ndgl), intent(in) :: nloen
    integer(kind=jpim) :: nptrfloff, my_region_ns, my_region_ew
    integer(kind=jpim) :: jglat, ioff, ilat, istlon, iendlon, jlon, jrof, ibl
    integer(kind=jpim), dimension(n_regions_ns) :: nfrstlat, nlstlat
    integer(kind=jpim), dimension(ndgl+n_regions_ns-1,n_regions_ew) :: nsta, nonl
    real(kind=jprd) :: zlat, zlon
  
    allocate(zmu(ndgl),nlatidxs(nproma,ngpblks),gelam(nproma,ngpblks),gelat(nproma,ngpblks),nmeng(ndgl))

    call trans_inq(kptrfloff=nptrfloff, &
       & kmy_region_ns=my_region_ns,    &
       & kmy_region_ew=my_region_ew,    &
       & kfrstlat=nfrstlat,             &
       & klstlat=nlstlat,               &
       & ksta=nsta,                     &
       & konl=nonl,                     &
       & pmu=zmu,                       &
       & knmeng=nmeng)
  
    ilat = nptrfloff
    ibl  = 1
    jrof = 1
    nfirstlat = nfrstlat(my_region_ns)
    nlastlat = nlstlat(my_region_ns)
    do jglat = nfirstlat, nlastlat
      zlat = asin(zmu(jglat))
      ilat = ilat + 1
      istlon = nsta(ilat,my_region_ew)
      iendlon = istlon - 1 + nonl(ilat,my_region_ew)
      do jlon = istlon, iendlon
        zlon = real(jlon-1,jprd)*2.0_jprd*z_pi/real(nloen(jglat),jprd)
        gelam(jrof,ibl) = zlon
        gelat(jrof,ibl) = zlat
        nlatidxs(jrof,ibl) = jglat
        jrof = jrof + 1
        if(jrof > nproma) then
          jrof = 1
          ibl  = ibl + 1
        end if
      end do
    end do
  
  end subroutine analytic_init

  !===================================================================================================  
  ! Subroutine analytic_end:
  ! Deallocate the helper arrays used for the analytic solutions.
  !===================================================================================================

  subroutine analytic_end()

    implicit none

    deallocate(zmu, nlatidxs, gelam, gelat)

    if(allocated(legpolys)) deallocate(legpolys)

  end subroutine analytic_end
  
  !===================================================================================================
  ! Subroutine buffer_legendre_polynomials_supolf:
  ! Compute Legendre polynomials with the help of supolf_test_mod and store them for all latitudes
  ! and zonal wavenumbers.
  !===================================================================================================

  subroutine buffer_legendre_polynomials_supolf(nsmax)
  
    use supolf_test_mod, only: supolf_test
    use tpm_pol_test, only: ini_pol_test, dfa
  
    implicit none
  
    integer(kind=jpim), intent(in) :: nsmax
    integer(kind=jpim) :: km, jgl
    integer(kind=jpim), dimension(nsmax) :: ndglu

    if(allocated(legpolys)) deallocate(legpolys)
    allocate(legpolys(nfirstlat:nlastlat, 0:nsmax+1, 0:nsmax+1))
    legpolys = 0.0

    call ini_pol_test(nsmax+1)
    do jgl=nfirstlat, nlastlat
      do km=0,nmeng(jgl)
        call supolf_test(km, nmeng(jgl)+1, zmu(jgl), legpolys(jgl,km,0:nmeng(jgl)+1))
        ! going only to nmeng(jgl) with km and inside the call makes sure that we skip the highest
        ! wavenumbers at high latitudes for the reduced grid
      end do
    end do
  
  end subroutine buffer_legendre_polynomials_supolf
  
  !===================================================================================================
  ! Subroutine buffer_legendre_polynomials_ectrans:
  ! Using ectrans via trans_inq to retrieve the Legendre polynomials. These are only used to check if
  ! ectrans is computing the same values like the copied supolf_test_mod. They should not be used to
  ! compute the analytic solutions since then we would not detect bugs introduced in the computation
  ! of the Legendre polynomials inside ectrans.
  ! Caution: in the single precision version ectrans only returns the Legendre polynomials in single
  !          precision! To compute the errors as accurately as possible the errors should always be
  !          calculated from the Legendre coefficients in double precision (which is the case as 
  !          long as supolf_test is used)
  !===================================================================================================

  subroutine buffer_legendre_polynomials_ectrans(nsmax, ndgl)
  
    use supolf_test_mod, only: supolf_test
    use tpm_pol_test, only: ini_pol_test, dfa
    use parkind1, only: jprb
  
    implicit none
  
    integer(kind=jpim), intent(in) :: nsmax, ndgl
    real(kind=jprb), dimension(ndgl,(nsmax+2)*(nsmax+3)/2) :: rpnm
    integer(kind=jpim) :: jnm, jn, jninv, jm, ilat
  
    if(allocated(legpolys_ectrans)) deallocate(legpolys_ectrans)
    allocate(legpolys_ectrans(nfirstlat:nlastlat, 0:nsmax, 0:nsmax))
    legpolys_ectrans = 0.0

    call trans_inq(prpnm=rpnm)
  
    do ilat=nfirstlat, nlastlat
      jnm = 0
      do jm=0,nsmax+1
        do jninv=jm,nsmax+1
          jn = nsmax+1-jninv+jm
          jnm = jnm + 1
          if(jm<=nsmax .and. jn<=nsmax) legpolys_ectrans(ilat,jm,jn) = rpnm(ilat-nfirstlat+1,jnm)
        end do
      end do
    end do
  
  end subroutine buffer_legendre_polynomials_ectrans
  
  !===================================================================================================
  ! Subroutine compute_analytic_solution:
  ! Compute analytic solution for a specific total wavenumber n and zonal wavenumber m by going through
  ! all points and using the point-wise function analytic_spherical_harmonic_point.
  ! For kimag==.true. the imaginary part is computed. For kimag==.false. the real part is computed.
  !===================================================================================================
  
  subroutine compute_analytic_solution(nproma, ngpblks, nsmax, ngptot, kzonal, ktotal, kimag, sph_analytic)
  
    implicit none
  
    integer(kind=jpim), intent(in) :: nproma, ngpblks, nsmax, ngptot
    real(kind=jprd), dimension(nproma,ngpblks), intent(out) :: sph_analytic
    integer(kind=jpim), intent(in) :: kzonal, ktotal
    logical, intent(in) :: kimag
    integer(kind=jpim) :: jkglo, iend, ioff, ibl, jrof
  
    do jkglo=1,ngptot,nproma
      iend = min(nproma,ngptot-jkglo+1)
      ioff = jkglo - 1
      ibl  = (jkglo-1)/nproma+1
      do jrof=1,iend
        sph_analytic(jrof,ibl) = analytic_spherical_harmonic_point( ktotal, kzonal, gelam(jrof,ibl), gelat(jrof,ibl), kimag, legpolys(nlatidxs(jrof,ibl), kzonal, ktotal))
      end do
    end do
  
  end subroutine compute_analytic_solution
  
  !===================================================================================================
  ! Function analytic_spherical_harmonic_point:
  ! Compute analytic solution for a single point with longitude lon and latitude lat and for a specific
  ! total wavenumber n and zonal wavenumber m by using the Legendre coefficient legPoly.
  ! For imag==.true. the imaginary part is computed. For imag==.false. the real part is computed.
  ! For an example of how to use this function see the subroutine compute_analytic_solution.
  !===================================================================================================
  
  function analytic_spherical_harmonic_point(n, m, lon, lat, imag, legPoly) result(pointValue)

    implicit none

    integer(jpim), intent(in), value :: n
    integer(jpim), intent(in), value :: m
    real(jprd), intent(in), value    :: lon
    real(jprd), intent(in), value    :: lat
    logical, intent(in), value       :: imag
    real(jprd), intent(in), value    :: legPoly
    real(jprd) :: pointValue
    integer(jpim) :: abs_m
    abs_m = abs(m)
    if(n<abs_m) call abor1("Error in analytic_spherical_harmonic_point: assertion (n >= abs_m) failed")
    if (m == 0) then
      if (imag) then
        pointValue = 0.0
      else
        pointValue = legPoly
      end if
    end if
  
    if (m > 0) then
      if (imag) then
          pointValue = (-2 * sin(m * lon) * legPoly)
      else
          pointValue = (2 * cos(m * lon) * legPoly)
      end if
    end if
  end function analytic_spherical_harmonic_point
  
  !===================================================================================================
  ! Function analytic_eastwest_derivative_point:
  ! Compute analytic solution of the east-west derivative for a single point with longitude lon and 
  ! latitude lat and for a specific total wavenumber n and zonal wavenumber m by using the Legendre
  ! coefficient legPoly. For imag==.true. the imaginary part is computed. For imag==.false. the real
  ! part is computed. For an example of how to use this function see the subroutine
  ! compute_analytic_eastwest_derivative.
  !===================================================================================================
  
  function analytic_eastwest_derivative_point(n, m, lon, lat, imag, legPoly) result(pointValue)

    use tpm_constants, only: ra

    implicit none

    integer(jpim), intent(in), value :: n
    integer(jpim), intent(in), value :: m
    real(jprd), intent(in), value    :: lon
    real(jprd), intent(in), value    :: lat
    logical, intent(in), value       :: imag
    real(jprd), intent(in), value    :: legPoly
    real(jprd) :: pointValue
    if (imag) then
      pointValue = (- m * analytic_spherical_harmonic_point(n, m, lon, lat, .false., legPoly) / (ra * cos(lat)));
    else
      pointValue = (m * analytic_spherical_harmonic_point(n, m, lon, lat, .true., legPoly) / (ra * cos(lat)));
    end if
  
  end function analytic_eastwest_derivative_point
  
  !===================================================================================================
  ! Function analytic_northsouth_derivative_point:
  ! Compute analytic solution of the north-south derivative for a single point with longitude lon and 
  ! latitude lat and for a specific total wavenumber n and zonal wavenumber m by using the Legendre
  ! coefficients legpolyp1 (for total wavenumber n+1) and legpolym1 (for total wavenumber n-1).
  ! For imag==.true. the imaginary part is computed. For imag==.false. the real part is computed.
  ! For an example of how to use this function see the subroutine compute_analytic_northsouth_derivative.
  !===================================================================================================
  
  function analytic_northsouth_derivative_point(n, m, lon, lat, imag, legpolyp1, legpolym1) result(pointValue)

    use tpm_constants, only: ra

    implicit none

    integer(jpim), intent(in), value :: n
    integer(jpim), intent(in), value :: m
    real(jprd), intent(in), value    :: lon
    real(jprd), intent(in), value    :: lat
    logical, intent(in), value       :: imag
    real(jprd), intent(in), value    :: legpolyp1, legpolym1
    real(jprd) :: pointValue
    real(jprd) :: sinlat, coslat, colat, coeff_a, coeff_b
    integer(jpim) :: abs_m
    sinlat = sin(lat)
    coslat = cos(lat)
    abs_m = abs(m)
    if(n<abs_m) call abor1("Error in analytic_spherical_harmonic_point_northsouth_derivative: assertion (n >= abs_m) failed")
    colat = z_pi/2.0 - lat
    coeff_a = (n+1)*sqrt(real(n*n-m*m,jprd)/real(4.0*n*n-1,jprd))
    coeff_b = -n*sqrt(real((n+1)*(n+1)-m*m,jprd)/real(4.0*(n+1)*(n+1)-1,jprd))
    if (m == 0) then
      if (imag) then
        pointValue = 0.0
      else
        if (n > m) then
          pointValue = (coeff_a * legpolym1 + coeff_b * legpolyp1) / (ra * coslat)
        else
          pointValue = (coeff_b * legpolyp1) / (ra * coslat)
        end if
      end if
    end if
  
    if (m > 0) then
      if (n > m) then
        if (imag) then
            pointValue = -2 * sin(m * lon) * (coeff_a * legpolym1 + coeff_b * legpolyp1) / (ra * coslat);
        else
            pointValue = 2 * cos(m * lon) * (coeff_a * legpolym1 + coeff_b * legpolyp1) / (ra * coslat);
        end if
      else
        if (imag) then
            pointValue = -2 * sin(m * lon) * (coeff_b * legpolyp1) / (ra * coslat);
        else
            pointValue = 2 * cos(m * lon) * (coeff_b * legpolyp1) / (ra * coslat);
        end if
      end if
    end if
  
  end function analytic_northsouth_derivative_point
  
  !===================================================================================================
  ! Subroutine compute_analytic_eastwest_derivative:
  ! Compute analytic solution of the east-west derivative for a specific total wavenumber n and zonal
  ! wavenumber m by going through all points and using the point-wise function
  ! analytic_eastwest_derivative_point. For kimag==.true. the imaginary part is computed.
  ! For kimag==.false. the real part is computed.
  !===================================================================================================
  
  subroutine compute_analytic_eastwest_derivative(nproma, ngpblks, nsmax, ngptot, kzonal, ktotal, kimag, sph_analytic)
  
    implicit none
  
    integer(kind=jpim), intent(in) :: nproma, ngpblks, nsmax, ngptot
    real(kind=jprd), dimension(nproma,ngpblks), intent(out) :: sph_analytic
    integer(kind=jpim), intent(in) :: kzonal, ktotal
    logical, intent(in) :: kimag
    integer(kind=jpim) :: jkglo, iend, ioff, ibl, jrof
  
    do jkglo=1,ngptot,nproma
      iend = min(nproma,ngptot-jkglo+1)
      ioff = jkglo - 1
      ibl  = (jkglo-1)/nproma+1
      do jrof=1,iend
        sph_analytic(jrof,ibl) = analytic_eastwest_derivative_point( ktotal, kzonal, gelam(jrof,ibl), gelat(jrof,ibl), kimag, legpolys(nlatidxs(jrof,ibl), kzonal, ktotal))
      end do
    end do
  
  end subroutine compute_analytic_eastwest_derivative
  
  !===================================================================================================
  ! Subroutine compute_analytic_northsouth_derivative:
  ! Compute analytic solution of the north-south derivative for a specific total wavenumber n and zonal
  ! wavenumber m by going through all points and using the point-wise function
  ! analytic_northsouth_derivative_point. For kimag==.true. the imaginary part is computed.
  ! For kimag==.false. the real part is computed.
  !===================================================================================================
  
  subroutine compute_analytic_northsouth_derivative(nproma, ngpblks, nsmax, ngptot, kzonal, ktotal, kimag, sph_analytic)
  
    implicit none
  
    integer(kind=jpim), intent(in) :: nproma, ngpblks, nsmax, ngptot
    real(kind=jprd), dimension(nproma,ngpblks), intent(out) :: sph_analytic
    integer(kind=jpim), intent(in) :: kzonal, ktotal
    logical, intent(in) :: kimag
    real(kind=jprd) :: legpolyp1, legpolym1
    integer(kind=jpim) :: jkglo, iend, ioff, ibl, jrof
  
    legpolyp1 = 0.0
    legpolym1 = 0.0
  
    do jkglo=1,ngptot,nproma
      iend = min(nproma,ngptot-jkglo+1)
      ioff = jkglo - 1
      ibl  = (jkglo-1)/nproma+1
      do jrof=1,iend
        if(kzonal<=nmeng(nlatidxs(jrof,ibl))) then
          if(ktotal<nsmax+1) then
            legpolyp1 = legpolys(nlatidxs(jrof,ibl), kzonal, ktotal+1)
          end if
          if(ktotal>0) then
            legpolym1 = legpolys(nlatidxs(jrof,ibl), kzonal, ktotal-1)
          end if
          sph_analytic(jrof,ibl) = analytic_northsouth_derivative_point( ktotal, kzonal, gelam(jrof,ibl), gelat(jrof,ibl), kimag, legpolyp1, legpolym1)
        else
          sph_analytic(jrof,ibl) = 0.0
        end if
      end do
    end do
  
  end subroutine compute_analytic_northsouth_derivative
  
  !===================================================================================================
  ! Subroutine compute_analytic_uv:
  ! Compute analytic solution of the wind speed components u and v for a specific total wavenumber n
  ! and zonal wavenumber m by going through all points and using the point-wise function
  ! analytic_spherical_harmonic_point. For kimag==.true. the imaginary part is computed.
  ! For kimag==.false. the real part is computed. For the maths behind this computation see the paper
  ! Temperton, 1991, MWR 119 p1303.
  !===================================================================================================
  
  subroutine compute_analytic_uv(nproma, ngpblks, nsmax, ngptot, m, n, kimag, u_analytic, v_analytic)
  
    use tpm_constants, only: ra

    implicit none
  
    integer(kind=jpim), intent(in) :: nproma, ngpblks, nsmax, ngptot
    real(kind=jprd), dimension(nproma,ngpblks), intent(out) :: u_analytic, v_analytic
    integer(kind=jpim), intent(in) :: m, n
    logical, intent(in) :: kimag
    real(kind=jprd) :: coeff0, coeff1, coeff2, coeff3, r1, r2, r3
    integer(kind=jpim) :: jkglo, iend, ioff, ibl, jrof
  
    if (kimag) then
      coeff0 = 1.0_jprd
    else
      coeff0 = -1.0_jprd
    end if
    if (n>0) then
      coeff1 = - real(ra,jprd)*real(ra,jprd)/real(n*(n+1),jprd)
      coeff2 = sqrt(real((n+1)*(n+1)-m*m,jprd)/real(4.0*(n+1)*(n+1)-1,jprd))
      coeff3 = sqrt(real(n*n-m*m,jprd)/real(4.0*n*n-1,jprd))
    else
      coeff1 = 0.0
      coeff2 = 0.0
      coeff3 = 0.0
    end if

    do jkglo=1,ngptot,nproma
      iend = min(nproma,ngptot-jkglo+1)
      ioff = jkglo - 1
      ibl  = (jkglo-1)/nproma+1
      do jrof=1,iend
        if(m<=nmeng(nlatidxs(jrof,ibl)) .and. n>0) then
          ! first terms in Temperton eq. (2.12) and (2.13):
          r1 = - real(m,jprd) * coeff0 * coeff1 * analytic_spherical_harmonic_point( n, m, gelam(jrof,ibl), gelat(jrof,ibl), .not. kimag, legpolys(nlatidxs(jrof,ibl), m, n)) / (ra * cos(gelat(jrof,ibl)))
          u_analytic(jrof,ibl) = r1
          v_analytic(jrof,ibl) = r1
          ! second terms psi_n-1^m in Temperton eq.(2.12) and (2.13) contribute for n+1:
          if(n<nsmax+1) then
            r2 = n * coeff1 * coeff2 * analytic_spherical_harmonic_point( n + 1, m, gelam(jrof,ibl), gelat(jrof,ibl), kimag, legpolys(nlatidxs(jrof,ibl), m, n + 1)) / (ra * cos(gelat(jrof,ibl)))
            u_analytic(jrof,ibl) = u_analytic(jrof,ibl) + r2
            v_analytic(jrof,ibl) = v_analytic(jrof,ibl) - r2
          end if
          ! third terms psi_n+1^m in Temperton eq.(2.12) and (2.13) contribute for n-1:
          if(n>m) then
            r3 = (n + 1) * coeff1 * coeff3 * analytic_spherical_harmonic_point( n - 1, m, gelam(jrof,ibl), gelat(jrof,ibl), kimag, legpolys(nlatidxs(jrof,ibl), m, n - 1)) / (ra * cos(gelat(jrof,ibl)))
            u_analytic(jrof,ibl) = u_analytic(jrof,ibl) - r3
            v_analytic(jrof,ibl) = v_analytic(jrof,ibl) + r3
          end if
        else
          u_analytic(jrof,ibl) = 0.0
          v_analytic(jrof,ibl) = 0.0
        end if
      end do
    end do
  
  end subroutine compute_analytic_uv
  
  !===================================================================================================
  ! Subroutine compute_analytic_uv_derivative_ew:
  ! Compute analytic solution of the east-west derivatives of the wind speed components u and v for a 
  ! specific total wavenumber n and zonal wavenumber m by going through all points and using the
  ! point-wise function analytic_eastwest_derivative_point. For kimag==.true. the imaginary part is computed.
  ! For kimag==.false. the real part is computed. For the maths behind this computation see the paper
  ! Temperton, 1991, MWR 119 p1303. Main difference compared to the subroutine compute_analytic_uv 
  ! is the use of analytic_eastwest_derivative_point instead of analytic_spherical_harmonic_point.
  !===================================================================================================
  
  subroutine compute_analytic_uv_derivative_ew(nproma, ngpblks, nsmax, ngptot, m, n, kimag, uder_analytic, vder_analytic)
  
    use tpm_constants, only: ra

    implicit none
  
    integer(kind=jpim), intent(in) :: nproma, ngpblks, nsmax, ngptot
    real(kind=jprd), dimension(nproma,ngpblks), intent(out) :: uder_analytic, vder_analytic
    integer(kind=jpim), intent(in) :: m, n
    logical, intent(in) :: kimag
    real(kind=jprd) :: coeff0, coeff1, coeff2, coeff3, r1, r2, r3
    integer(kind=jpim) :: jkglo, iend, ioff, ibl, jrof
  
    if (kimag) then
      coeff0 = 1.0_jprd
    else
      coeff0 = -1.0_jprd
    end if
    if (n>0) then
      coeff1 = - real(ra,jprd)*real(ra,jprd)/real(n*(n+1),jprd)
      coeff2 = sqrt(real((n+1)*(n+1)-m*m,jprd)/real(4.0*(n+1)*(n+1)-1,jprd))
      coeff3 = sqrt(real(n*n-m*m,jprd)/real(4.0*n*n-1,jprd))
    else
      coeff1 = 0.0
      coeff2 = 0.0
      coeff3 = 0.0
    end if

    do jkglo=1,ngptot,nproma
      iend = min(nproma,ngptot-jkglo+1)
      ioff = jkglo - 1
      ibl  = (jkglo-1)/nproma+1
      do jrof=1,iend
        if(m<=nmeng(nlatidxs(jrof,ibl)) .and. n>0) then
          ! first terms in Temperton eq. (2.12) and (2.13):
          r1 = - real(m,jprd) * coeff0 * coeff1 * analytic_eastwest_derivative_point( n, m, gelam(jrof,ibl), gelat(jrof,ibl), .not. kimag, legpolys(nlatidxs(jrof,ibl), m, n)) / (ra * cos(gelat(jrof,ibl)))
          uder_analytic(jrof,ibl) = r1
          vder_analytic(jrof,ibl) = r1
          ! second terms psi_n-1^m in Temperton eq.(2.12) and (2.13) contribute for n+1:
          if(n<nsmax+1) then
            r2 = n * coeff1 * coeff2 * analytic_eastwest_derivative_point( n + 1, m, gelam(jrof,ibl), gelat(jrof,ibl), kimag, legpolys(nlatidxs(jrof,ibl), m, n + 1)) / (ra * cos(gelat(jrof,ibl)))
            uder_analytic(jrof,ibl) = uder_analytic(jrof,ibl) + r2
            vder_analytic(jrof,ibl) = vder_analytic(jrof,ibl) - r2
          end if
          ! third terms psi_n+1^m in Temperton eq.(2.12) and (2.13) contribute for n-1:
          if(n>m) then
            r3 = (n + 1) * coeff1 * coeff3 * analytic_eastwest_derivative_point( n - 1, m, gelam(jrof,ibl), gelat(jrof,ibl), kimag, legpolys(nlatidxs(jrof,ibl), m, n - 1)) / (ra * cos(gelat(jrof,ibl)))
            uder_analytic(jrof,ibl) = uder_analytic(jrof,ibl) - r3
            vder_analytic(jrof,ibl) = vder_analytic(jrof,ibl) + r3
          end if
        else
          uder_analytic(jrof,ibl) = 0.0
          vder_analytic(jrof,ibl) = 0.0
        end if
      end do
    end do
  
  end subroutine compute_analytic_uv_derivative_ew
  
  !===================================================================================================
  
end module analytic_solutions_mod