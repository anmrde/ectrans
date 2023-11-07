module analytic_solutions_mod

  use parkind1, only: jpim, jprd

  real(kind=jprd) :: z_pi = 4.d0*datan(1.d0)
  real(kind=jprd), allocatable :: zmu(:), legpolys(:,:,:), legpolys_ectrans(:,:,:), gelam(:,:), gelat(:,:)
  integer(kind=jpim), allocatable :: nlatidxs(:,:), nmeng(:)
  integer(kind=jpim) :: nfirstlat, nlastlat

  #include "trans_inq.h"

  contains

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
  
  subroutine analytic_end()

    implicit none

    deallocate(zmu, nlatidxs, gelam, gelat)

    if(allocated(legpolys)) deallocate(legpolys)

  end subroutine analytic_end
  
  !===================================================================================================
  
  subroutine buffer_legendre_polynomials(nsmax)
  
    implicit none
  
    integer(kind=jpim), intent(in) :: nsmax
    real(kind=jprd) :: x
    integer(kind=jpim) :: n, m, ilat
  
    if(allocated(legpolys)) deallocate(legpolys)
    allocate(legpolys(nfirstlat:nlastlat, 0:nsmax+1, 0:nsmax+1))
    legpolys = 0.0

    do n = 0,nsmax+1
      do ilat = nfirstlat, nlastlat
        x = zmu(ilat)
        legpolys(ilat, n, n) = (double_factorial(2 * n - 1) * sqrt(1. - x * x) ** n)
      end do
    end do
  
    do n = 1,nsmax+1
      do ilat = nfirstlat, nlastlat
        x = zmu(ilat)
        legpolys(ilat, n-1, n) = x * (2 * n - 1) * legpolys(ilat, n-1, n-1)
      end do
    end do
  
    do n = 2,nsmax+1
      do m = 0,n-2
        do ilat = nfirstlat, nlastlat
          x = zmu(ilat)
          legpolys(ilat, m, n) = (x * (2 * n - 1) * legpolys(ilat, m, n-1) - (n + m - 1) * legpolys(ilat, m, n-2)) / (n - m)
        end do
      end do
    end do
  
    do n = 0,nsmax+1
      do m = 0,n
        do ilat = nfirstlat, nlastlat
          legpolys(ilat, m, n) = K(n, m) * legpolys(ilat, m, n)
        end do
      end do
    end do
  
    call legendre_polynomials_adjust_reduced_grid(nsmax)

  end subroutine buffer_legendre_polynomials
  
  !===================================================================================================
  
  subroutine buffer_legendre_polynomials_belusov(nsmax)
  
    use supol_test_mod, only: supol_test
    use tpm_pol_test, only: ini_pol_test

    implicit none
  
    integer(kind=jpim), intent(in) :: nsmax
    real(kind=jprd), dimension(0:nsmax+1, 0:nsmax+1) :: zfn
    real(kind=jprd) :: zfnn
    integer(kind=jpim) :: jn, jgl, iodd
  
    if(allocated(legpolys)) deallocate(legpolys)
    allocate(legpolys(nfirstlat:nlastlat, 0:nsmax+1, 0:nsmax+1))
    legpolys = 0.0

    zfn(0,0) = 2._jprd
    do jn = 1,nsmax+1
      zfnn = zfn(0,0)
      do jgl = 1,jn
        zfnn = zfnn * sqrt(1._jprd - 0.25_jprd / real(jgl**2,jprd))
      end do
  
      iodd = mod(jn,2)
      zfn(jn,jn) = zfnn
      do jgl = 2,jn-iodd,2
        zfn(jn,jn-jgl) = zfn(jn,jn-jgl+2) * real((jgl-1)*(2*jn-jgl+2),jprd) / real(jgl*(2*jn-jgl+1),jprd)
      end do
    end do
  
    call ini_pol_test(nsmax+1)
    do jgl=nfirstlat, nlastlat
      call supol_test(nsmax+1, zmu(jgl), zfn, legpolys(jgl,:,:))
    end do
  
    call legendre_polynomials_adjust_reduced_grid(nsmax)
  
  end subroutine buffer_legendre_polynomials_belusov
  
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
      do km=0,nsmax+1
        call supolf_test(km, nsmax+1, zmu(jgl), legpolys(jgl,km,:))
      end do
    end do
  
    call legendre_polynomials_adjust_reduced_grid(nsmax)
  
  end subroutine buffer_legendre_polynomials_supolf
  
  !===================================================================================================
  
  ! Using ectrans via trans_inq to retrieve the Legendre polynomials
  ! Caution: ectrans only returns the Legendre polynomials in single precision!!! (supol and supolf return them in double precision)
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
  
  subroutine legendre_polynomials_adjust_reduced_grid(nsmax)

    implicit none

    integer(kind=jpim), intent(in) :: nsmax
    integer(kind=jpim) :: ilat, jm, jn

    do ilat=nfirstlat, nlastlat
      do jm=nmeng(ilat)+1,nsmax
        do jn=jm,nsmax
          legpolys(ilat,jm,jn) = 0.0
        end do
      end do
    end do

  end subroutine legendre_polynomials_adjust_reduced_grid

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
        !sph_analytic(jrof,ibl) = ectrans_init_spherical_harmonic( ktotal, kzonal, gelam(jrof,ibl), gelat(jrof,ibl), kimag)
        !print *,"C: result=",sph_analytic(jrof,ibl)
        sph_analytic(jrof,ibl) = analytic_spherical_harmonic_point( ktotal, kzonal, gelam(jrof,ibl), gelat(jrof,ibl), kimag, legpolys(nlatidxs(jrof,ibl), kzonal, ktotal))
        !print *,"Fortran: result=",sph_analytic(jrof,ibl)
        !stop "debugging"
        !sph_analytic(jrof,ibl) = ectrans_init_spherical_harmonic_hardcoded( ktotal, kzonal, gelam(jrof,ibl), gelat(jrof,ibl), kimag)
      end do
    end do
  
  end subroutine compute_analytic_solution
  
  !===================================================================================================
  
  function analytic_spherical_harmonic_point(n, m, lon, lat, imag, legpoly) result(pointValue)

    implicit none

    integer(jpim), intent(in), value :: n
    integer(jpim), intent(in), value :: m
    real(jprd), intent(in), value    :: lon
    real(jprd), intent(in), value    :: lat
    logical, intent(in), value       :: imag
    real(jprd), intent(in), value    :: legpoly
    real(jprd) :: pointValue
    real(jprd) :: sinlat, coslat, colat
    integer(jpim) :: abs_m
    sinlat = sin(lat)
    coslat = cos(lat)
    abs_m = abs(m)
    if(n<abs_m) call abor1("Error in analytic_spherical_harmonic_point: assertion (n >= abs_m) failed")
    colat = z_pi/2.0 - lat
    if (m == 0) then
      if (imag) then
        pointValue = 0.0
      else
        pointValue = legpoly
      end if
    end if
  
    if (m > 0) then
      if (imag) then
          pointValue = (-2 * sin(m * lon) * legpoly)
      else
          pointValue = (2 * cos(m * lon) * legpoly)
      end if
    end if
  !print*,"Fortran: K=",K(n,m)," P=",P(n, m, sinlat)," result=",pointValue
  end function analytic_spherical_harmonic_point
  
  function factorial(n) result(fact)

    implicit none

    integer(jpim) :: n, i
    real(jprd) :: fact
    if (n < 0) call abor1("factorial of negative number not defined!")
    !fact = int(gamma(n+1), kind=jpim)
    fact = 1.0
    do i = 2,n
      fact = fact * i
    end do
  end function factorial
  
  function double_factorial(x) result(fact)

    implicit none

    integer(jpim) :: x, y
    real(jprd) :: fact
    y = x
    if (y == 0 .or. y == -1) then
      fact = 1
    else
      fact = y
      do while (y > 2)
        y = y - 2
        fact = fact * y
      end do
    end if
  end function double_factorial
  
  recursive function P(n, m, x) result(legPoly)

    implicit none

    integer(jpim) :: n, m
    real(jprd) :: x, legPoly
    ! No recursive calculation needed
    if (n == m) then
        legPoly = (double_factorial(2 * m - 1) * sqrt(1. - x * x) ** m)
    end if
  
    if (n == m + 1) then
      legPoly = x * (2 * m + 1) * P(m, m, x)
    end if
  
    ! Formula 1
    if (n > m + 1) then
      legPoly =  (x * (2 * n - 1) * P(n - 1, m, x) - (n + m - 1) * P(n - 2, m, x)) / (n - m)
    end if
  end function P
  
  ! Pn: associated Legendre polynomials with normalization used in ectrans (see Belousov 1962 p.5, in particular eq.5)
  function Pn(n, m, x) result(legPolyNorm)

    implicit none

    integer(jpim) :: n, m
    real(jprd) :: x, legPolyNorm
  
    legPolyNorm = K(n, m) * P(n, m, x)
  
  end function Pn
  
  function K(n, m) result(kFactor)

    implicit none

    integer(jpim) :: n, m
    real(jprd) :: kFactor
    kFactor = sqrt(((2 * n + 1) * factorial(n - m)) / (factorial(n + m)))
  end function K
  
  function analytic_eastwest_derivative_point(n, m, lon, lat, imag, legpoly) result(pointValue)

    use tpm_constants, only: ra

    implicit none

    integer(jpim), intent(in), value :: n
    integer(jpim), intent(in), value :: m
    real(jprd), intent(in), value    :: lon
    real(jprd), intent(in), value    :: lat
    logical, intent(in), value       :: imag
    real(jprd), intent(in), value    :: legpoly
    real(jprd) :: pointValue
    if (imag) then
      pointValue = (- m * analytic_spherical_harmonic_point(n, m, lon, lat, .false., legpoly) / (ra * cos(lat)));
    else
      pointValue = (m * analytic_spherical_harmonic_point(n, m, lon, lat, .true., legpoly) / (ra * cos(lat)));
    end if
  
  end function analytic_eastwest_derivative_point
  
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
  
  subroutine check_legendre_polynomials(nsmax, ndgl)
  
    implicit none
  
    integer(kind=jpim), intent(in) :: nsmax, ndgl
    integer(kind=jpim) :: ilat, jm, jn
    logical :: lprint
  
  
    write(33320,'("Legendre coefficients")')
    write(33320,'("ilat   m   n ┃     supolf    ectrans ")')
    write(33320,'("━━━━━━━━━━━━━╋━━━━━━━━━━━━━━━━━━━━━━┿")')
    do ilat=nfirstlat, nlastlat
      if(ilat<ndgl/2) then
        do jm=0,nsmax
          do jn=jm,nsmax
            lprint = .false.
            !print*,ilat,jm,jn,abs(legpolys(ilat,jm,jn)),abs(legpolys_ectrans(ilat,jm,jn))
            if(legpolys(ilat,jm,jn)==legpolys(ilat,jm,jn) .and. legpolys_ectrans(ilat,jm,jn)==legpolys_ectrans(ilat,jm,jn)) then ! just to be sure that there are no nans
              if(max(abs(legpolys(ilat,jm,jn)),abs(legpolys_ectrans(ilat,jm,jn)))>0) then
                if(abs(legpolys(ilat,jm,jn)-legpolys_ectrans(ilat,jm,jn))/max(abs(legpolys(ilat,jm,jn)),abs(legpolys_ectrans(ilat,jm,jn)))>1E-5) then
                  lprint = .true.
                end if
              end if
            else
              !lprint = .true.
            end if
            if(lprint) write(33320,'(i4,i4,i4," ┃",e11.3,e11.3,e11.3,e11.3,e11.3)') ilat, jm, jn, legpolys(ilat,jm,jn), legpolys_ectrans(ilat,jm,jn),abs(legpolys(ilat,jm,jn)-legpolys_ectrans(ilat,jm,jn)),max(abs(legpolys(ilat,jm,jn)),abs(legpolys_ectrans(ilat,jm,jn))),abs(legpolys(ilat,jm,jn)-legpolys_ectrans(ilat,jm,jn))/max(abs(legpolys(ilat,jm,jn)),abs(legpolys_ectrans(ilat,jm,jn)))
          end do
        end do
      end if
    end do
    
  end subroutine check_legendre_polynomials

  !===================================================================================================

  subroutine init_check_fields(lwrite_errors, nsmax, myproc, cgrid)

    use parkind1, only: jprb, jprd

    implicit none

    logical, intent(in) :: lwrite_errors
    integer(kind=jpim), intent(in) :: nsmax
    integer(kind=jpim), intent(in) :: myproc
    character(len=16), intent(in)  :: cgrid
    character(len=100) :: filename
    character(len=2)  :: precision

    if (myproc == 1) then
      if (jprb == jprd) then
        precision = 'dp'
      else
        precision = 'sp'
      end if
      write(filename,'(3a,i0,3a)')'errors-',precision,'-gridpoint_T',nsmax,'_',trim(cgrid),'.txt'
      open(40000+nsmax, file = filename, status='replace')
      write(filename,'(3a,i0,3a)')'errors-',precision,'-spectral_T',nsmax,'_',trim(cgrid),'.txt'
      open(60000+nsmax, file = filename, status='replace')
      if(lwrite_errors) then
        write(40000+nsmax,'("lmax-error in grid point space")')
        write(40000+nsmax,'("                 ┃           grid point data        │      north-south derivative      │       east-west derivative       |       wind speed      |  east-west derivative ")')
        write(40000+nsmax,'("   m   n it. im. ┃        pgp       pgp2      pgp3a │        pgp       pgp2      pgp3a │        pgp       pgp2      pgp3a │          u          v │          u          v ")')
        write(40000+nsmax,'("━━━━━━━━━━━━━━━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┿━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┿━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┿━━━━━━━━━━━━━━━━━━━━━━━┿━━━━━━━━━━━━━━━━━━━━━━━")')
        write(60000+nsmax,'("lmax-error in spectral space")')
        write(60000+nsmax,'("                 ┃             lmax-error           │             l2 - error           │      location of max    ")')
        write(60000+nsmax,'("   m   n it. im. ┃     zspsc2    zspsc2b    zspsc3a │     zspsc2    zspsc2b    zspsc3a │ initial index_max  rank")')
        write(60000+nsmax,'("━━━━━━━━━━━━━━━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┿━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┿━━━━━━━━━━━━━━━━━━━━━━━━━━━")')
        call flush(40000+nsmax)
        call flush(60000+nsmax)
      end if
    end if
    
  end subroutine init_check_fields

  !===================================================================================================

  subroutine close_check_fields(lwrite_errors, nsmax, myproc)

    implicit none

    logical, intent(in) :: lwrite_errors
    integer(kind=jpim), intent(in) :: nsmax
    integer(kind=jpim), intent(in) :: myproc

    if (myproc == 1) then
      close(40000+nsmax)
      close(60000+nsmax)
    end if

  end subroutine close_check_fields

  !===================================================================================================

  function check_gp_fields(rtolerance, lwrite_errors, nflevg, nfld, jstep, nzonal, ntotal, &
    & limag, zreel, zgp2, zgp3a, zgpuv, zsph_analytic, znsde_analytic, zewde_analytic, &
    & zu_analytic, zv_analytic, zuder_analytic, zvder_analytic, nout, nsmax, luse_mpi, ngptotg, lscders, lvordiv, luvders, myproc) result(rlmax_error)

    use parkind1, only: jprb
    use mpl_module

    implicit none

    real(kind=jprd), intent(in) :: rtolerance ! jprb would be enough
    logical, intent(in) :: lwrite_errors
    integer(kind=jpim), intent(in) :: jstep, nzonal, ntotal, nflevg, nfld
    logical, intent(in) :: limag
    real(kind=jprd), intent(in) :: zreel(:,:,:), zgp2(:,:,:), zgp3a(:,:,:,:), zgpuv(:,:,:,:) ! should be jprb
    real(kind=jprd), intent(in) :: zsph_analytic(:,:), znsde_analytic(:,:), zewde_analytic(:,:), zu_analytic(:,:), zv_analytic(:,:), zuder_analytic(:,:), zvder_analytic(:,:)
    real(kind=jprd), intent(out) :: rlmax_error
    integer(kind=jpim), intent(in) :: nout, nsmax, ngptotg, myproc
    logical, intent(in) :: luse_mpi, lscders, lvordiv, luvders
    real(kind=jprd) :: rlmax_quo, rlmax_nsde_quo, rlmax_ewde_quo, rlmax_u_quo, rlmax_v_quo, rlmax_uder_quo, rlmax_vder_quo, rlmax_fac, rlmax_nsde_fac, rlmax_ewde_fac, rlmax_u_fac, rlmax_v_fac, rlmax_uder_fac, rlmax_vder_fac
    real(kind=jprd) :: rlmax_errors(nflevg*nfld+2), rlmax_errors_nsde(nflevg*nfld+2), rlmax_errors_ewde(nflevg*nfld+2), rlmax_errors_uv(2*nflevg), rlmax_errors_uvder(2*nflevg)
    real(kind=jprd) :: rl2_errors(nflevg*nfld+2), rl2_errors_nsde(nflevg*nfld+2), rl2_errors_ewde(nflevg*nfld+2), rl2_errors_uv(2*nflevg), rl2_errors_uvder(2*nflevg)
    logical :: lpassed(nflevg*nfld+2), lpassed_nsde(nflevg*nfld+2), lpassed_ewde(nflevg*nfld+2), lpassed_all
    integer :: i, j, ntests

    ntests = nflevg*nfld+2
    rlmax_quo = maxval(abs( zsph_analytic(:,:)))
    rlmax_fac = 1.0_jprd
    if (luse_mpi) then
      call mpl_allreduce(rlmax_quo, 'max', ldreprod=.false.)
    end if
    if(rlmax_quo>0.0) rlmax_fac = 1.0_jprd/ rlmax_quo
    rlmax_errors(1) = maxval(abs(zreel(:,1,:)- zsph_analytic(:,:)))*rlmax_fac
    rlmax_errors(2) = maxval(abs( zgp2(:,1,:)- zsph_analytic(:,:)))*rlmax_fac
    rl2_errors(1) = sum((zreel(:,1,:)- zsph_analytic(:,:))**2)
    rl2_errors(2) = sum(( zgp2(:,1,:)- zsph_analytic(:,:))**2)
    do j=1,nflevg
      do i=1,nfld
        rlmax_errors(i*j+2) = maxval(abs(zgp3a(:,j,i,:)- zsph_analytic(:,:)))*rlmax_fac
        rl2_errors(i*j+2)   = sum((zgp3a(:,j,i,:)- zsph_analytic(:,:))**2)
      end do
    end do
    if (luse_mpi) then
      call mpl_allreduce(rlmax_errors, 'max', ldreprod=.false.)
      call mpl_allreduce(rl2_errors,   'sum', ldreprod=.false.)
    end if
    rl2_errors = sqrt(rl2_errors/ngptotg)
    rlmax_error = maxval(rlmax_errors)
    if (lscders) then
      rlmax_nsde_quo = maxval(abs(znsde_analytic(:,:)))
      rlmax_ewde_quo = maxval(abs(zewde_analytic(:,:)))
      rlmax_nsde_fac = 1.0_jprd
      rlmax_ewde_fac = 1.0_jprd
      if (luse_mpi) then
        call mpl_allreduce(rlmax_nsde_quo, 'max', ldreprod=.false.)
        call mpl_allreduce(rlmax_ewde_quo, 'max', ldreprod=.false.)
      end if
      if(rlmax_nsde_quo>0.0) rlmax_nsde_fac = 1.0_jprd/rlmax_nsde_quo
      if(rlmax_ewde_quo>0.0) rlmax_ewde_fac = 1.0_jprd/rlmax_ewde_quo
      rlmax_errors_nsde(1) = maxval(abs(zreel(:,2,:)-znsde_analytic(:,:)))*rlmax_nsde_fac
      rlmax_errors_nsde(2) = maxval(abs( zgp2(:,2,:)-znsde_analytic(:,:)))*rlmax_nsde_fac
      rl2_errors_nsde(1) = sum((zreel(:,2,:)- znsde_analytic(:,:))**2)
      rl2_errors_nsde(2) = sum(( zgp2(:,2,:)- znsde_analytic(:,:))**2)
      do j=1,nflevg
        do i=1,nfld
          rlmax_errors_nsde(i*j+2) = maxval(abs(zgp3a(:,j,nfld+i,:)-znsde_analytic(:,:)))*rlmax_nsde_fac
          rl2_errors_nsde(i*j+2)   = sum((zgp3a(:,j,nfld+i,:)- znsde_analytic(:,:))**2)
        end do
      end do
      rlmax_errors_ewde(1) = maxval(abs(zreel(:,3,:)-zewde_analytic(:,:)))*rlmax_ewde_fac
      rlmax_errors_ewde(2) = maxval(abs( zgp2(:,3,:)-zewde_analytic(:,:)))*rlmax_ewde_fac
      rl2_errors_ewde(1) = sum((zreel(:,3,:)- zewde_analytic(:,:))**2)
      rl2_errors_ewde(2) = sum(( zgp2(:,3,:)- zewde_analytic(:,:))**2)
      do j=1,nflevg
        do i=1,nfld
          rlmax_errors_ewde(i*j+2) = maxval(abs(zgp3a(:,j,2*nfld+i,:)-zewde_analytic(:,:)))*rlmax_ewde_fac
          rl2_errors_ewde(i*j+2)   = sum((zgp3a(:,j,2*nfld+i,:)- zewde_analytic(:,:))**2)
        end do
      end do
      if (luse_mpi) then
        call mpl_allreduce(rlmax_errors_nsde, 'max', ldreprod=.false.)
        call mpl_allreduce(rlmax_errors_ewde, 'max', ldreprod=.false.)
        call mpl_allreduce(rl2_errors_nsde,   'sum', ldreprod=.false.)
        call mpl_allreduce(rl2_errors_ewde,   'sum', ldreprod=.false.)
      end if
      rl2_errors_nsde = sqrt(rl2_errors_nsde/ngptotg)
      rl2_errors_ewde = sqrt(rl2_errors_ewde/ngptotg)
      rlmax_error = max(rlmax_error, maxval(rlmax_errors_nsde), maxval(rlmax_errors_ewde))
    else
      rlmax_errors_nsde = 0.0
      rlmax_errors_ewde = 0.0
    end if
    if (lvordiv) then
      rlmax_u_quo = maxval(abs(zu_analytic(:,:)))
      rlmax_v_quo = maxval(abs(zv_analytic(:,:)))
      rlmax_u_fac = 1.0_jprd
      rlmax_v_fac = 1.0_jprd
      if (luse_mpi) then
        call mpl_allreduce(rlmax_u_quo, 'max', ldreprod=.false.)
        call mpl_allreduce(rlmax_v_quo, 'max', ldreprod=.false.)
      end if
      if(rlmax_u_quo>0.0) rlmax_u_fac = 1.0_jprd/rlmax_u_quo
      if(rlmax_v_quo>0.0) rlmax_v_fac = 1.0_jprd/rlmax_v_quo
      do j=1,nflevg
        rlmax_errors_uv(2*j-1) = maxval(abs(zgpuv(:,j,1,:) - zu_analytic(:,:)))*rlmax_u_fac
        rlmax_errors_uv(2*j)   = maxval(abs(zgpuv(:,j,2,:) - zv_analytic(:,:)))*rlmax_v_fac
        rl2_errors_uv(2*j-1) = sum((zgpuv(:,j,1,:) - zu_analytic(:,:))**2)
        rl2_errors_uv(2*j)   = sum((zgpuv(:,j,2,:) - zv_analytic(:,:))**2)
      end do
      if (luse_mpi) then
        call mpl_allreduce(rlmax_errors_uv, 'max', ldreprod=.false.)
        call mpl_allreduce(rl2_errors_uv,   'sum', ldreprod=.false.)
      end if
      rl2_errors_uv = sqrt(rl2_errors_uv/ngptotg)
      rlmax_error = max(rlmax_error, maxval(rlmax_errors_uv))
    else
      rlmax_errors_uv = 0.0
    end if
    if (luvders) then
      rlmax_uder_quo = maxval(abs(zuder_analytic(:,:)))
      rlmax_vder_quo = maxval(abs(zvder_analytic(:,:)))
      rlmax_uder_fac = 1.0_jprd
      rlmax_vder_fac = 1.0_jprd
      if (luse_mpi) then
        call mpl_allreduce(rlmax_uder_quo, 'max', ldreprod=.false.)
        call mpl_allreduce(rlmax_vder_quo, 'max', ldreprod=.false.)
      end if
      if(rlmax_uder_quo>0.0) rlmax_uder_fac = 1.0_jprd/rlmax_uder_quo
      if(rlmax_vder_quo>0.0) rlmax_vder_fac = 1.0_jprd/rlmax_vder_quo
      do j=1,nflevg
        rlmax_errors_uvder(2*j-1) = maxval(abs(zgpuv(:,j,3,:) - zuder_analytic(:,:)))*rlmax_uder_fac
        rlmax_errors_uvder(2*j)   = maxval(abs(zgpuv(:,j,4,:) - zvder_analytic(:,:)))*rlmax_vder_fac
        rl2_errors_uvder(2*j-1) = sum((zgpuv(:,j,3,:) - zuder_analytic(:,:))**2)
        rl2_errors_uvder(2*j)   = sum((zgpuv(:,j,4,:) - zvder_analytic(:,:))**2)
      end do
      if (luse_mpi) then
        call mpl_allreduce(rlmax_errors_uvder, 'max', ldreprod=.false.)
        call mpl_allreduce(rl2_errors_uvder,   'sum', ldreprod=.false.)
      end if
      rl2_errors_uvder = sqrt(rl2_errors_uvder/ngptotg)
      rlmax_error = max(rlmax_error, maxval(rlmax_errors_uvder))
    else
      rlmax_errors_uvder = 0.0
    end if
    if(lwrite_errors .and. myproc == 1) then
      write(40000+nsmax,'(3i4,L4" ┃",3e11.3," │",3e11.3," │",3e11.3," │",2e11.3," │",2e11.3)') nzonal,ntotal,jstep,limag, &
        & rlmax_errors(1), &
        & rlmax_errors(2), &
        & rlmax_errors(3), &
        & rlmax_errors_nsde(1), &
        & rlmax_errors_nsde(2), &
        & rlmax_errors_nsde(3), &
        & rlmax_errors_ewde(1), &
        & rlmax_errors_ewde(2), &
        & rlmax_errors_ewde(3), &
        & rlmax_errors_uv(1), &
        & rlmax_errors_uv(2), &
        & rlmax_errors_uvder(1), &
        & rlmax_errors_uvder(2)
    !   & sqrt(sum((zgp3a(:,1,3,:)-zewde_analytic(:,:))**2)/ngptot)*lmaxewde_fac
      call flush(40000+nsmax)
    end if

    do i=1,ntests
      lpassed(i) = (rlmax_errors(i) < rtolerance)
      lpassed_nsde(i) = (rlmax_errors_nsde(i) < rtolerance)
      lpassed_ewde(i) = (rlmax_errors_ewde(i) < rtolerance)
    end do

    lpassed_all = rlmax_error < rtolerance

    if(.not. lpassed_all) then
      if(lwrite_errors) then
        write(50000+myproc,'("myproc=",i4," m=",i4," n=",i4," jstep=",i4," imag=",L1)')myproc,nzonal,ntotal,jstep,limag
        write(50000+myproc,'(a,a)')"                     ┃                                                        │                 north-south derivatives                │", &
        & "                 east-west derivatives                  │       wind speed component u     │       wind speed component v     │    east-west derivative of u     │    east-west derivative of v     "
        write(50000+myproc,'(a,a)')"     lat/°     lon/° ┃ analytical        pgp       pgp2      pgp3a  max-error │ analytical        pgp       pgp2      pgp3a  max-error │", &
        & " analytical        pgp       pgp2      pgp3a  max-error │ analytical      zgpuv  max-error │ analytical      zgpuv  max-error │ analytical      zgpuv  max-error │ analytical      zgpuv  max-error "
        write(50000+myproc,'(a,a)')"━━━━━━━━━━━━━━━━━━━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┿━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┿", &
        & "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┿━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┿━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┿━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┿━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
        do j=1,ubound(gelat,2)
          do i=1,ubound(gelat,1)
            if(max(abs(zreel(i,1,j)-zsph_analytic(i,j))*rlmax_fac,abs(zgp2(i,1,j)-zsph_analytic(i,j) )*rlmax_fac, &
              & maxval(abs((zgp3a(i,:,1,j)-zsph_analytic(i,j) )))*rlmax_fac, abs(zreel(i,2,j)-znsde_analytic(i,j))*rlmax_nsde_fac, &
              & abs(zgp2(i,2,j)-znsde_analytic(i,j))*rlmax_nsde_fac, maxval(abs((zgp3a(i,:,2,j)-znsde_analytic(i,j))))*rlmax_nsde_fac, &
              & abs(zreel(i,3,j)-zewde_analytic(i,j))*rlmax_ewde_fac, abs(zgp2(i,3,j)-zewde_analytic(i,j))*rlmax_ewde_fac, &
              & maxval(abs((zgp3a(i,:,3,j)-zewde_analytic(i,j))))*rlmax_ewde_fac, maxval(abs(zgpuv(i,:,1,j)-zu_analytic(i,j)))*rlmax_u_fac, &
              & maxval(abs(zgpuv(i,:,2,j)-zv_analytic(i,j)))*rlmax_v_fac) > rtolerance) then
              write(50000+myproc,'(2f10.3," ┃",5e11.3," │",5e11.3," │",5e11.3," │",3e11.3," │",3e11.3," │",3e11.3," │",3e11.3)') &
                & gelat(i,  j)*180/z_pi,gelam(i,j)*180/z_pi, &
                & abs(zsph_analytic(i,j)), abs(zreel(i,1,j)), abs(zgp2(i,1,j)), maxval(abs(zgp3a(i,:,1,j))), &
                & max(         (zreel(i,  1,j)-zsph_analytic(i,j) )*rlmax_fac, &
                &               (zgp2(i,  1,j)-zsph_analytic(i,j) )*rlmax_fac, &
                &   maxval(abs((zgp3a(i,:,1,j)-zsph_analytic(i,j) )))*rlmax_fac)     , &
                & znsde_analytic(i,j), zreel(i,2,j), zgp2(i,2,j), maxval(abs(zgp3a(i,:,2,j))), &
                & max(           (zreel(i,2,j)-znsde_analytic(i,j))*rlmax_nsde_fac, &
                &                 (zgp2(i,2,j)-znsde_analytic(i,j))*rlmax_nsde_fac, &
                &   maxval(abs((zgp3a(i,:,2,j)-znsde_analytic(i,j))))*rlmax_nsde_fac), &
                & zewde_analytic(i,j), zreel(i,3,j), zgp2(i,3,j), maxval(abs(zgp3a(i,:,3,j))), &
                & max(           (zreel(i,3,j)-zewde_analytic(i,j))*rlmax_ewde_fac, &
                &                 (zgp2(i,3,j)-zewde_analytic(i,j))*rlmax_ewde_fac, &
                &   maxval(abs((zgp3a(i,:,3,j)-zewde_analytic(i,j))))*rlmax_ewde_fac), &
                & zu_analytic(i,j), zgpuv(i,1,1,j), maxval(abs(zgpuv(i,:,1,j)-zu_analytic(i,j))*rlmax_u_fac), &
                & zv_analytic(i,j), zgpuv(i,1,2,j), maxval(abs(zgpuv(i,:,2,j)-zv_analytic(i,j))*rlmax_v_fac), & !maxval(abs(zgpuv(i,:,2,j)))
                & zuder_analytic(i,j), zgpuv(i,1,3,j), maxval(abs(zgpuv(i,:,3,j)-zuder_analytic(i,j))*rlmax_uder_fac), &
                & zvder_analytic(i,j), zgpuv(i,1,4,j), maxval(abs(zgpuv(i,:,4,j)-zvder_analytic(i,j))*rlmax_vder_fac) !maxval(abs(zgpuv(i,:,2,j)))
            end if
          end do
        end do
      end if
      if(sum(abs(int(not(lpassed))))+sum(abs(int(not(lpassed_nsde))))+sum(abs(int(not(lpassed_ewde)))) == 3*ntests) then
        write(nout,'("Error: all tests (fields and derivatives) fail for m=",i4," n=",i4)')nzonal,ntotal
      else if((sum(abs(int(not(lpassed)))) == ntests) .and. (sum(abs(int(not(lpassed_nsde))))+sum(abs(int(not(lpassed_ewde)))) == 0)) then
        write(nout,'("Error: all derivatives are correct but all fields zreel, zgp2 and zgp3a are wrong for m=",i4," n=",i4)')nzonal,ntotal
      else if((sum(abs(int(not(lpassed)))) == 0) .and. (sum(abs(int(not(lpassed_nsde))))+sum(abs(int(not(lpassed_ewde)))) == 2*ntests)) then
        write(nout,'("Error: all fields zreel, zgp2 and zgp3a are correct but all derivatives are wrong for m=",i4," n=",i4)')nzonal,ntotal
      else if((sum(abs(int(not(lpassed)))) == 0) .and. (sum(abs(int(not(lpassed_nsde)))) == ntests) .and. (sum(abs(int(not(lpassed_ewde)))) == 0)) then
        write(nout,'("Error: all fields and east-west derivatives are correct but all north-south derivatives are wrong for m=",i4," n=",i4)')nzonal,ntotal
      else if(lpassed(1) .and. lpassed_nsde(1) .and. lpassed_ewde(1) .and. (sum(abs(int(not(lpassed(2:)))))+sum(abs(int(not(lpassed_nsde(2:)))))+sum(abs(int(not(lpassed_ewde(2:))))) == 3*(ntests-1))) then
        write(nout,'("Error: first invtrans is correct but all values are wrong in second invtrans call for m=",i4," n=",i4)')nzonal,ntotal
      else
        if(sum(abs(int(not(lpassed))))>0) write(nout,'("Error: "i3," fields out of ",i4," are wrong for m=",i4," n=",i4)')sum(abs(int(not(lpassed)))),ntests,nzonal,ntotal
        if(sum(abs(int(not(lpassed_nsde))))>0) write(nout,'("Error: "i3," north-south derivatives out of ",i4," are wrong for m=",i4," n=",i4)')sum(abs(int(not(lpassed_nsde)))),ntests,nzonal,ntotal
        if(sum(abs(int(not(lpassed_ewde))))>0) write(nout,'("Error: "i3," east-west derivatives out of ",i4," are wrong for m=",i4," n=",i4)')sum(abs(int(not(lpassed_ewde)))),ntests,nzonal,ntotal
      end if
      call flush(nout)
      stop "error check in grid point space not passed"
    end if
    
  end function check_gp_fields
  
  !===================================================================================================

  function check_sp_fields(rtolerance, lwrite_errors, nflevg, nfld, jstep, nzonal, ntotal, nindex, limag, zspsc2, &
    & zspsc2b, zspsc3a, nout, luse_mpi, nsmax, myproc) result(rlmax_error)

    use parkind1, only: jprb
    use mpl_module

    implicit none

    real(kind=jprd), intent(in) :: rtolerance ! jprb would be enough
    logical, intent(in) :: lwrite_errors
    integer(kind=jpim), intent(in) :: jstep, nzonal, ntotal, nflevg, nfld, nindex, nout, nsmax, myproc
    logical, intent(in) :: limag, luse_mpi
    real(kind=jprd), intent(in) :: zspsc2(:,:), zspsc2b(:,:), zspsc3a(:,:,:) ! should be kind jprb
    real(kind=jprd), intent(out) :: rlmax_error
    real(kind=jprd) :: rlmax_error_zspsc2, rlmax_error_zspsc2b, rlmax_error_zspsc3a
    real(kind=jprd) :: rl2_error_zspsc2, rl2_error_zspsc2b, rl2_error_zspsc3a
    real(kind=jprd) :: rlmax_error_zspsc2_local
    integer(kind=jpim) :: index_max, rank_max, nindex_max, j
    real(kind=jprd), allocatable :: zindex(:)

    logical :: lpassed_zspsc2, lpassed_zspsc2b, lpassed_zspsc3a, lpassed_all

    allocate(zindex(size(zspsc2,2)))
    zindex = 0.0
    if(nindex>0) zindex(nindex) = 1.0
    rlmax_error_zspsc2 = 0.0
    ! checking maximum in a loop to also find the index of the maximum:
    do j = 1,size(zspsc2,2)
      if(abs(zspsc2(1,j)-zindex(j))>rlmax_error_zspsc2) then
        rlmax_error_zspsc2 = abs(zspsc2(1,j)-zindex(j))
        index_max = j
      end if
    end do
    rl2_error_zspsc2 = sum((zspsc2(1,:)-zindex(:))**2)
    rlmax_error_zspsc2_local = rlmax_error_zspsc2
    if(nindex>0) then
      rlmax_error_zspsc2b = max(maxval(abs(zspsc2b(:,1:nindex-1  ))),maxval(abs(zspsc2b(:,nindex)  -1.0_jprd)),maxval(abs(zspsc2b(:,nindex+1:  ))))
      rlmax_error_zspsc3a = max(maxval(abs(zspsc3a(:,1:nindex-1,:))),maxval(abs(zspsc3a(:,nindex,:)-1.0_jprd)),maxval(abs(zspsc3a(:,nindex+1:,:))))
      rl2_error_zspsc2b = sum((zspsc2b(:,1:nindex-1  ))**2)+sum((zspsc2b(:,nindex)  -1.0_jprd)**2)+sum((zspsc2b(:,nindex+1:  ))**2)
      rl2_error_zspsc3a = sum((zspsc3a(:,1:nindex-1,:))**2)+sum((zspsc3a(:,nindex,:)-1.0_jprd)**2)+sum((zspsc3a(:,nindex+1:,:))**2)
    else
      rlmax_error_zspsc2b = maxval(abs(zspsc2b(:,:  )))
      rlmax_error_zspsc3a = maxval(abs(zspsc3a(:,:,:)))
      rl2_error_zspsc2b   = sum((zspsc2b(:,:  ))**2)
      rl2_error_zspsc3a   = sum((zspsc3a(:,:,:))**2)
    end if
    if (luse_mpi) then
      call mpl_allreduce(rlmax_error_zspsc2 , 'max', ldreprod=.false.)
      call mpl_allreduce(rlmax_error_zspsc2b, 'max', ldreprod=.false.)
      call mpl_allreduce(rlmax_error_zspsc3a, 'max', ldreprod=.false.)
      call mpl_allreduce(rl2_error_zspsc2   , 'sum', ldreprod=.false.)
      call mpl_allreduce(rl2_error_zspsc2b  , 'sum', ldreprod=.false.)
      call mpl_allreduce(rl2_error_zspsc3a  , 'sum', ldreprod=.false.)
      ! find rank and index of maximum error for zspsc2:
      if(rlmax_error_zspsc2==rlmax_error_zspsc2_local) then
        rank_max = myproc
        nindex_max = nindex
      else
        rank_max = 0
        index_max = 0
        nindex_max = 0
      end if
      call mpl_allreduce(rank_max  , 'max', ldreprod=.false.)
      call mpl_allreduce(index_max , 'max', ldreprod=.false.)
      call mpl_allreduce(nindex_max , 'max', ldreprod=.false.)
    end if
    rl2_error_zspsc2  = rl2_error_zspsc2 /(nsmax**2)*2
    rl2_error_zspsc2b = rl2_error_zspsc2b/(nsmax**2)*2
    rl2_error_zspsc3a = rl2_error_zspsc3a/(nsmax**2)*2
  
    rlmax_error = max(rlmax_error_zspsc2, rlmax_error_zspsc2, rlmax_error_zspsc3a)
    lpassed_zspsc2  = (rlmax_error_zspsc2  < rtolerance)
    lpassed_zspsc2b = (rlmax_error_zspsc2b < rtolerance)
    lpassed_zspsc3a = (rlmax_error_zspsc3a < rtolerance)
    lpassed_all = lpassed_zspsc2 .and. lpassed_zspsc2b .and. lpassed_zspsc3a
    if(lwrite_errors .and. myproc==1) then
      write(60000+nsmax,'(3i4,L4," ┃",3e11.3," |",3e11.3," |",3i8)') nzonal,ntotal,jstep,limag, &
      & rlmax_error_zspsc2,  &
      & rlmax_error_zspsc2b, &
      & rlmax_error_zspsc3a, &
      & rl2_error_zspsc2,    &
      & rl2_error_zspsc2b,   &
      & rl2_error_zspsc3a,   &
      & nindex_max, index_max, rank_max
      call flush(60000+nsmax)
    end if
    if(.not. lpassed_all) then
      if((.not. lpassed_zspsc2) .and. (.not. lpassed_zspsc2b) .and. (.not. lpassed_zspsc3a)) write(nout,'("Error: all spectral fields are wrong for m=",i4," n=",i4)')nzonal,ntotal
      write(nout,'("there must be something wrong in the direct transform")')
      call flush(nout)
      stop "error check in spectral space not passed"
    end if
  end function check_sp_fields

  !===================================================================================================
  
end module analytic_solutions_mod