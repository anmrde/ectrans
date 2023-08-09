module analytic_solutions_mod

  use parkind1, only: jpim, jprd, jprb

  real(kind=jprb) :: z_pi = 4.d0*datan(1.d0)
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
    integer(kind=jpim), dimension(ndgl+n_regions_ew-1,n_regions_ew) :: nsta, nonl
    real(kind=jprd) :: zlat, zlon
    real(kind=jprd) :: rpi = 2.0_jprd*asin(1.0_jprd)
  
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
      iendlon = nonl(ilat,my_region_ew)
      do jlon = istlon, iendlon
        zlon = real(jlon-1,jprd)*2.0_jprd*rpi/real(nloen(jglat),jprd)
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
    allocate(legpolys(nfirstlat:nlastlat, 0:nsmax, 0:nsmax))
    legpolys = 0.0

    do n = 0,nsmax
      do ilat = nfirstlat, nlastlat
        x = zmu(ilat)
        legpolys(ilat, n, n) = (double_factorial(2 * n - 1) * sqrt(1. - x * x) ** n)
      end do
    end do
  
    do n = 1,nsmax
      do ilat = nfirstlat, nlastlat
        x = zmu(ilat)
        legpolys(ilat, n-1, n) = x * (2 * n - 1) * legpolys(ilat, n-1, n-1)
      end do
    end do
  
    do n = 2,nsmax
      do m = 0,n-2
        do ilat = nfirstlat, nlastlat
          x = zmu(ilat)
          legpolys(ilat, m, n) = (x * (2 * n - 1) * legpolys(ilat, m, n-1) - (n + m - 1) * legpolys(ilat, m, n-2)) / (n - m)
        end do
      end do
    end do
  
    do n = 0,nsmax
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
    real(kind=jprd), dimension(0:nsmax, 0:nsmax) :: zfn
    real(kind=jprd) :: zfnn
    integer(kind=jpim) :: jn, jgl, iodd
  
    if(allocated(legpolys)) deallocate(legpolys)
    allocate(legpolys(nfirstlat:nlastlat, 0:nsmax, 0:nsmax))
    legpolys = 0.0

    zfn(0,0) = 2._jprd
    do jn = 1,nsmax
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
  
    call ini_pol_test(nsmax)
    do jgl=nfirstlat, nlastlat
      call supol_test(nsmax, zmu(jgl), zfn, legpolys(jgl,:,:))
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
    allocate(legpolys(nfirstlat:nlastlat, 0:nsmax, 0:nsmax))
    legpolys = 0.0

    call ini_pol_test(nsmax)
    do jgl=nfirstlat, nlastlat
      do km=0,nsmax
        call supolf_test(km, nsmax, zmu(jgl), legpolys(jgl,km,:))
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
  
    implicit none
  
    integer(kind=jpim), intent(in) :: nsmax, ndgl
    real(kind=jprb), dimension(ndgl,(nsmax+2)*(nsmax+3)/2) :: rpnm
    integer(kind=jpim) :: jnm, jn, jninv, jm, ilat
  
    if(allocated(legpolys_ectrans)) deallocate(legpolys_ectrans)
    allocate(legpolys_ectrans(nfirstlat:nlastlat, 0:nsmax, 0:nsmax))
    legpolys_ectrans = 0.0

    call trans_inq(prpnm=rpnm)
    !print*,"rpnm(1,:)=",rpnm(1,:)
    print*,"rpnm: dimension 1: ",lbound(rpnm,1)," to ",ubound(rpnm,1)," dimension 2: ",lbound(rpnm,2)," to ",ubound(rpnm,2)
  
    !stop "debugging in buffer_legendre_polynomials_ectrans"
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
  
    !jgl = 1
    !print*,"buffer-supolf: sinlat=",sinlats(jgl)," ZLFPOL=",legpolys(jgl,0,:)
  
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
    coeff_a = (n+1)*sqrt((n*n-m*m)/(4.0*n*n-1))
    coeff_b = -n*sqrt(((n+1)*(n+1)-m*m)/(4.0*(n+1)*(n+1)-1))
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
  !      sph_analytic(jrof,ibl) = ectrans_init_spherical_harmonic_eastwest_derivative( ktotal, kzonal, gelam(jrof,ibl), gelat(jrof,ibl), kimag)
        sph_analytic(jrof,ibl) = analytic_eastwest_derivative_point( ktotal, kzonal, gelam(jrof,ibl), gelat(jrof,ibl), kimag, legpolys(nlatidxs(jrof,ibl), kzonal, ktotal))
        !sph_analytic(jrof,ibl) = ectrans_init_spherical_harmonic_hardcoded( ktotal, kzonal, gelam(jrof,ibl), gelat(jrof,ibl), kimag)
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
        if(ktotal<nsmax) then
          legpolyp1 = legpolys(nlatidxs(jrof,ibl), kzonal, ktotal+1)
        end if
        if(ktotal>0) then
          legpolym1 = legpolys(nlatidxs(jrof,ibl), kzonal, ktotal-1)
        end if
        !sph_analytic(jrof,ibl) = ectrans_init_spherical_harmonic_northsouth_derivative_hardcoded( ktotal, kzonal, gelam(jrof,ibl), gelat(jrof,ibl), kimag)
        !sph_analytic(jrof,ibl) = ectrans_init_spherical_harmonic_northsouth_derivative( ktotal, kzonal, gelam(jrof,ibl), gelat(jrof,ibl), kimag)
        sph_analytic(jrof,ibl) = analytic_northsouth_derivative_point( ktotal, kzonal, gelam(jrof,ibl), gelat(jrof,ibl), kimag, legpolyp1, legpolym1)
        !sph_analytic(jrof,ibl) = ectrans_init_spherical_harmonic_hardcoded( ktotal, kzonal, gelam(jrof,ibl), gelat(jrof,ibl), kimag)
      end do
    end do
  
  end subroutine compute_analytic_northsouth_derivative
  
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

  function check_lmax_all_fields(rtolerance, lwrite_errors, nzonal, ntotal, zreel, zgp2, zgp3a, zsph_analytic, znsde_analytic, zewde_analytic) result(lpassed)

    implicit none

    real(kind=jprb), intent(in) :: rtolerance
    logical, intent(in) :: lwrite_errors
    integer(kind=jpim), intent(in) :: nzonal, ntotal
    real(kind=jprb), intent(in) :: zreel(:,:,:), zgp2(:,:,:), zgp3a(:,:,:,:)
    real(kind=jprd), intent(in) :: zsph_analytic(:,:), znsde_analytic(:,:), zewde_analytic(:,:)
    real(kind=jprd) :: lmaxrelquo, lmaxnsderelquo, lmaxewderelquo, lmaxrelfac, lmaxnsderelfac, lmaxewderelfac
    real(kind=jprd) :: lmax_error1, lmax_nsde_error1, lmax_ewde_error1
    real(kind=jprd) :: lmax_error2, lmax_nsde_error2, lmax_ewde_error2
    real(kind=jprd) :: lmax_error3, lmax_nsde_error3, lmax_ewde_error3
    logical :: lpassed
    logical :: lpassed1, lpassed_nsde1, lpassed_ewde1
    logical :: lpassed2, lpassed_nsde2, lpassed_ewde2
    logical :: lpassed3, lpassed_nsde3, lpassed_ewde3

    print*,"DEBUGGING: checking errors with tolerance ", rtolerance
    lmaxrelquo = maxval(abs( zsph_analytic(:,:)))
    lmaxnsderelquo = maxval(abs(znsde_analytic(:,:)))
    lmaxewderelquo = maxval(abs(zewde_analytic(:,:)))
    lmaxrelfac = 1.0_jprd
    lmaxnsderelfac = 1.0_jprd
    lmaxewderelfac = 1.0_jprd
    if(    lmaxrelquo>0.0)     lmaxrelfac = 1.0_jprd/    lmaxrelquo
    if(lmaxnsderelquo>0.0) lmaxnsderelfac = 1.0_jprd/lmaxnsderelquo
    if(lmaxewderelquo>0.0) lmaxewderelfac = 1.0_jprd/lmaxewderelquo
    lmax_error1      = maxval(abs(  zreel(:,1,:)- zsph_analytic(:,:)))*lmaxrelfac
    lmax_error2      = maxval(abs(   zgp2(:,1,:)- zsph_analytic(:,:)))*lmaxrelfac
    lmax_error3      = maxval(abs(zgp3a(:,1,1,:)- zsph_analytic(:,:)))*lmaxrelfac
    lmax_nsde_error1 = maxval(abs(  zreel(:,2,:)-znsde_analytic(:,:)))*lmaxnsderelfac
    lmax_nsde_error2 = maxval(abs(   zgp2(:,2,:)-znsde_analytic(:,:)))*lmaxnsderelfac
    lmax_nsde_error3 = maxval(abs(zgp3a(:,1,2,:)-znsde_analytic(:,:)))*lmaxnsderelfac
    lmax_ewde_error1 = maxval(abs(  zreel(:,3,:)-zewde_analytic(:,:)))*lmaxewderelfac
    lmax_ewde_error2 = maxval(abs(   zgp2(:,3,:)-zewde_analytic(:,:)))*lmaxewderelfac
    lmax_ewde_error3 = maxval(abs(zgp3a(:,1,3,:)-zewde_analytic(:,:)))*lmaxewderelfac
    if(lwrite_errors) then
      write(33342,'(i4,i4," ┃",e11.3,e11.3,e11.3," │",e11.3,e11.3,e11.3," │",e11.3,e11.3,e11.3)') nzonal,ntotal, &
      & lmax_error1, &
      & lmax_error2, &
      & lmax_error3, &
      & lmax_nsde_error1, &
      & lmax_nsde_error2, &
      & lmax_nsde_error3, &
      & lmax_ewde_error1, &
      & lmax_ewde_error2, &
      & lmax_ewde_error3
  !   & sqrt(sum((zgp3a(:,1,3,:)-zewde_analytic(:,:))**2)/ngptot)*lmaxewderelfac
    end if

    lpassed1      = (lmax_error1      < rtolerance)
    lpassed_nsde1 = (lmax_error2      < rtolerance)
    lpassed_ewde1 = (lmax_error3      < rtolerance)
    lpassed2      = (lmax_nsde_error1 < rtolerance)
    lpassed_nsde2 = (lmax_nsde_error2 < rtolerance)
    lpassed_ewde2 = (lmax_nsde_error3 < rtolerance)
    lpassed3      = (lmax_ewde_error1 < rtolerance)
    lpassed_nsde3 = (lmax_ewde_error2 < rtolerance)
    lpassed_ewde3 = (lmax_ewde_error3 < rtolerance)
    if((.not. lpassed1) .and. (.not. lpassed2) .and. (.not. lpassed3)) then
    end if
    
    lpassed = (lpassed1 .and. lpassed_nsde1 .and. lpassed_ewde1 .and. &
    &          lpassed2 .and. lpassed_nsde2 .and. lpassed_ewde2 .and. &
    &          lpassed3 .and. lpassed_nsde3 .and. lpassed_ewde3)
  
  end function check_lmax_all_fields
  
  !===================================================================================================
  
end module analytic_solutions_mod