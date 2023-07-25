module analytic_solutions_mod

  use parkind1, only: jpim, jprd, jprb

  real(kind=jprb) :: zra = 6371229._jprb
  real(kind=jprb) :: z_pi = 4.d0*datan(1.d0)

  #include "trans_inq.h"

  contains

  subroutine calc_gelam_gelat(nproma, ngpblks, ndgl, n_regions_ns, n_regions_ew, nloen, gelam, gelat, nlatidxs, sinlats, nlats)

    implicit none
  
    integer(kind=jpim), intent(in) :: nproma, ngpblks, ndgl, n_regions_ns, n_regions_ew
    integer(kind=jpim), dimension(ndgl), intent(in) :: nloen
    real(kind=jprd), dimension(nproma,ngpblks), intent(out) :: gelam, gelat
    real(kind=jprd), dimension(ndgl), intent(out) :: sinlats
    integer(kind=jpim), dimension(nproma,ngpblks), intent(out) :: nlatidxs
    integer(kind=jpim), intent(out) :: nlats
    integer(kind=jpim) :: nptrfloff, my_region_ns, my_region_ew
    integer(kind=jpim) :: jglat, ioff, ilat, istlon, iendlon, jlon, jrof, ibl, ilat0
    integer(kind=jpim), dimension(n_regions_ns) :: nfrstlat, nlstlat
    integer(kind=jpim), dimension(ndgl+n_regions_ew-1,n_regions_ew) :: nsta, nonl
    real(kind=jprd), dimension(ndgl) :: zmu
    real(kind=jprd) :: zlat, zlon
    real(kind=jprd) :: rpi = 2.0_jprd*asin(1.0_jprd)
  
    call trans_inq(kptrfloff=nptrfloff, &
       & kmy_region_ns=my_region_ns,    &
       & kmy_region_ew=my_region_ew,    &
       & kfrstlat=nfrstlat,             &
       & klstlat=nlstlat,               &
       & ksta=nsta,                     &
       & konl=nonl,                     &
       & pmu=zmu)
  
       !DEBUGGING:
       !print*,"calc_gelam: nfrstlat=",nfrstlat(my_region_ns)," nlstlat=",nlstlat(my_region_ns)," istlon=",istlon," iendlon=",iendlon," num=",(nlstlat(my_region_ns)-nfrstlat(my_region_ns)+1)*(iendlon-istlon)," ngptot=",ngptot
  !write(33350,'("calc_gelam: num=",i10," ngptot=",i10)'),sum(nonl(nptrfloff+1:nptrfloff+nlstlat(my_region_ns)-nfrstlat(my_region_ns)+1,my_region_ew)-nsta(nptrfloff+1:nptrfloff+nlstlat(my_region_ns)-nfrstlat(my_region_ns)+1,my_region_ew)+1),ngptot
    ilat = nptrfloff
    ilat0 = 0
    ibl  = 1
    jrof = 1
    do jglat = nfrstlat(my_region_ns), nlstlat(my_region_ns)
      zlat = asin(zmu(jglat))
      ilat = ilat + 1
      ilat0 = ilat0 + 1
      istlon = nsta(ilat,my_region_ew)
      iendlon = nonl(ilat,my_region_ew)
      do jlon = istlon, iendlon
        zlon = real(jlon-1,jprd)*2.0_jprd*rpi/real(nloen(jglat),jprd)
        gelam(jrof,ibl) = zlon
        gelat(jrof,ibl) = zlat
        nlatidxs(jrof,ibl) = ilat0
        sinlats(ilat0) = zmu(jglat)
        jrof = jrof + 1
        if(jrof > nproma) then
          jrof = 1
          ibl  = ibl + 1
        end if
      end do
    end do
    nlats = ilat0
    print*,"calc_gelam: jrof=",jrof
  
    print*,"calc_gelam_gelat finished"
  
  end subroutine calc_gelam_gelat
  
  !===================================================================================================
  
  subroutine buffer_legendre_polynomials(nsmax, nlats, sinlats, legpolys)
  
    implicit none
  
    integer(kind=jpim), intent(in) :: nsmax, nlats
    real(kind=jprd), dimension(nlats), intent(in) :: sinlats
    real(kind=jprd), dimension(nlats, 0:nsmax, 0:nsmax), intent(out) :: legpolys
    real(kind=jprd) :: x
    integer(kind=jpim) :: n, m, ilat
  
    do n = 0,nsmax
      do ilat = 1,nlats
        x = sinlats(ilat)
        legpolys(ilat, n, n) = (double_factorial(2 * n - 1) * sqrt(1. - x * x) ** n)
      end do
    end do
  
    do n = 1,nsmax
      do ilat = 1,nlats
        x = sinlats(ilat)
        legpolys(ilat, n-1, n) = x * (2 * n - 1) * legpolys(ilat, n-1, n-1)
      end do
    end do
  
    do n = 2,nsmax
      do m = 0,n-2
        do ilat = 1,nlats
          x = sinlats(ilat)
          legpolys(ilat, m, n) = (x * (2 * n - 1) * legpolys(ilat, m, n-1) - (n + m - 1) * legpolys(ilat, m, n-2)) / (n - m)
        end do
      end do
    end do
  
    do n = 0,nsmax
      do m = 0,n
        do ilat = 1,nlats
          legpolys(ilat, m, n) = K(n, m) * legpolys(ilat, m, n)
        end do
      end do
    end do
  
  end subroutine buffer_legendre_polynomials
  
  !===================================================================================================
  
  subroutine buffer_legendre_polynomials_belusov(nsmax, nlats, sinlats, legpolys)
  
    use supol_test_mod, only: supol_test
    use tpm_pol_test, only: ini_pol_test

    implicit none
  
    integer(kind=jpim), intent(in) :: nsmax, nlats
    real(kind=jprd), dimension(nlats), intent(in) :: sinlats
    real(kind=jprd), dimension(nlats, 0:nsmax, 0:nsmax), intent(out) :: legpolys
    real(kind=jprd), dimension(0:nsmax, 0:nsmax) :: zfn
    real(kind=jprd) :: zfnn
    integer(kind=jpim) :: jn, jgl, iodd
  
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
    do jgl=1,nlats
      call supol_test(nsmax, sinlats(jgl), zfn, legpolys(jgl,:,:))
    end do
  
    !jgl = 1
    !print*,"buffer-belusov: sinlat=",sinlats(jgl)," ZFN=",zfn(:,0)," ZLFPOL=",legpolys(jgl,0,:)
  
  end subroutine buffer_legendre_polynomials_belusov
  
  !===================================================================================================
  
  subroutine buffer_legendre_polynomials_supolf(nsmax, nlats, sinlats, legpolys)
  
    use supolf_test_mod, only: supolf_test
    use tpm_pol_test, only: ini_pol_test, dfa
  
    implicit none
  
    integer(kind=jpim), intent(in) :: nsmax, nlats
    real(kind=jprd), dimension(nlats), intent(in) :: sinlats
    real(kind=jprd), dimension(nlats, 0:nsmax, 0:nsmax), intent(out) :: legpolys
    integer(kind=jpim) :: km, jgl
    integer(kind=jpim), dimension(nsmax) :: ndglu
    !stop "debugging in buffer_legendre_polynomials_supolf"
    call ini_pol_test(nsmax)
    do jgl=1,nlats
      do km=0,nsmax
        call supolf_test(km, nsmax, sinlats(jgl), legpolys(jgl,km,:))
      end do
    end do
  
    ! remove unused high wavenumbers at high latitudes
    !call trans_inq(kdglu=ndglu)
    !do jm=0,nsmax
    !  do ilat=
    !do ilat=1,nlats
    !  do jm=0,nsmax+1
    !    do jn=jm,nsmax+1
    !      if() legpolys(ilat,jm,jn) = rpnm(ilat,jnm)
    !    end do
    !  end do
    !end do
    !jgl = 1
    !print*,"buffer-supolf: sinlat=",sinlats(jgl)," ZLFPOL=",legpolys(jgl,0,:)
  
  end subroutine buffer_legendre_polynomials_supolf
  
  !===================================================================================================
  
  ! Using ectrans via trans_inq to retrieve the Legendre polynomials
  ! Caution: ectrans only returns the Legendre polynomials in single precision!!! (supol and supolf return them in double precision)
  subroutine buffer_legendre_polynomials_ectrans(nsmax, nlats, ndgl, sinlats, legpolys)
  
    use supolf_test_mod, only: supolf_test
    use tpm_pol_test, only: ini_pol_test, dfa
  
    implicit none
  
    integer(kind=jpim), intent(in) :: nsmax, nlats, ndgl
    real(kind=jprd), dimension(nlats), intent(in) :: sinlats
    real(kind=jprd), dimension(nlats, 0:nsmax, 0:nsmax), intent(out) :: legpolys
    real(kind=jprb), dimension(ndgl,(nsmax+2)*(nsmax+3)/2) :: rpnm
    integer(kind=jpim) :: jnm, jn, jninv, jm, ilat
  
    call trans_inq(prpnm=rpnm)
    !print*,"rpnm(1,:)=",rpnm(1,:)
    print*,"rpnm: dimension 1: ",lbound(rpnm,1)," to ",ubound(rpnm,1)," dimension 2: ",lbound(rpnm,2)," to ",ubound(rpnm,2)
  
    !stop "debugging in buffer_legendre_polynomials_ectrans"
    do ilat=1,nlats
      jnm = 0
      do jm=0,nsmax+1
        do jninv=jm,nsmax+1
          jn = nsmax+1-jninv+jm
          jnm = jnm + 1
          if(jm<=nsmax .and. jn<=nsmax) legpolys(ilat,jm,jn) = rpnm(ilat,jnm)
        end do
      end do
    end do
  
    !jgl = 1
    !print*,"buffer-supolf: sinlat=",sinlats(jgl)," ZLFPOL=",legpolys(jgl,0,:)
  
  end subroutine buffer_legendre_polynomials_ectrans
  
  !===================================================================================================
  
  subroutine check_legendre_polynomials(nsmax, nlats, legpolys, legpolys2)
  
    implicit none
  
    integer(kind=jpim), intent(in) :: nsmax, nlats
    real(kind=jprd), dimension(nlats, 0:nsmax, 0:nsmax), intent(in) :: legpolys, legpolys2
    integer(kind=jpim) :: ilat, jm, jn
    logical :: lprint
  
  
    write(33320,'("Legendre coefficients")')
    write(33320,'("ilat   m   n ┃     supolf    ectrans ")')
    write(33320,'("━━━━━━━━━━━━━╋━━━━━━━━━━━━━━━━━━━━━━┿")')
    do ilat=1,nlats/2
      do jm=0,nsmax
        do jn=jm,nsmax
          lprint = .false.
          !print*,ilat,jm,jn,abs(legpolys(ilat,jm,jn)),abs(legpolys2(ilat,jm,jn))
          if(legpolys(ilat,jm,jn)==legpolys(ilat,jm,jn) .and. legpolys2(ilat,jm,jn)==legpolys2(ilat,jm,jn)) then ! just to be sure that there are no nans
            if(max(abs(legpolys(ilat,jm,jn)),abs(legpolys2(ilat,jm,jn)))>0) then
              if(abs(legpolys(ilat,jm,jn)-legpolys2(ilat,jm,jn))/max(abs(legpolys(ilat,jm,jn)),abs(legpolys2(ilat,jm,jn)))>1E-5) then
                lprint = .true.
              end if
            end if
          else
            !lprint = .true.
          end if
          if(lprint) write(33320,'(i4,i4,i4," ┃",e11.3,e11.3,e11.3,e11.3,e11.3)') ilat, jm, jn, legpolys(ilat,jm,jn), legpolys2(ilat,jm,jn),abs(legpolys(ilat,jm,jn)-legpolys2(ilat,jm,jn)),max(abs(legpolys(ilat,jm,jn)),abs(legpolys2(ilat,jm,jn))),abs(legpolys(ilat,jm,jn)-legpolys2(ilat,jm,jn))/max(abs(legpolys(ilat,jm,jn)),abs(legpolys2(ilat,jm,jn)))
        end do
      end do
    end do
    
  end subroutine check_legendre_polynomials
  
  !===================================================================================================
  
  subroutine compute_analytic_solution(nproma, ngpblks, nlats, nsmax, ngptot, gelam, gelat, kzonal, ktotal, kimag, legpolys, klatidxs, sph_analytic)
  
    implicit none
  
    integer(kind=jpim), intent(in) :: nproma, ngpblks, nlats, nsmax, ngptot
    real(kind=jprd), dimension(nproma,ngpblks), intent(in) :: gelam, gelat
    real(kind=jprd), dimension(nproma,ngpblks), intent(out) :: sph_analytic
    integer(kind=jpim), intent(in) :: kzonal, ktotal
    logical, intent(in) :: kimag
    real(kind=jprd), dimension(nlats, 0:nsmax, 0:nsmax), intent(in) :: legpolys
    integer(kind=jpim), dimension(nproma,ngpblks), intent(in) :: klatidxs
    integer(kind=jpim) :: jkglo, iend, ioff, ibl, jrof
  
    do jkglo=1,ngptot,nproma
      iend = min(nproma,ngptot-jkglo+1)
      ioff = jkglo - 1
      ibl  = (jkglo-1)/nproma+1
      do jrof=1,iend
        !sph_analytic(jrof,ibl) = ectrans_init_spherical_harmonic( ktotal, kzonal, gelam(jrof,ibl), gelat(jrof,ibl), kimag)
        !print *,"C: result=",sph_analytic(jrof,ibl)
        sph_analytic(jrof,ibl) = analytic_spherical_harmonic_point( ktotal, kzonal, gelam(jrof,ibl), gelat(jrof,ibl), kimag, legpolys(klatidxs(jrof,ibl), kzonal, ktotal))
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

    implicit none

    integer(jpim), intent(in), value :: n
    integer(jpim), intent(in), value :: m
    real(jprd), intent(in), value    :: lon
    real(jprd), intent(in), value    :: lat
    logical, intent(in), value       :: imag
    real(jprd), intent(in), value    :: legpoly
    real(jprd) :: pointValue
    if (imag) then
      pointValue = (- m * analytic_spherical_harmonic_point(n, m, lon, lat, .false., legpoly) / (zra * cos(lat)));
    else
      pointValue = (m * analytic_spherical_harmonic_point(n, m, lon, lat, .true., legpoly) / (zra * cos(lat)));
    end if
  
  end function analytic_eastwest_derivative_point
  
  function analytic_northsouth_derivative_point(n, m, lon, lat, imag, legpolyp1, legpolym1) result(pointValue)

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
          pointValue = (coeff_a * legpolym1 + coeff_b * legpolyp1) / (zra * coslat)
        else
          pointValue = (coeff_b * legpolyp1) / (zra * coslat)
        end if
      end if
    end if
  
    if (m > 0) then
      if (n > m) then
        if (imag) then
            pointValue = -2 * sin(m * lon) * (coeff_a * legpolym1 + coeff_b * legpolyp1) / (zra * coslat);
        else
            pointValue = 2 * cos(m * lon) * (coeff_a * legpolym1 + coeff_b * legpolyp1) / (zra * coslat);
        end if
      else
        if (imag) then
            pointValue = -2 * sin(m * lon) * (coeff_b * legpolyp1) / (zra * coslat);
        else
            pointValue = 2 * cos(m * lon) * (coeff_b * legpolyp1) / (zra * coslat);
        end if
      end if
    end if
  
  end function analytic_northsouth_derivative_point
  
  !===================================================================================================
  
  subroutine compute_analytic_eastwest_derivative(nproma, ngpblks, nlats, nsmax, ngptot, gelam, gelat, kzonal, ktotal, kimag, legpolys, klatidxs, sph_analytic)
  
    implicit none
  
    integer(kind=jpim), intent(in) :: nproma, ngpblks, nlats, nsmax, ngptot
    real(kind=jprd), dimension(nproma,ngpblks), intent(in) :: gelam, gelat
    real(kind=jprd), dimension(nproma,ngpblks), intent(out) :: sph_analytic
    integer(kind=jpim), intent(in) :: kzonal, ktotal
    logical, intent(in) :: kimag
    real(kind=jprd), dimension(nlats, 0:nsmax, 0:nsmax), intent(in) :: legpolys
    integer(kind=jpim), dimension(nproma,ngpblks), intent(in) :: klatidxs
    integer(kind=jpim) :: jkglo, iend, ioff, ibl, jrof
  
    do jkglo=1,ngptot,nproma
      iend = min(nproma,ngptot-jkglo+1)
      ioff = jkglo - 1
      ibl  = (jkglo-1)/nproma+1
      do jrof=1,iend
  !      sph_analytic(jrof,ibl) = ectrans_init_spherical_harmonic_eastwest_derivative( ktotal, kzonal, gelam(jrof,ibl), gelat(jrof,ibl), kimag)
        sph_analytic(jrof,ibl) = analytic_eastwest_derivative_point( ktotal, kzonal, gelam(jrof,ibl), gelat(jrof,ibl), kimag, legpolys(klatidxs(jrof,ibl), kzonal, ktotal))
        !sph_analytic(jrof,ibl) = ectrans_init_spherical_harmonic_hardcoded( ktotal, kzonal, gelam(jrof,ibl), gelat(jrof,ibl), kimag)
      end do
    end do
  
  end subroutine compute_analytic_eastwest_derivative
  
  !===================================================================================================
  
  subroutine compute_analytic_northsouth_derivative(nproma, ngpblks, nlats, nsmax, ngptot, gelam, gelat, kzonal, ktotal, kimag, legpolys, klatidxs, sph_analytic)
  
    implicit none
  
    integer(kind=jpim), intent(in) :: nproma, ngpblks, nlats, nsmax, ngptot
    real(kind=jprd), dimension(nproma,ngpblks), intent(in) :: gelam, gelat
    real(kind=jprd), dimension(nproma,ngpblks), intent(out) :: sph_analytic
    integer(kind=jpim), intent(in) :: kzonal, ktotal
    logical, intent(in) :: kimag
    real(kind=jprd), dimension(nlats, 0:nsmax, 0:nsmax), intent(in) :: legpolys
    real(kind=jprd) :: legpolyp1, legpolym1
    integer(kind=jpim), dimension(nproma,ngpblks), intent(in) :: klatidxs
    integer(kind=jpim) :: jkglo, iend, ioff, ibl, jrof
  
    legpolyp1 = 0.0
    legpolym1 = 0.0
  
    do jkglo=1,ngptot,nproma
      iend = min(nproma,ngptot-jkglo+1)
      ioff = jkglo - 1
      ibl  = (jkglo-1)/nproma+1
      do jrof=1,iend
        if(ktotal<nsmax) then
          legpolyp1 = legpolys(klatidxs(jrof,ibl), kzonal, ktotal+1)
        end if
        if(ktotal>0) then
          legpolym1 = legpolys(klatidxs(jrof,ibl), kzonal, ktotal-1)
        end if
        !sph_analytic(jrof,ibl) = ectrans_init_spherical_harmonic_northsouth_derivative_hardcoded( ktotal, kzonal, gelam(jrof,ibl), gelat(jrof,ibl), kimag)
        !sph_analytic(jrof,ibl) = ectrans_init_spherical_harmonic_northsouth_derivative( ktotal, kzonal, gelam(jrof,ibl), gelat(jrof,ibl), kimag)
        sph_analytic(jrof,ibl) = analytic_northsouth_derivative_point( ktotal, kzonal, gelam(jrof,ibl), gelat(jrof,ibl), kimag, legpolyp1, legpolym1)
        !sph_analytic(jrof,ibl) = ectrans_init_spherical_harmonic_hardcoded( ktotal, kzonal, gelam(jrof,ibl), gelat(jrof,ibl), kimag)
      end do
    end do
  
  end subroutine compute_analytic_northsouth_derivative
  
  !===================================================================================================
  
end module analytic_solutions_mod