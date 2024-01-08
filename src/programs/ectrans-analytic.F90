! (C) Copyright 2023- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

program transform_test

!
! Spectral transform analytic test
!
! This test checks the results of the inverse and direct spectral transform computed
! by the ectrans library for specific wavenumbers by using analytically computed
! spherical harmonics. This driver is based on the benchmark driver written by Sam Hatfield et al..
! The main differences to the benchmark driver are:
! - removed code used for measuring the runtime
! - using analytically computed spherical harmonics to check the correctness of inverse and direct
!   transform separately for individual wavenumbers
! - loop over multiple wavenumbers. The loop over multiple iterations (replicating the time-steps
!   from the IFS) still exists inside this loop over the wavenumbers
! - always checking two different calls to the library to make sure that no configuration specific
!   code has been introduced which would break consecutive ectrans calls in the IFS
! - check correctness of scalar fields, vector fields, uv-fields and all the derivatives
! - different default values (in particular using full grid by default!)
! - by default only the highest wavenumber is tested to make the default test as quick as possible
! - option to test all wavenumbers by using flag --test-all
! - summary of all errors is always written to files. Extensive error file showing all latitudes
!   where errors exceed tolerance is optional with flag --errorfiles
! - attempt to indicate in which part of the library an error originates (written to stdout in
!   the function check_gp_fields in analytic_solutions_mod.F90)
! - error-files include in the filename information about precision, truncation, grid and number of
!   MPI processes. The purpose of this is to allow comparing the accuracy of different configurations.
!
! The following comments from the benchmark driver are still correct for this version:
!
! 1) One "surface" field is always transformed:
!      zspsc2(1,1:nspec2) <-> zgmvs(1:nproma,1:1,1:ngbplk)
!
! 2) A Multiple "3d" fields are transformed and can be disabled with "--nfld 0"
!
!      zspsc3a(1:nlev,1:nspec2,1:nfld) <-> zgp3a(1:nproma,1:nlev,1:nfld,1:ngpblk)
!
! 3) Optionally a "3d" vorticity/divergence field is transformed to uv (wind) and
!   can be enabled with "--vordiv"
!
!      zspvor(1:nlev,1:nspec2) / zspdiv(1:nlev,1:nspec2) <-> zgpuv(1:nproma,1:nlev,1:2,1:ngpblk)
!
! 4) Optionally scalar derivatives can be computed for the fields described in 1) and 2)
!    This must be enabled with "--scders"
!
! 5) Optionally uv East-West derivate can be computed from vorticity/divergence.
!    This must be enabled with "--vordiv --uvders"
!
!
! Authors of the benchmark version:
!           George Mozdzynski
!           Willem Deconinck
!           Ioan Hadade
!           Sam Hatfield
!
! Author of the analytic solution modifications:
!           Andreas Müller
!

use parkind1, only: jpim, jprb, jprd, jprm
use oml_mod ,only : oml_max_threads
use mpl_module
use yomhook, only : dr_hook_init
use analytic_solutions_mod, only: analytic_init, analytic_end, &
& buffer_legendre_polynomials_supolf, &
& buffer_legendre_polynomials_ectrans, check_legendre_polynomials, &
& compute_analytic_solution, compute_analytic_eastwest_derivative, &
& compute_analytic_northsouth_derivative, gelam, gelat, init_check_fields, &
& close_check_fields, check_gp_fields, check_sp_fields, compute_analytic_uv, &
& compute_analytic_uv_derivative_ew

implicit none

! Number of points in top/bottom latitudes
integer(kind=jpim), parameter :: min_octa_points = 20

real(kind=jprb), dimension(1) :: zmaxerr(5), zerr(5)
real(kind=jprb) :: zmaxerrg

! Output unit numbers
integer(kind=jpim), parameter :: nerr     = 0 ! Unit number for STDERR
integer(kind=jpim), parameter :: nout     = 6 ! Unit number for STDOUT
integer(kind=jpim), parameter :: noutdump = 7 ! Unit number for field output

! Default parameters
integer(kind=jpim) :: nsmax   = 21  ! Spectral truncation
integer(kind=jpim) :: iters   = 2  ! Number of iterations for transform test
integer(kind=jpim) :: nfld    = 1   ! Number of scalar fields 
integer(kind=jpim) :: nlev    = 1   ! Number of vertical levels
integer(kind=jpim) :: nzonal  = -1   ! zonal wavenumber to be tested
integer(kind=jpim) :: ntotal  = -1   ! total wavenumber to be tested
logical            :: limag = .false. ! test imaginary part of spectral data
integer(kind=jpim) :: nflevg
integer(kind=jpim) :: ndgl ! Number of latitudes
integer(kind=jpim) :: nspec2
integer(kind=jpim) :: ngptot
integer(kind=jpim) :: ngptotg
integer(kind=jpim) :: ifld
integer(kind=jpim) :: jroc
integer(kind=jpim) :: jb
integer(kind=jpim) :: nspec2g
integer(kind=jpim) :: i
integer(kind=jpim) :: ja
integer(kind=jpim) :: ib
integer(kind=jpim) :: jprtrv
integer(kind=jpim) :: n_regions_ns
integer(kind=jpim) :: n_regions_ew

integer(kind=jpim), allocatable :: nloen(:), nprcids(:), nlatidxs(:,:)
integer(kind=jpim) :: myproc, jj, jf, ilf, igp
integer :: jstep

real(kind=jprd) :: ztinit, ztloop, timef, ztstepmax, ztstepmin, ztstepavg, ztstepmed
real(kind=jprd) :: ztstepmax1, ztstepmin1, ztstepavg1, ztstepmed1
real(kind=jprd) :: ztstepmax2, ztstepmin2, ztstepavg2, ztstepmed2
real(kind=jprd), allocatable :: ztstep(:), ztstep1(:), ztstep2(:)

! Grid-point space data structures
real(kind=jprb), allocatable, target :: zgmv   (:,:,:,:) ! Multilevel fields at t and t-dt
real(kind=jprb), allocatable, target :: zgmvs  (:,:,:)   ! Single level fields at t and t-dt
real(kind=jprb), pointer :: zgp3a (:,:,:,:) ! Multilevel fields at t and t-dt
real(kind=jprb), pointer :: zgpuv   (:,:,:,:) ! Multilevel fields at t and t-dt
real(kind=jprb), pointer :: zgp2 (:,:,:) ! Single level fields at t and t-dt

! Spectral space data structures
real(kind=jprb), allocatable, target :: sp3d(:,:,:)
real(kind=jprb), pointer :: zspvor(:,:) => null()
real(kind=jprb), pointer :: zspdiv(:,:) => null()
real(kind=jprb), pointer :: zspsc3a(:,:,:) => null()
real(kind=jprb), allocatable :: zspsc2(:,:), zspsc2b(:,:)
real(kind=jprb), allocatable :: zreel(:,:,:)
real(kind=jprd), allocatable :: zsph_analytic(:,:)
real(kind=jprd), allocatable :: zu_analytic(:,:)
real(kind=jprd), allocatable :: zv_analytic(:,:)
real(kind=jprd), allocatable :: zuder_analytic(:,:)
real(kind=jprd), allocatable :: zvder_analytic(:,:)
real(kind=jprd), allocatable :: zewde_analytic(:,:)
real(kind=jprd), allocatable :: znsde_analytic(:,:)
real(kind=jprd), allocatable :: zsinlats(:)
real(kind=jprd), allocatable :: zlegpolys(:,:,:)
real(kind=jprd), allocatable :: zlegpolys2(:,:,:)

logical :: luserpnm = .true.
logical :: lkeeprpnm = .true.
logical :: luseflt = .false. ! Use fast legendre transforms
logical :: lfftw = .true. ! Use FFTW for Fourier transforms
logical :: lvordiv = .false.
logical :: lscders = .false.
logical :: luvders = .false.
logical :: ltest_all = .false.
logical :: lwrite_errors = .true. ! Write error files for all tested wavenumbers
real(kind=jprd) :: rtolerance = 1e-9 ! maximum relative lmax error tolerance for
                                     ! passing analytyic solution tests in double precision
                                     ! Value for single precision is set at the beginning
                                     ! of the execution.
integer(kind=jpim) :: nstats_mem = 0
integer(kind=jpim) :: ntrace_stats = 0
integer(kind=jpim) :: nprnt_stats = 1

logical :: lmpoff = .false. ! Message passing switch

! Verbosity level (-1, 0 or 1)
integer :: verbosity = -1

real(kind=jprb) :: zra = 6371229._jprb
real(kind=jprb) :: z_pi = 4.d0*datan(1.d0)

integer(kind=jpim) :: nmax_resol = 37 ! Max number of resolutions
integer(kind=jpim) :: npromatr = 0 ! nproma for trans lib
integer(kind=jpim) :: ncombflen = 1800000 ! Size of comm buffer

integer(kind=jpim) :: nproc ! Number of procs
integer(kind=jpim) :: nthread
integer(kind=jpim) :: nprgpns ! Grid-point decomp
integer(kind=jpim) :: nprgpew ! Grid-point decomp
integer(kind=jpim) :: nprtrv = 1 ! Spectral decomp
integer(kind=jpim) :: nprtrw = 0 ! Spectral decomp
integer(kind=jpim) :: nspecresmin = 80 ! Minimum spectral resolution, for controlling nprtrw
integer(kind=jpim) :: mysetv
integer(kind=jpim) :: mysetw
integer(kind=jpim) :: mp_type = 2 ! Message passing type
integer(kind=jpim) :: mbx_size = 150000000 ! Mailbox size
integer(kind=jpim) :: nindex

integer(kind=jpim), allocatable :: numll(:), ivset(:)
integer(kind=jpim) :: ivsetsc(1), ivsetsc1(1)

integer(kind=jpim) :: nflevl
integer(kind=jpim) :: nlats

! sumpini
integer(kind=jpim) :: isqr
logical :: lsync_trans = .true. ! Activate barrier sync
logical :: leq_regions = .true. ! Eq regions flag


integer(kind=jpim) :: nproma = 0
integer(kind=jpim) :: ngpblks
! locals
integer(kind=jpim) :: iprtrv
integer(kind=jpim) :: iprtrw
integer(kind=jpim) :: iprused, ilevpp, irest, ilev, jlev, iprev

integer(kind=jpim) :: ndimgmv  = 0 ! Third dim. of gmv "(nproma,nflevg,ndimgmv,ngpblks)"
integer(kind=jpim) :: ndimgmvs = 0 ! Second dim. gmvs "(nproma,ndimgmvs,ngpblks)"

integer(kind=jpim) :: jbegin_uv = 0
integer(kind=jpim) :: jend_uv   = 0
integer(kind=jpim) :: jbegin_sc = 0
integer(kind=jpim) :: jend_sc   = 0
integer(kind=jpim) :: jbegin_scder_NS = 0
integer(kind=jpim) :: jend_scder_NS = 0
integer(kind=jpim) :: jbegin_scder_EW = 0
integer(kind=jpim) :: jend_scder_EW = 0
integer(kind=jpim) :: jbegin_uder_EW = 0
integer(kind=jpim) :: jend_uder_EW = 0
integer(kind=jpim) :: jbegin_vder_EW = 0
integer(kind=jpim) :: jend_vder_EW = 0
real(kind=jprd) :: rlmax_error_leg, rlmax_error_inv, rlmax_error_dir
logical :: li
integer(kind=jpim) :: m, n, imag_idx

logical :: ldump_values = .false.

integer, external :: ec_mpirank
logical :: luse_mpi = .true.

character(len=16) :: cgrid = ''
character(len=50) :: cBinID = ''  ! binary ID (normally everything after "-analytic-"). This
                                  ! is used to create unique file names for the error files
                                  ! which include information about usage of CPU/GPU, OpenACC/OpenMP, ...

!===================================================================================================

#include "setup_trans0.h"
#include "setup_trans.h"
#include "inv_trans.h"
#include "dir_trans.h"
#include "trans_inq.h"
#include "abor1.intfb.h"

!===================================================================================================

luse_mpi = detect_mpirun()
if(jprb==jprm) rtolerance = 1e-3 ! tolerance for single precision
! Setup
call get_command_line_arguments(nsmax, cgrid, iters, nfld, nlev, lvordiv, lscders, luvders, &
  & luseflt, nproma, verbosity, ltest_all, lwrite_errors, nprtrv, nprtrw, ntotal, nzonal, limag, &
  & rtolerance, cBinID)
if (cgrid == '') cgrid = cubic_full_grid(nsmax)!cubic_octahedral_gaussian_grid(nsmax)
call parse_grid(cgrid, ndgl, nloen)
nflevg = nlev
if(nzonal<0) nzonal = nsmax
if(ntotal<0) ntotal = nsmax

!===================================================================================================

if (luse_mpi) then
  call mpl_init(ldinfo=(verbosity>=1))
  nproc  = mpl_nproc()
  myproc = mpl_myrank()
else
  nproc = 1
  myproc = 1
  mpl_comm = -1
endif
nthread = oml_max_threads()

!===================================================================================================

! only output to stdout on pe 1
if (nproc > 1) then
  if (myproc /= 1) then
    open(unit=nout, file='/dev/null')
  endif
endif

!===================================================================================================

allocate(nprcids(nproc))
do jj = 1, nproc
  nprcids(jj) = jj
enddo

if (nproc <= 1) then
  lmpoff = .true.
endif

! Compute nprgpns and nprgpew
! This version selects most square-like distribution
! These will change if leq_regions=.true.
if (nproc == 0) nproc = 1
isqr = int(sqrt(real(nproc,jprb)))
do ja = isqr, nproc
  ib = nproc/ja
  if (ja*ib == nproc) then
    nprgpns = max(ja,ib)
    nprgpew = min(ja,ib)
    exit
  endif
enddo

! From sumpini, although this should be specified in namelist
if (nspecresmin == 0) nspecresmin = nproc

! Compute nprtrv and nprtrw if not provided on the command line
if (nprtrv > 0 .or. nprtrw > 0) then
  if (nprtrv == 0) nprtrv = nproc/nprtrw
  if (nprtrw == 0) nprtrw = nproc/nprtrv
  if (nprtrw*nprtrv /= nproc) call abor1('transform_test:nprtrw*nprtrv /= nproc')
  if (nprtrw > nspecresmin) call abor1('transform_test:nprtrw > nspecresmin')
else
  do jprtrv = 4, nproc
    nprtrv = jprtrv
    nprtrw = nproc/nprtrv
    if (nprtrv*nprtrw /= nproc) cycle
    if (nprtrv > nprtrw) exit
    if (nprtrw > nspecresmin) cycle
    if (nprtrw <= nspecresmin/(2*oml_max_threads())) exit
  enddo
  ! Go for approx square partition for backup
  if (nprtrv*nprtrw /= nproc .or. nprtrw > nspecresmin .or. nprtrv > nprtrw) then
    isqr = int(sqrt(real(nproc,jprb)))
    do ja = isqr, nproc
      ib = nproc/ja
      if (ja*ib == nproc) then
        nprtrw = max(ja, ib)
        nprtrv = min(ja, ib)
        if (nprtrw > nspecresmin ) then
          call abor1('transform_test:nprtrw (approx square value) > nspecresmin')
        endif
        exit
      endif
    enddo
  endif
endif

! Create communicators for mpi groups
if (.not.lmpoff) then
  call mpl_groups_create(nprtrw, nprtrv)
endif

if (lmpoff) then
  mysetw = (myproc - 1)/nprtrv + 1
  mysetv = mod(myproc - 1, nprtrv) + 1
else
  call mpl_cart_coords(myproc, mysetw, mysetv)

  ! Just checking for now...
  iprtrv = mod(myproc - 1, nprtrv) + 1
  iprtrw = (myproc - 1)/nprtrv + 1
  if (iprtrv /= mysetv .or. iprtrw /= mysetw) then
    call abor1('transform_test:inconsistency when computing mysetw and mysetv')
  endif
endif

if (.not. lmpoff) then
  call mpl_buffer_method(kmp_type=mp_type, kmbx_size=mbx_size, kprocids=nprcids, ldinfo=(verbosity>=1))
endif

! Determine number of local levels for fourier and legendre calculations
! based on the values of nflevg and nprtrv
allocate(numll(nprtrv+1))

! Calculate remainder
iprused = min(nflevg+1, nprtrv)
ilevpp = nflevg/nprtrv
irest = nflevg -ilevpp*nprtrv
do jroc = 1, nprtrv
  if (jroc <= irest) then
    numll(jroc) = ilevpp+1
  else
    numll(jroc) = ilevpp
  endif
enddo
numll(iprused+1:nprtrv+1) = 0

nflevl = numll(mysetv)

ivsetsc(1)  = iprused
ivsetsc1(1) = nprtrv
ifld = 0

!===================================================================================================
! Call ecTrans setup routines
!===================================================================================================

if (verbosity >= 1) write(nout,'(a)')'======= Setup ecTrans ======='

call setup_trans0(kout=nout, kerr=nerr, kprintlev=merge(2, 0, verbosity == 1),                &
  &               kmax_resol=nmax_resol, kpromatr=npromatr, kprgpns=nprgpns, kprgpew=nprgpew, &
  &               kprtrw=nprtrw, kcombflen=ncombflen, ldsync_trans=lsync_trans,               &
  &               ldeq_regions=leq_regions, prad=zra, ldalloperm=.true.,                      &
  &               ldmpoff=.not.luse_mpi, k_regions_ns=n_regions_ns, k_regions_ew=n_regions_ew)

call set_ectrans_gpu_nflev(nflevl)
  ! We pass nflevl via environment variable in order not to change API
  ! In long run, ectrans should grow its internal buffers automatically
call setup_trans(ksmax=nsmax, kdgl=ndgl, kloen=nloen, ldsplit=.true.,          &
  &                 ldusefftw=lfftw, lduserpnm=luserpnm, ldkeeprpnm=lkeeprpnm, &
  &                 lduseflt=luseflt)

call trans_inq(kspec2=nspec2, kspec2g=nspec2g, kgptot=ngptot, kgptotg=ngptotg)

if (nproma == 0) then ! no blocking (default when not specified)
  nproma = ngptot
endif

! Calculate number of NPROMA blocks
ngpblks = (ngptot - 1)/nproma+1

!===================================================================================================
! Print information before starting
!===================================================================================================

! Print configuration details
if (verbosity >= 0) then
  write(nout,'(" ")')
  write(nout,'(a)')'======= Start of runtime parameters ======='
  write(nout,'(" ")')
  write(nout,'("nsmax     ",i0)') nsmax
  write(nout,'("grid      ",a)') trim(cgrid)
  write(nout,'("ndgl      ",i0)') ndgl
  write(nout,'("nproc     ",i0)') nproc
  write(nout,'("nthread   ",i0)') nthread
  write(nout,'("nprgpns   ",i0)') nprgpns
  write(nout,'("nprgpew   ",i0)') nprgpew
  write(nout,'("nprtrw    ",i0)') nprtrw
  write(nout,'("nprtrv    ",i0)') nprtrv
  write(nout,'("ngptot    ",i0)') ngptot
  write(nout,'("ngptotg   ",i0)') ngptotg
  write(nout,'("nfld      ",i0)') nfld
  write(nout,'("nlev      ",i0)') nlev
  write(nout,'("nproma    ",i0)') nproma
  write(nout,'("ngpblks   ",i0)') ngpblks
  write(nout,'("nspec2    ",i0)') nspec2
  write(nout,'("nspec2g   ",i0)') nspec2g
  write(nout,'("nzonal    ",i0)') nzonal
  write(nout,'("ntotal    ",i0)') ntotal
  write(nout,'("iters     ",i0)') iters
  write(nout,'("limag     ",l)') limag
  write(nout,'("luseflt   ",l)') luseflt
  write(nout,'("lvordiv   ",l)') lvordiv
  write(nout,'("lscders   ",l)') lscders
  write(nout,'("luvders   ",l)') luvders
  write(nout,'(" ")')
  write(nout,'(a)') '======= End of runtime parameters ======='
  write(nout,'(" ")')
end if

!===================================================================================================
! Allocate and Initialize spectral arrays
!===================================================================================================

! Allocate spectral arrays
! Try to mimick IFS layout as much as possible
nullify(zspvor)
nullify(zspdiv)
nullify(zspsc3a)
allocate(sp3d(nflevl,nspec2,2+nfld))
allocate(zspsc2(1,nspec2),zspsc2b(1,nspec2))

! Point convenience variables to storage variable sp3d
zspvor  => sp3d(:,:,1)
zspdiv  => sp3d(:,:,2)
zspsc3a => sp3d(:,:,3:3+(nfld-1))

!===================================================================================================
! Allocate gridpoint arrays
!===================================================================================================

allocate(ivset(nflevg))

! Compute spectral distribution
ilev = 0
do jb = 1, nprtrv
  do jlev=1, numll(jb)
    ilev = ilev + 1
    ivset(ilev) = jb
  enddo
enddo

! Allocate grid-point arrays
if (lvordiv) then
  jbegin_uv = 1
  jend_uv = 2
endif
if (luvders) then
  jbegin_uder_EW  = jend_uv + 1
  jend_uder_EW    = jbegin_uder_EW
  jbegin_vder_EW  = jend_uder_EW + 1
  jend_vder_EW    = jbegin_vder_EW
else
  jbegin_uder_EW = jend_uv
  jend_uder_EW   = jend_uv
  jbegin_vder_EW = jend_uv
  jend_vder_EW   = jend_uv
endif

jbegin_sc = jend_vder_EW + 1
jend_sc   = jend_vder_EW + nfld

if (lscders) then
  ndimgmvs = 3
  jbegin_scder_NS = jend_sc + 1
  jend_scder_NS   = jend_sc + nfld
  jbegin_scder_EW = jend_scder_NS + 1
  jend_scder_EW   = jend_scder_NS + nfld
else
  ndimgmvs = 1
  jbegin_scder_NS = jend_sc
  jend_scder_NS   = jend_sc
  jbegin_scder_EW = jend_sc
  jend_scder_EW   = jend_sc
endif

ndimgmv = jend_scder_EW
allocate(zgmv(nproma,nflevg,ndimgmv,ngpblks))
allocate(zgmvs(nproma,ndimgmvs,ngpblks))
allocate(zreel(nproma,3,ngpblks))

zgpuv => zgmv(:,:,1:jend_vder_EW,:)
zgp3a => zgmv(:,:,jbegin_sc:jend_scder_EW,:)
zgp2  => zgmvs(:,:,:)

! Allocate arrays for analytic solutions
allocate(zsph_analytic(nproma,ngpblks),zu_analytic(nproma,ngpblks), &
  & zv_analytic(nproma,ngpblks),zuder_analytic(nproma,ngpblks), &
  & zvder_analytic(nproma,ngpblks),zewde_analytic(nproma,ngpblks), &
  & znsde_analytic(nproma,ngpblks),nlatidxs(nproma,ngpblks),zsinlats(ndgl))
! Compute geographic longitude gelam and latitude gelat and index array nlatidxs(nproma,ngpblks):
call analytic_init(nproma, ngpblks, ndgl, n_regions_ns, n_regions_ew, nloen)
call init_check_fields(lwrite_errors, nsmax, myproc, nproc, cgrid, cBinID)
call buffer_legendre_polynomials_supolf(nsmax)
! Check correctness of Legendre coefficients:
if(nproc==1) then
  ! This test is currently only done for nproc==1. For nproc>1 the ectrans library does not correctly
  ! return all the Legendre coefficients. This should be fixed in the future.
  call buffer_legendre_polynomials_ectrans(nsmax, ndgl)
  rlmax_error_leg = check_legendre_polynomials(rtolerance, lwrite_errors, nsmax, myproc, nproc, ndgl, cgrid, cBinID, nout)
end if

if (iters <= 0) call abor1('transform_test:iters <= 0')

if (verbosity >= 0) then
  write(nout,'(a)') '======= Start of spectral transforms  ======='
  write(nout,'(" ")')
end if

ilf = 0
if(nprtrv == mysetv) then
  ilf = 1
endif
rlmax_error_inv = 0.0
rlmax_error_dir = 0.0

!===================================================================================================
! Perform tests
!===================================================================================================

! Loop over all wavenumbers (check actually tested wavenumber inside)
do n = 0,nsmax
  do m = 0,n
    do imag_idx = 0,1
      li = (imag_idx == 1) ! test imaginary part
      if(((.not. li) .and. (m == 0)) .or. (m>0)) then ! there is no imaginary part for m==0
        if(ltest_all .or. (m == nzonal .and. n == ntotal .and. li == limag)) then ! check if this wavenumber should be tested

          ! Initialize arrays
          zreel(:,:,:)  = 0._jprb
          zgmv(:,:,:,:) = 0._jprb
          zgmvs(:,:,:)  = 0._jprb
          call initialize_spectral_arrays(nsmax, zspsc2, sp3d, m, n, li, nindex)
          zspsc2b = zspsc2

          !=================================================================================================
          ! Compute analytic solutions
          !=================================================================================================

          call compute_analytic_solution(nproma, ngpblks, nsmax, ngptot, m, n, li, zsph_analytic)
          if (lscders) then
            call compute_analytic_eastwest_derivative(nproma, ngpblks, nsmax, ngptot, m, n, li, zewde_analytic)
            call compute_analytic_northsouth_derivative(nproma, ngpblks, nsmax, ngptot, m, n, li, znsde_analytic)
          end if
          if (lvordiv) call compute_analytic_uv(nproma, ngpblks, nsmax, ngptot, m, n, li, zu_analytic, zv_analytic)
          if (luvders) call compute_analytic_uv_derivative_ew(nproma, ngpblks, nsmax, ngptot, m, n, li, zuder_analytic, zvder_analytic)

          !=================================================================================================
          ! Loop over multiple iterations (to see how the errors grow over multiple timesteps)
          !=================================================================================================

          do jstep = 1, iters

            !=================================================================================================
            ! Do inverse transform
            !=================================================================================================

            ! single transform
            call inv_trans(kresol=1, kproma=nproma,   &
              & pspscalar=zspsc2b(1:ilf,:),           & ! spectral scalar
              & ldscders=.true.,                      & ! scalar derivatives
              & kvsetsc=ivsetsc1,                     &
              & pgp=zreel)
            if (lvordiv) then
              ! full time step
              call inv_trans(kresol=1, kproma=nproma, &
                & pspsc2=zspsc2,                      & ! spectral surface pressure
                & pspvor=zspvor,                      & ! spectral vorticity
                & pspdiv=zspdiv,                      & ! spectral divergence
                & pspsc3a=zspsc3a,                    & ! spectral scalars
                & ldscders=lscders,                   &
                & ldvorgp=.false.,                    & ! no gridpoint vorticity
                & lddivgp=.false.,                    & ! no gridpoint divergence
                & lduvder=luvders,                    &
                & kvsetuv=ivset,                      &
                & kvsetsc2=ivsetsc,                   &
                & kvsetsc3a=ivset,                    &
                & pgp2=zgp2,                          &
                & pgpuv=zgpuv,                        &
                & pgp3a=zgp3a)
            else
              call inv_trans(kresol=1, kproma=nproma, &
                  & pspsc2=zspsc2,                    & ! spectral surface pressure
                  & pspsc3a=zspsc3a,                  & ! spectral scalars
                  & ldscders=lscders,                 & ! scalar derivatives
                  & kvsetsc2=ivsetsc,                 &
                  & kvsetsc3a=ivset,                  &
                  & pgp2=zgp2,                        &
                  & pgp3a=zgp3a)
            endif

            ! Compute errors by comparing results from inv_trans with analytic solutions:
            rlmax_error_inv = max(rlmax_error_inv, check_gp_fields(rtolerance, lwrite_errors, nflevg, &
              & nfld, jstep, m, n, li, real(zreel,kind=jprd), real(zgp2,kind=jprd), &
              & real(zgp3a,kind=jprd), real(zgpuv,kind=jprd), zsph_analytic, znsde_analytic, &
              & zewde_analytic, zu_analytic, zv_analytic, zuder_analytic, zvder_analytic, nout, &
              & nsmax, luse_mpi, ngptotg, lscders, lvordiv, luvders, myproc, nproc, cgrid, cBinID))

            !=================================================================================================
            ! Do direct transform
            !=================================================================================================

            call dir_trans(kresol=1, kproma=nproma,   &
              & pgp=zreel(:,:,:),                     &
              & pspscalar=zspsc2b(1:ilf,:),           & ! spectral scalar
              & kvsetsc=ivsetsc1)
            if (lvordiv) then
              call dir_trans(kresol=1, kproma=nproma, &
                & pgp2=zgmvs(:,1:1,:),                &
                & pgpuv=zgpuv(:,:,1:2,:),             &
                & pgp3a=zgp3a(:,:,1:nfld,:),          &
                & pspvor=zspvor,                      &
                & pspdiv=zspdiv,                      &
                & pspsc2=zspsc2,                      &
                & pspsc3a=zspsc3a,                    &
                & kvsetuv=ivset,                      &
                & kvsetsc2=ivsetsc,                   &
                & kvsetsc3a=ivset)
            else
              call dir_trans(kresol=1, kproma=nproma, &
                & pgp2=zgmvs(:,1:1,:),                &
                & pgp3a=zgp3a(:,:,1:nfld,:),          &
                & pspsc2=zspsc2,                      &
                & pspsc3a=zspsc3a,                    &
                & kvsetsc2=ivsetsc,                   &
                & kvsetsc3a=ivset)
            endif

            ! Compute errors by comparing results from dir_trans with analytic initial values:
            rlmax_error_dir = max(rlmax_error_dir, check_sp_fields(rtolerance, lwrite_errors, nflevg, &
              & nfld, jstep, m, n, nindex, li, real(zspsc2,kind=jprd), real(zspsc2b,kind=jprd), &
              & real(zspsc3a,kind=jprd), nout, luse_mpi, nsmax, myproc, nproc))            
          enddo
          write(nout,'("m=",i4," n=",i4," imag=",l1," done")')m,n,li
        end if
      end if
    end do
  end do
end do
write(nout,'("All tests finished.")')
! Write maximum of all errors to stdout:
if(nproc==1) then
  ! This test is currently only done for nproc==1. For nproc>1 the ectrans library does not correctly
  ! return all the Legendre coefficients. This should be fixed in the future.
  write(nout,'("Maximum relative error of Legendre coefficients: ",e11.3)')rlmax_error_leg
end if
write(nout,'("Maximum relative error after invtrans: ",e11.3)')rlmax_error_inv
write(nout,'("Maximum relative error after dirtrans: ",e11.3)')rlmax_error_dir
call flush(nout)
call close_check_fields(lwrite_errors, nsmax, myproc)

!===================================================================================================
! Cleanup
!===================================================================================================

deallocate(zreel)
deallocate(zgmv)
deallocate(zgmvs)

!===================================================================================================
! Finalize MPI
!===================================================================================================

if (luse_mpi) then
  call mpl_end(ldmeminfo=.false.)
endif

!===================================================================================================
! Close file
!===================================================================================================

if (nproc > 1) then
  if (myproc /= 1) then
    close(unit=nout)
  endif
endif

!===================================================================================================

contains

!===================================================================================================

subroutine parse_grid(cgrid,ndgl,nloen)

  character(len=*) :: cgrid
  integer, intent(inout) :: ndgl
  integer, intent(inout), allocatable :: nloen(:)
  integer :: ios
  integer :: gaussian_number
  read(cgrid(2:len_trim(cgrid)),*,IOSTAT=ios) gaussian_number
  if (ios==0) then
    ndgl = 2 * gaussian_number
    allocate(nloen(ndgl))
    if (cgrid(1:1) == 'F') then ! Regular Gaussian grid
      nloen(:) = gaussian_number * 4
      return
    endif
    if (cgrid(1:1) == 'O') then ! Octahedral Gaussian grid
      do i = 1, ndgl / 2
        nloen(i) = 20 + 4 * (i - 1)
        nloen(ndgl - i + 1) = nloen(i)
      end do
      return
    endif
  endif
  call parsing_failed("ERROR: Unsupported grid specified: "// trim(cgrid))

end subroutine

!===================================================================================================

function get_real_value(cname, iarg) result(value)

  real :: value
  character(len=*), intent(in) :: cname
  integer, intent(inout) :: iarg
  character(len=128) :: carg
  integer :: stat

  carg = get_str_value(cname, iarg)
  call str2real(carg, value, stat)

  if (stat /= 0) then
    call parsing_failed("Invalid argument for " // trim(cname) // ": " // trim(carg))
  end if

end function

!===================================================================================================

function get_int_value(cname, iarg) result(value)

  integer :: value
  character(len=*), intent(in) :: cname
  integer, intent(inout) :: iarg
  character(len=128) :: carg
  integer :: stat

  carg = get_str_value(cname, iarg)
  call str2int(carg, value, stat)

  if (stat /= 0) then
    call parsing_failed("Invalid argument for " // trim(cname) // ": " // trim(carg))
  end if

end function

!===================================================================================================

function get_str_value(cname, iarg) result(value)

  character(len=128) :: value
  character(len=*), intent(in) :: cname
  integer, intent(inout) :: iarg

  iarg = iarg + 1
  call get_command_argument(iarg, value)

  if (value == "") then
    call parsing_failed("Invalid argument for " // trim(cname) // ": no value provided")
  end if

end function

!===================================================================================================

subroutine parsing_failed(message)

  character(len=*), intent(in) :: message
  if (luse_mpi) call mpl_init(ldinfo=.false.)
  if (ec_mpirank() == 0) then
    write(nerr,"(a)") trim(message)
    call print_help(unit=nerr)
  endif
  if (luse_mpi) call mpl_end(ldmeminfo=.false.)
  stop

end subroutine

!===================================================================================================

subroutine get_command_line_arguments(nsmax, cgrid, iters, nfld, nlev, lvordiv, lscders, luvders, &
  &                                   luseflt, nproma, verbosity, ltest_all, lwrite_errors, &
  &                                   nprtrv, nprtrw, ntotal, nzonal, limag, rtolerance, cBinID)

  integer, intent(inout) :: nsmax           ! Spectral truncation
  character(len=16), intent(inout) :: cgrid ! Grid
  integer, intent(inout) :: iters           ! Number of iterations for transform test
  integer, intent(inout) :: nfld            ! Number of scalar fields
  integer, intent(inout) :: nlev            ! Number of vertical levels
  logical, intent(inout) :: lvordiv         ! Also transform vorticity/divergence
  logical, intent(inout) :: lscders         ! Compute scalar derivatives
  logical, intent(inout) :: luvders         ! Compute uv East-West derivatives
  logical, intent(inout) :: luseflt         ! Use fast Legendre transforms
  integer, intent(inout) :: nproma          ! NPROMA
  integer, intent(inout) :: verbosity       ! Level of verbosity
  logical, intent(inout) :: ltest_all       ! Test all wavenumbers up to truncation
  logical, intent(inout) :: lwrite_errors   ! Write error files for all tested wavenumbers
  integer, intent(inout) :: nprtrv          ! Size of V set (spectral decomposition)
  integer, intent(inout) :: nprtrw          ! Size of W set (spectral decomposition)
  integer, intent(inout) :: ntotal          ! total wavenumber to be tested
  integer, intent(inout) :: nzonal          ! zonal wavenumber to be tested
  logical, intent(inout) :: limag           ! test imaginary part
  real(jprd), intent(inout) :: rtolerance      ! relative error tolerance for analytic solutions
  character(len=50), intent(out) :: cBinID  ! binary ID (normally everything after "-analytic-"). This
                                            ! is used to create unique file names for the error files
                                            ! which include information about usage of CPU/GPU, OpenACC/OpenMP, ...

  character(len=128) :: carg          ! Storage variable for command line arguments
  integer            :: iarg = 1      ! Argument index
  integer            :: stat          ! For storing success status of string->integer conversion
  integer            :: myproc
  integer            :: idx
  character(len=20) :: searchstring

  ! Extract information from the binary name and return it in cBinID:
  call get_command_argument(0, carg)
  searchstring = "-analytic-"
  idx = index(carg, trim(searchstring), .true.)
  if (idx==0) then ! did not find "-analytic-" => take the entire binary name
    searchstring = "/"
    idx = index(carg, trim(searchstring), .true.)
    ! if there is no "/" => fine, we can take the entire binary name
  end if
  cBinID = trim(carg(idx+len(trim(searchstring)):))
  ! And just in case someone added something like ".exe" to the binary name:
  searchstring = "."
  idx = index(cBinID, trim(searchstring), .false.)
  if (idx>2) cBinID = cBinID(1:idx-1)

  do while (iarg <= command_argument_count())
    call get_command_argument(iarg, carg)

    select case(carg)
      ! Parse help argument
      case('-h', '--help')
        if (luse_mpi) call mpl_init(ldinfo=.false.)
        if (ec_mpirank()==0) call print_help()
        if (luse_mpi) call mpl_end(ldmeminfo=.false.)
        stop
      ! Parse verbosity argument
      case('-v')
        verbosity = 0
      ! Parse number of iterations argument
      case('-n', '--niter')
        iters = get_int_value('-n', iarg)
        if (iters < 1) then
          call parsing_failed("Invalid argument for -n: must be > 0")
        end if
      ! Parse spectral truncation argument
      case('-t', '--truncation')
        nsmax = get_int_value('-t', iarg)
        if (nsmax < 1) then
          call parsing_failed("Invalid argument for -t: must be > 0")
        end if
      case('-g', '--grid'); cgrid = get_str_value('-g', iarg)
      case('-f', '--nfld'); nfld = get_int_value('-f', iarg)
      case('-l', '--nlev'); nlev = get_int_value('-l', iarg)
      case('--nzonal'); nzonal = get_int_value('--nzonal', iarg)
      case('--ntotal'); ntotal = get_int_value('--ntotal', iarg)
      case('--imaginary'); limag = .true.
      case('--vordiv'); lvordiv = .True.
      case('--scders'); lscders = .True.
      case('--uvders'); luvders = .True.
      case('--flt'); luseflt = .True.
      case('--nproma'); nproma = get_int_value('--nproma', iarg)
      case('--dump-values'); ldump_values = .true.
      case('--test-all'); ltest_all = .true.
      case('--errorfiles'); lwrite_errors = .true.
      case('--nprtrv'); nprtrv = get_int_value('--nprtrv', iarg)
      case('--nprtrw'); nprtrw = get_int_value('--nprtrw', iarg)
      case('--tolerance'); rtolerance = get_real_value('--tolerance', iarg)
      case default
        call parsing_failed("Unrecognised argument: " // trim(carg))

    end select
    iarg = iarg + 1
  end do

  if (.not. lvordiv) then
    luvders = .false.
  endif

end subroutine get_command_line_arguments

!===================================================================================================

function cubic_octahedral_gaussian_grid(nsmax) result(cgrid)

  character(len=16) :: cgrid
  integer, intent(in) :: nsmax
  write(cgrid,'(a,i0)') 'O',nsmax+1

end function

!===================================================================================================

function cubic_full_grid(nsmax) result(cgrid)

  character(len=16) :: cgrid
  integer, intent(in) :: nsmax
  write(cgrid,'(a,i0)') 'F',nsmax+1

end function

!===================================================================================================

subroutine str2int(str, int, stat)

  character(len=*), intent(in) :: str
  integer, intent(out) :: int
  integer, intent(out) :: stat
  read(str, *, iostat=stat) int

end subroutine str2int

!===================================================================================================

subroutine str2real(str, real, stat)

  character(len=*), intent(in) :: str
  real, intent(out) :: real
  integer, intent(out) :: stat
  read(str, *, iostat=stat) real

end subroutine str2real

!===================================================================================================

subroutine print_help(unit)

  integer, optional :: unit
  integer :: nout = 6
  if (present(unit)) then
    nout = unit
  endif

  write(nout, "(a)") ""

  if (jprb == jprd) then
    write(nout, "(a)") "NAME    ectrans-analytic-dp"
  else
    write(nout, "(a)") "NAME    ectrans-analytic-sp"
  end if
  write(nout, "(a)") ""

  write(nout, "(a)") "DESCRIPTION"
  write(nout, "(a)") "        This program tests ecTrans by using analytic solutions."
  write(nout, "(a)") ""

  write(nout, "(a)") "USAGE"
  if (jprb == jprd) then
    write(nout, "(a)") "        ectrans-analytic-dp [options]"
  else
    write(nout, "(a)") "        ectrans-analytic-sp [options]"
  end if
  write(nout, "(a)") ""

  write(nout, "(a)") "OPTIONS"
  write(nout, "(a)") "    -h, --help          Print this message"
  write(nout, "(a)") "    -v                  Run with verbose output"
  write(nout, "(a)") "    -t, --truncation T  Run with this triangular spectral truncation"
  write(nout, "(a)") "                        (default = 21)"
  write(nout, "(a)") "    -g, --grid GRID     Run with this grid. Possible values: O<N>, F<N>"
  write(nout, "(a)") "                        If not specified, O<N> is used with N=truncation+1"
  write(nout, "(a)") "                        (cubic relation)"
  write(nout, "(a)") "    -n, --niter NITER   Run for this many inverse/direct transform"
  write(nout, "(a)") "                        iterations (default = 2)"
  write(nout, "(a)") "    -f, --nfld NFLD     Number of scalar fields (default = 1)"
  write(nout, "(a)") "    -l, --nlev NLEV     Number of vertical levels (default = 1)"
  write(nout, "(a)") "    --nzonal NZONAL     Zonal wavenumber that is tested (default = truncation)"
  write(nout, "(a)") "    --ntotal NTOTAL     Total wavenumber that is tested (default = truncation)"
  write(nout, "(a)") "    --imaginary         Test imaginary part (default = false)"
  write(nout, "(a)") "    --test-all          Test all wavenumbers up to (including) truncation"
  write(nout, "(a)") "                        This overwrites --nzonal, --ntotal and --imaginary"
  write(nout, "(a)") "    --tolerance         Test is passed if largest relative lmax-error is"
  write(nout, "(a)") "                        smaller than this tolerance (real value)"
  write(nout, "(a)") "    --vordiv            Also transform vorticity-divergence to wind (default off)"
  write(nout, "(a)") "    --scders            Compute scalar derivatives (default off)"
  write(nout, "(a)") "    --uvders            Compute uv East-West derivatives (default off). Only"
  write(nout, "(a)") "                        when also --vordiv is given"
  write(nout, "(a)") "    --flt               Run with fast Legendre transforms (default off)"
  write(nout, "(a)") "    --nproma NPROMA     Run with NPROMA (default no blocking: NPROMA=ngptot)"
  write(nout, "(a)") "    --nprtrv            Size of V set in spectral decomposition"
  write(nout, "(a)") "    --nprtrw            Size of W set in spectral decomposition"
  write(nout, "(a)") ""

end subroutine print_help

!===================================================================================================

subroutine initialize_spectral_arrays(nsmax, zsp, sp3d, kzonal, ktotal, kimag, kindex)

  integer,            intent(in)    :: nsmax       ! Spectral truncation
  real(kind=jprb),    intent(inout) :: zsp(:,:)    ! Surface pressure
  real(kind=jprb),    intent(inout) :: sp3d(:,:,:) ! 3D fields
  integer,            intent(in)    :: kzonal      ! Zonal wavenumber
  integer,            intent(in)    :: ktotal      ! Total wavenumber
  logical,            intent(in)    :: kimag       ! test imaginary part
  integer(kind=jpim), intent(out)   :: kindex      ! return index that was set to 1

  integer(kind=jpim) :: nflevl
  integer(kind=jpim) :: nfield

  integer :: i, j

  integer :: index, num_my_zon_wns
  integer, allocatable :: my_zon_wns(:), nasm0(:)

  ! Get zonal wavenumbers this rank is responsible for
  call trans_inq(knump=num_my_zon_wns)
  allocate(my_zon_wns(num_my_zon_wns))
  call trans_inq(kmyms=my_zon_wns)

  ! If rank is responsible for the chosen zonal wavenumber...
  if (any(my_zon_wns == kzonal) ) then
    ! Get array of spectral array addresses (this maps (m, n=m) to array index)
    allocate(nasm0(0:nsmax))
    call trans_inq(kasm0=nasm0)

    ! Find out local array index of chosen spherical harmonic
    index = nasm0(kzonal) + 2 * (ktotal - kzonal)
    if(kimag) index = index + 1

    kindex = index
  else
    kindex = -1
  end if

  nflevl = size(sp3d, 1)
  nfield = size(sp3d, 3)

  ! First initialize surface pressure
  call initialize_2d_spectral_field(zsp(1,:), kindex)

  ! Then initialize all of the 3D fields
  do i = 1, nflevl
    do j = 1, nfield
      call initialize_2d_spectral_field(sp3d(i,:,j), kindex)
    end do
  end do

end subroutine initialize_spectral_arrays

!===================================================================================================

subroutine initialize_2d_spectral_field(field, kindex)

  real(kind=jprb), intent(inout) :: field(:) ! Field to initialize
  integer,         intent(out)   :: kindex   ! return index that is set to one for testing result

  ! First initialise all spectral coefficients to zero
  field(:) = 0.0

  ! If rank is responsible for the chosen zonal wavenumber...
  if (kindex > 0) then
    ! Set just that element to a constant value
    field(kindex) = 1.0
  else
    return
  end if

end subroutine initialize_2d_spectral_field

!===================================================================================================

function detect_mpirun() result(lmpi_required)
  logical :: lmpi_required
  integer :: ilen
  integer, parameter :: nvars = 5
  character(len=32), dimension(nvars) :: cmpirun_detect
  character(len=4) :: clenv_dr_hook_assert_mpi_initialized
  integer :: ivar

  ! Environment variables that are set when mpirun, srun, aprun, ... are used
  cmpirun_detect(1) = 'OMPI_COMM_WORLD_SIZE'  ! openmpi
  cmpirun_detect(2) = 'ALPS_APP_PE'           ! cray pe
  cmpirun_detect(3) = 'PMI_SIZE'              ! intel
  cmpirun_detect(4) = 'SLURM_NTASKS'          ! slurm
  cmpirun_detect(5) = 'ECTRANS_USE_MPI'       ! forced

  lmpi_required = .false.
  do ivar = 1, nvars
    call get_environment_variable(name=trim(cmpirun_detect(ivar)), length=ilen)
    if (ilen > 0) then
      lmpi_required = .true.
      exit ! break
    endif
  enddo
end function

!===================================================================================================

subroutine set_ectrans_gpu_nflev(kflev)
  use ec_env_mod, only : ec_putenv
  integer(kind=jpim), intent(in) :: kflev
  character(len=32) :: ECTRANS_GPU_NFLEV
  write(ECTRANS_GPU_NFLEV,'(A,I0)') "ECTRANS_GPU_NFLEV=",kflev
  call ec_putenv(ECTRANS_GPU_NFLEV, overwrite=.true.)
end subroutine

end program transform_test

!===================================================================================================
