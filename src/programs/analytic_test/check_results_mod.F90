module check_results_mod

  use parkind1, only: jpim, jprd

  real(kind=jprd) :: z_pi = 4.d0*datan(1.d0)

  contains

  !===================================================================================================
  ! Function check_legendre_polynomials:
  ! Computes the difference between the Legendre coefficients returned by the ectrans library
  ! with the subroutine buffer_legendre_polynomials_ectrans in the array legpolys_ectrans and
  ! the array legpolys computed with the subroutine supolf_test by buffer_legendre_polynomials_supolf.
  ! These arrays need to be called before this function is called! This function returns the maximum
  ! relative difference between the Legendre coefficients computed with the two code paths. For
  ! lwrite_errors==.true. it also writes the coefficients where the error exceeds the tolerance to file.
  ! This function works currently only for nproc==1. For nproc>1 the ectrans library does not correctly
  ! return all Legendre coefficients.
  !
  ! Even if the errors of the Legendre coefficients exceed the tolerance the code currently continues
  ! to run. As long as the errors of inv_trans and dir_trans are small it should be fine. For nproc>1
  ! ectrans returns only some of the coefficients correctly even though the overall results are correct.
  !===================================================================================================
  
  function check_legendre_polynomials(rtolerance, lwrite_errors, nsmax, myproc, nproc, ndgl, &
    & cgrid, cBinID, nout) result(rlmax_error)
  
    use parkind1, only: jprb, jprd
    use analytic_solutions_mod, only: legpolys, legpolys_ectrans, nfirstlat, nlastlat, nmeng

    implicit none
  
    real(kind=jprd), intent(in) :: rtolerance
    logical, intent(in) :: lwrite_errors
    integer(kind=jpim), intent(in) :: nsmax, ndgl
    integer(kind=jpim), intent(in) :: myproc, nproc
    character(len=16), intent(in)  :: cgrid
    real(kind=jprd) :: rlmax_error
    integer(kind=jpim), intent(in) :: nout
    character(len=50), intent(in) :: cBinID   ! binary ID (normally everything in the name of the binary
        ! after "-analytic-"). This is used to create unique file names for the error files
        ! which include information about usage of CPU/GPU, OpenACC/OpenMP, ...
    real(kind=jprd) :: rlmax_quo, rlmax_fac
    character(len=100) :: filename
    integer(kind=jpim) :: ilat, jm, jn
    logical :: lprint, linitialized

    if (myproc == 1) then
      linitialized = .false.
      rlmax_error = 0.0_jprd
      rlmax_quo = maxval(abs(legpolys))
      if(rlmax_quo == 0.0) rlmax_quo = maxval(abs(legpolys_ectrans))
      rlmax_fac = 1.0_jprd
      if(rlmax_quo>0.0) rlmax_fac = 1.0_jprd/rlmax_quo
      do ilat=nfirstlat, nlastlat
        if(ilat<ndgl/2) then
          do jm=0,nmeng(ilat)
            do jn=jm,nmeng(ilat)
              lprint = .false.
              if(legpolys(ilat,jm,jn)==legpolys(ilat,jm,jn) .and. legpolys_ectrans(ilat,jm,jn)==legpolys_ectrans(ilat,jm,jn)) then ! just to be sure that there are no nans
                if(abs(legpolys(ilat,jm,jn)-legpolys_ectrans(ilat,jm,jn))*rlmax_fac > rtolerance) then
                  lprint = .true.
                end if
                if(lprint .and. lwrite_errors) then
                  if(.not. linitialized) then
                    write(filename,'(3a,i0,3a,i0,a)')'errors-',trim(cBinID),'-legendre_T',nsmax,'_',trim(cgrid),'_nproc',nproc,'.txt'
                    open( 30, file = filename, status='replace')
                    write(30,'("Legendre coefficients")')
                    write(30,'("ilat   m   n ┃     supolf    ectrans      error ")')
                    write(30,'("━━━━━━━━━━━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┿")')
                    linitialized = .true.
                  endif            
                  write(30,'(i4,i4,i4," ┃",e11.3,e11.3,e11.3)') ilat, jm, jn, legpolys(ilat,jm,jn), legpolys_ectrans(ilat,jm,jn),abs(legpolys(ilat,jm,jn)-legpolys_ectrans(ilat,jm,jn))*rlmax_fac
                end if
                rlmax_error = max(rlmax_error, abs(legpolys(ilat,jm,jn)-legpolys_ectrans(ilat,jm,jn))*rlmax_fac)
              else
                if(legpolys_ectrans(ilat,jm,jn)==legpolys_ectrans(ilat,jm,jn)) write(nout,'("check_legendre_polynomials: ilat=",i4," jm=",i4," jn=",i4," legpolys_ectrans=",e11.3," Issue with Legendre coefficients computed in analytic_solutions_mod.F90 but will continue.")') ilat, jm, jn, legpolys_ectrans(ilat,jm,jn)
                if(legpolys(ilat,jm,jn)==legpolys(ilat,jm,jn)) write(nout,'("check_legendre_polynomials: ilat=",i4," jm=",i4," jn=",i4," legpolys=",e11.3," Issue with Legendre coefficients computed in ectrans library but will continue.")') ilat, jm, jn, legpolys(ilat,jm,jn)
                flush(nout)
              end if
            end do
          end do
        end if
      end do
    end if
    
  end function check_legendre_polynomials

  !===================================================================================================
  ! Subroutine init_check_fields:
  ! Initialize the files in which the errors are written.
  !===================================================================================================

  subroutine init_check_fields(lwrite_errors, nsmax, myproc, nproc, cgrid, cBinID)

    use parkind1, only: jprb, jprd

    implicit none

    logical, intent(in) :: lwrite_errors
    integer(kind=jpim), intent(in) :: nsmax
    integer(kind=jpim), intent(in) :: myproc, nproc
    character(len=16), intent(in)  :: cgrid
    character(len=50), intent(in) :: cBinID   ! binary ID (normally everything in the name of the binary
        ! after "-analytic-"). This is used to create unique file names for the error files
        ! which include information about usage of CPU/GPU, OpenACC/OpenMP, ...
    character(len=100) :: filename

    if (myproc == 1) then
      write(filename,'(3a,i0,3a,i0,a)')'errors-',trim(cBinID),'-gridpoint_T',nsmax,'_',trim(cgrid),'_nproc',nproc,'.txt'
      open(40, file = filename, status='replace')
      write(filename,'(3a,i0,3a,i0,a)')'errors-',trim(cBinID),'-spectral_T',nsmax,'_',trim(cgrid),'_nproc',nproc,'.txt'
      open(60, file = filename, status='replace')
      if(lwrite_errors) then
        write(40,'("lmax-error in grid point space")')
        write(40,'("                 ┃           grid point data        │      north-south derivative      │       east-west derivative       |       wind speed      |  east-west derivative ")')
        write(40,'("   m   n it. im. ┃        pgp       pgp2      pgp3a │        pgp       pgp2      pgp3a │        pgp       pgp2      pgp3a │          u          v │          u          v ")')
        write(40,'("━━━━━━━━━━━━━━━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┿━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┿━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┿━━━━━━━━━━━━━━━━━━━━━━━┿━━━━━━━━━━━━━━━━━━━━━━━")')
        write(60,'("lmax-error in spectral space")')
        write(60,'("                 ┃             lmax-error           │             l2 - error           │      location of max    ")')
        write(60,'("   m   n it. im. ┃     zspsc2    zspsc2b    zspsc3a │     zspsc2    zspsc2b    zspsc3a │ initial index_max  rank")')
        write(60,'("━━━━━━━━━━━━━━━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┿━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┿━━━━━━━━━━━━━━━━━━━━━━━━━━━")')
        call flush(40)
        call flush(60)
      end if
    end if
    
  end subroutine init_check_fields

  !===================================================================================================
  ! Subroutine close_check_fields:
  ! Close the files in which the errors are written.
  !===================================================================================================

  subroutine close_check_fields(lwrite_errors, nsmax, myproc)

    implicit none

    logical, intent(in) :: lwrite_errors
    integer(kind=jpim), intent(in) :: nsmax
    integer(kind=jpim), intent(in) :: myproc

    if (myproc == 1) then
      close(40)
      close(60)
    end if

  end subroutine close_check_fields

  !===================================================================================================
  ! Function check_gp_fields:
  ! Compute the errors in gridpoint space by comparing the analytically computed values in the
  ! *_analytic variables with the values computed by the ectrans library. This function returns the
  ! maximum relative error. A summary of the errors is written to file. For lwrite_errors this function
  ! also writes a file with all the latitudes where the errors exceed the tolerance. An attempt is made
  ! to analyze in which part of the code the errors originate (written to standard output unit nout).
  !===================================================================================================

  function check_gp_fields(rtolerance, lwrite_errors, nflevg, nfld, jstep, nzonal, ntotal, &
    & limag, zreel, zgp2, zgp3a, zgpuv, zsph_analytic, znsde_analytic, zewde_analytic, &
    & zu_analytic, zv_analytic, zuder_analytic, zvder_analytic, nout, nsmax, luse_mpi, &
    & ngptotg, lscders, lvordiv, luvders, myproc, nproc, cgrid, cBinID) result(rlmax_error)

    use parkind1, only: jprb, jprd
    use analytic_solutions_mod, only: gelam, gelat
    use mpl_module

    implicit none

    real(kind=jprd), intent(in) :: rtolerance
    logical, intent(in) :: lwrite_errors
    integer(kind=jpim), intent(in) :: jstep, nzonal, ntotal, nflevg, nfld
    logical, intent(in) :: limag
    real(kind=jprd), intent(in) :: zreel(:,:,:), zgp2(:,:,:), zgp3a(:,:,:,:), zgpuv(:,:,:,:) ! should be jprb
    real(kind=jprd), intent(in) :: zsph_analytic(:,:), znsde_analytic(:,:), zewde_analytic(:,:), zu_analytic(:,:), zv_analytic(:,:), zuder_analytic(:,:), zvder_analytic(:,:)
    real(kind=jprd) :: rlmax_error
    integer(kind=jpim), intent(in) :: nout, nsmax, ngptotg, myproc, nproc
    character(len=16), intent(in)  :: cgrid
    character(len=50), intent(in) :: cBinID   ! binary ID (normally everything in the name of the binary
        ! after "-analytic-"). This is used to create unique file names for the error files
        ! which include information about usage of CPU/GPU, OpenACC/OpenMP, ...
    logical, intent(in) :: luse_mpi, lscders, lvordiv, luvders
    real(kind=jprd) :: rlmax_quo, rlmax_nsde_quo, rlmax_ewde_quo, rlmax_u_quo, rlmax_v_quo, rlmax_uder_quo, rlmax_vder_quo, rlmax_fac, rlmax_nsde_fac, rlmax_ewde_fac, rlmax_u_fac, rlmax_v_fac, rlmax_uder_fac, rlmax_vder_fac
    real(kind=jprd) :: rlmax_errors(nflevg*nfld+2), rlmax_errors_nsde(nflevg*nfld+2), rlmax_errors_ewde(nflevg*nfld+2), rlmax_errors_uv(2*nflevg), rlmax_errors_uvder(2*nflevg)
    real(kind=jprd) :: rl2_errors(nflevg*nfld+2), rl2_errors_nsde(nflevg*nfld+2), rl2_errors_ewde(nflevg*nfld+2), rl2_errors_uv(2*nflevg), rl2_errors_uvder(2*nflevg)
    logical :: lpassed(nflevg*nfld+2), lpassed_nsde(nflevg*nfld+2), lpassed_ewde(nflevg*nfld+2), lpassed_all
    integer :: i, j, ntests
    character(len=100) :: filename

    ntests = nflevg*nfld+2

    ! --------------------------------------------------------------
    ! Scalar fields zreel(:,1,:), zpg2(:,1,:) and zgp3a(:,:,1:nfld,:):
    ! --------------------------------------------------------------

    ! compute chosen normalization factor:
    rlmax_quo = maxval(abs( zsph_analytic(:,:)))
    rlmax_fac = 1.0_jprd
    if (luse_mpi) then
      call mpl_allreduce(rlmax_quo, 'max', ldreprod=.false.)
    end if
    if(rlmax_quo>0.0) rlmax_fac = 1.0_jprd/ rlmax_quo
    ! compute lmax- and l2-errors:
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

    ! --------------------------------------------------------------
    ! Scalar derivatives zreel(:,2:3,:), zpg2(:,2:3,:) and zgp3a(:,:,nfld:2*nfld,:):
    ! --------------------------------------------------------------

    if (lscders) then
      ! compute chosen normalization factor:
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
    ! compute lmax- and l2-errors of north-south-derivatives:
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
    ! compute lmax- and l2-errors of east-west-derivatives:
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

    ! --------------------------------------------------------------
    ! Wind speed components zgpuv:
    ! --------------------------------------------------------------

    if (lvordiv) then
      ! compute chosen normalization factor:
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
    ! compute lmax- and l2-errors of u and v:
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

    ! --------------------------------------------------------------
    ! East-west derivatives of wind speed components zgpuv:
    ! --------------------------------------------------------------

    if (luvders) then
      ! compute chosen normalization factor:
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
      ! compute lmax- and l2-errors of u and v east-west derivatives:
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

    ! --------------------------------------------------------------
    ! Write summary of all errors:
    ! --------------------------------------------------------------
   
    if(lwrite_errors .and. myproc == 1) then
      write(40,'(3i4,L4" ┃",3e11.3," │",3e11.3," │",3e11.3," │",2e11.3," │",2e11.3)') nzonal,ntotal,jstep,limag, &
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
      call flush(40)
    end if

    ! --------------------------------------------------------------
    ! Write detailed errors where they exceed the tolerance:
    ! --------------------------------------------------------------
   
    do i=1,ntests
      lpassed(i) = (rlmax_errors(i) < rtolerance)
      lpassed_nsde(i) = (rlmax_errors_nsde(i) < rtolerance)
      lpassed_ewde(i) = (rlmax_errors_ewde(i) < rtolerance)
    end do

    lpassed_all = rlmax_error < rtolerance

    if(.not. lpassed_all) then
      if(lwrite_errors) then
        write(filename,'(3a,i0,3a,i0,a,i0,a,i0,a,i0,a)')'test-failed-',trim(cBinID),'-gridpoint_T',nsmax,'_',trim(cgrid),'_proc',myproc,'_',nproc,'_n',ntotal,'_m',nzonal,'.txt'
        open(50, file = filename, status='replace')
        write(50,'("myproc=",i4," m=",i4," n=",i4," jstep=",i4," imag=",L1)')myproc,nzonal,ntotal,jstep,limag
        write(50,'(a,a)')"                     ┃                                                        │                 north-south derivatives                │", &
        & "                 east-west derivatives                  │       wind speed component u     │       wind speed component v     │    east-west derivative of u     │    east-west derivative of v     "
        write(50,'(a,a)')"     lat/°     lon/° ┃ analytical        pgp       pgp2      pgp3a  max-error │ analytical        pgp       pgp2      pgp3a  max-error │", &
        & " analytical        pgp       pgp2      pgp3a  max-error │ analytical      zgpuv  max-error │ analytical      zgpuv  max-error │ analytical      zgpuv  max-error │ analytical      zgpuv  max-error "
        write(50,'(a,a)')"━━━━━━━━━━━━━━━━━━━━━╋━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┿━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┿", &
        & "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┿━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┿━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┿━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┿━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
        do j=1,ubound(gelat,2)
          do i=1,ubound(gelat,1)
            if(max(abs(zreel(i,1,j)-zsph_analytic(i,j))*rlmax_fac,abs(zgp2(i,1,j)-zsph_analytic(i,j) )*rlmax_fac, &
              & maxval(abs((zgp3a(i,:,1,j)-zsph_analytic(i,j) )))*rlmax_fac, abs(zreel(i,2,j)-znsde_analytic(i,j))*rlmax_nsde_fac, &
              & abs(zgp2(i,2,j)-znsde_analytic(i,j))*rlmax_nsde_fac, maxval(abs((zgp3a(i,:,2,j)-znsde_analytic(i,j))))*rlmax_nsde_fac, &
              & abs(zreel(i,3,j)-zewde_analytic(i,j))*rlmax_ewde_fac, abs(zgp2(i,3,j)-zewde_analytic(i,j))*rlmax_ewde_fac, &
              & maxval(abs((zgp3a(i,:,3,j)-zewde_analytic(i,j))))*rlmax_ewde_fac, maxval(abs(zgpuv(i,:,1,j)-zu_analytic(i,j)))*rlmax_u_fac, &
              & maxval(abs(zgpuv(i,:,2,j)-zv_analytic(i,j)))*rlmax_v_fac) > rtolerance) then
              write(50,'(2i4,2f10.3," ┃",5e15.7," │",5e11.3," │",5e11.3," │",3e11.3," │",3e11.3," │",3e11.3," │",3e11.3)') &
                & i, j, gelat(i,  j)*180/z_pi,gelam(i,j)*180/z_pi, &
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
      ! Attempt to analyze where the errors originate:
      if(count(.not. lpassed)+count(.not. lpassed_nsde)+count(.not. lpassed_ewde) == 3*ntests) then
        write(nout,'("Error: all tests (fields and derivatives) fail for m=",i4," n=",i4)')nzonal,ntotal
      else if((count(.not. lpassed) == ntests) .and. (count(.not. lpassed_nsde)+count(.not. lpassed_ewde) == 0)) then
        write(nout,'("Error: all derivatives are correct but all fields zreel, zgp2 and zgp3a are wrong for m=",i4," n=",i4)')nzonal,ntotal
      else if((count(.not. lpassed) == 0) .and. (count(.not. lpassed_nsde)+count(.not. lpassed_ewde) == 2*ntests)) then
        write(nout,'("Error: all fields zreel, zgp2 and zgp3a are correct but all derivatives are wrong for m=",i4," n=",i4)')nzonal,ntotal
      else if((count(.not. lpassed) == 0) .and. (count(.not. lpassed_nsde) == ntests) .and. (count(.not. lpassed_ewde) == 0)) then
        write(nout,'("Error: all fields and east-west derivatives are correct but all north-south derivatives are wrong for m=",i4," n=",i4)')nzonal,ntotal
      else if(lpassed(1) .and. lpassed_nsde(1) .and. lpassed_ewde(1) .and. (count(.not. lpassed(2:))+count(.not. lpassed_nsde(2:))+count(.not. lpassed_ewde(2:)) == 3*(ntests-1))) then
        write(nout,'("Error: first invtrans is correct but all values are wrong in second invtrans call for m=",i4," n=",i4)')nzonal,ntotal
      else
        if(count(.not. lpassed)>0) write(nout,'("Error: "i3," fields out of ",i4," are wrong for m=",i4," n=",i4)')count(.not. lpassed),ntests,nzonal,ntotal
        if(count(.not. lpassed_nsde)>0) write(nout,'("Error: "i3," north-south derivatives out of ",i4," are wrong for m=",i4," n=",i4)')count(.not. lpassed_nsde),ntests,nzonal,ntotal
        if(count(.not. lpassed_ewde)>0) write(nout,'("Error: "i3," east-west derivatives out of ",i4," are wrong for m=",i4," n=",i4)')count(.not. lpassed_ewde),ntests,nzonal,ntotal
      end if
      call flush(nout)
      stop "error check in grid point space not passed"
    end if
    
  end function check_gp_fields
  
  !===================================================================================================
  ! Function check_sp_fields:
  ! Compute the errors in spectral space by comparing the values computed by the ectrans library with
  ! the known initial values (which is 1 for index nindex and zero everywhere else). This function
  ! returns the maximum relative error. A summary of the errors is written to file.
  !===================================================================================================

  function check_sp_fields(rtolerance, lwrite_errors, nflevg, nfld, jstep, nzonal, ntotal, nindex, limag, zspsc2, &
    & zspsc2b, zspsc3a, nout, luse_mpi, nsmax, myproc, nproc) result(rlmax_error)

    use parkind1, only: jprb
    use mpl_module

    implicit none

    real(kind=jprd), intent(in) :: rtolerance
    logical, intent(in) :: lwrite_errors
    integer(kind=jpim), intent(in) :: jstep, nzonal, ntotal, nflevg, nfld, nindex, nout, nsmax, myproc, nproc
    logical, intent(in) :: limag, luse_mpi
    real(kind=jprd), intent(in) :: zspsc2(:,:), zspsc2b(:,:), zspsc3a(:,:,:) ! should be kind jprb
    real(kind=jprd) :: rlmax_error
    real(kind=jprd) :: rlmax_error_zspsc2, rlmax_error_zspsc2b, rlmax_error_zspsc3a
    real(kind=jprd) :: rl2_error_zspsc2, rl2_error_zspsc2b, rl2_error_zspsc3a
    real(kind=jprd) :: rlmax_error_zspsc2_local
    integer(kind=jpim) :: index_max, rank_max, nindex_max, j
    real(kind=jprd), allocatable :: zindex(:)

    logical :: lpassed_zspsc2, lpassed_zspsc2b, lpassed_zspsc3a, lpassed_all

    ! all spectral arrays are initialized in the subroutine initialize_spectral_arrays
    ! in ectrans-analytic.F90 with zero except for one index with the number nindex.
    ! The index nindex is -1 if the index does not exist on the current MPI process.
    ! This initial field is constructed here in the array zindex:
    allocate(zindex(size(zspsc2,2)))
    zindex = 0.0
    if(nindex>0) zindex(nindex) = 1.0
    rlmax_error_zspsc2 = 0.0
    ! compute lmax-error of zspsc2 in a loop to also find the index of the maximum:
    do j = 1,size(zspsc2,2)
      if(abs(zspsc2(1,j)-zindex(j))>rlmax_error_zspsc2) then
        rlmax_error_zspsc2 = abs(zspsc2(1,j)-zindex(j))
        index_max = j
      end if
    end do
    ! compute l2-error of zspsc2:
    rl2_error_zspsc2 = sum((zspsc2(1,:)-zindex(:))**2)
    rlmax_error_zspsc2_local = rlmax_error_zspsc2
    ! compute lmax- and l2-errors for zspsc2b and zspsc3a:
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
      ! find MPI rank and index of maximum error for zspsc2:
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
    ! write summary of all errors to file:
    if(lwrite_errors .and. myproc==1) then
      write(60,'(3i4,L4," ┃",3e11.3," |",3e11.3," |",3i8)') nzonal,ntotal,jstep,limag, &
      & rlmax_error_zspsc2,  &
      & rlmax_error_zspsc2b, &
      & rlmax_error_zspsc3a, &
      & rl2_error_zspsc2,    &
      & rl2_error_zspsc2b,   &
      & rl2_error_zspsc3a,   &
      & nindex_max, index_max, rank_max
      call flush(60)
    end if
    if(.not. lpassed_all) then
      if((.not. lpassed_zspsc2) .and. (.not. lpassed_zspsc2b) .and. (.not. lpassed_zspsc3a)) write(nout,'("Error: all spectral fields are wrong for m=",i4," n=",i4)')nzonal,ntotal
      write(nout,'("there must be something wrong in the direct transform")')
      call flush(nout)
      stop "error check in spectral space not passed"
    end if
  end function check_sp_fields

  !===================================================================================================
  
end module check_results_mod