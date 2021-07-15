module CalcVelocityDistribution
  use SimulationParameters, only: NTIME, NX, cv, NS, NP,  &
                                  NG => iV_HISTO, &
                                  NUM_V => HISTOGRAM_NUM_V, &
                                  X_POS => HISTO_X_POS, &
                                  HXW => HISTO_XW, &
                                  DAMPING_REGION_ND, &
                                  ISKIP
  use SupplementalParameters, only: q, mass, IT_INTEGRAL_STEPS
  use Slps, only: PI
  use CalcForwardBackwardWaves, only:  ey_fwd, ez_fwd, ey_bwd, ez_bwd, &
                                       by_fwd, bz_fwd, by_bwd, bz_bwd
  use OutputHDF5, only: outputBoxDiag
  !$ use omp_lib
  ! DEBUG mode unables MPI (Compiler option)------------------------------------
  ! -D_DEBUG    => debug mode ON
  ! Not defined => debug mode OFF
  !-----------------------------------------------------------------------------

#ifndef _DEBUG
  use MPI
  use MPIParameters
#endif

  implicit none

#ifdef _DEBUG
  integer,parameter :: irank_mpi = 0
#endif

  integer,parameter :: VX_H = 2*NUM_V + 1
  integer,parameter :: VY_H =   NUM_V + 1

  double precision, save, dimension(VX_H+7, VY_H, NG, NS) :: histogram
  !DEC$ATTRIBUTES ALIGN: 64:: histogram
  integer(8), save, dimension(NG, NS) :: npg

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! Access permission
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
  private
  public vdistribution
contains
  subroutine vdistribution(it,x,vx,vy,vz,iq)
    integer         ,intent(in)                 :: it
    double precision,dimension(:),intent(in)    :: x
    double precision,dimension(:),intent(in)    :: vx, vy, vz
    integer(1)      ,dimension(:),intent(in)    :: iq

#ifndef _DEBUG
    integer,parameter :: histogram_size = (VX_H+7)*VY_H*NG*NS
#endif

    double precision :: v_norm
    double precision :: dv

    double precision :: v_para, v_perp
    double precision :: xmin, xmax

    integer :: is, ip, ig, ivx, ivy
    integer :: ivxp, ivyp
    integer :: ifirst, ilast

    integer :: itgroup
    integer,save :: istep

    if (mod(it, ISKIP) == ISKIP - IT_INTEGRAL_STEPS/2 + 1 .or. it == 1) then
      istep = 1
      histogram = 0d0
      npg = 0
    else
      istep = istep + 1
    end if

    v_norm = dble(NUM_V)/cv
    dv     = cv/dble(NUM_V)

    ifirst = 0
    ilast  = 0

    do is = 1, NS
      ifirst = ilast + 1
      ilast  = ilast + NP(is)

      !$omp parallel do private(ig, xmin, xmax, v_para, v_perp, ivxp, ivyp), &
      !$omp reduction(+:histogram), reduction(+:npg)
      do ip = ifirst, ilast
        if ( iq(ip) == 0 ) cycle

        do ig = 1, NG
          xmin = dble(NX/2) + X_POS(ig) - HXW(ig)/2d0
          xmax = dble(NX/2) + X_POS(ig) + HXW(ig)/2d0 
          if( x(ip) < xmin .or. x(ip) >= xmax ) cycle

          v_para = vx(ip)  + cv
          v_perp = dsqrt(vy(ip)**2 + vz(ip)**2)
          v_perp = max(v_perp, 0d0)

          v_para = v_para * v_norm + 1.5d0
          v_perp = v_perp * v_norm + 1.5d0

          ivxp = floor(v_para)
          ivyp = floor(v_perp)

          histogram(ivxp,ivyp,ig,is) = histogram(ivxp,ivyp,ig,is) + dble(iq(ip))
          npg(ig,is) = npg(ig,is) + iq(ip)
          exit
        end do
      end do
      !$omp end parallel do
    end do

    if (mod(it, ISKIP) == IT_INTEGRAL_STEPS/2) then
      itgroup = it - IT_INTEGRAL_STEPS/2 + 1

#ifndef _DEBUG
      if (irank_mpi == 0) then
        call MPI_REDUCE(MPI_IN_PLACE,histogram,histogram_size,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr_mpi)
        call MPI_REDUCE(MPI_IN_PLACE,npg      ,NG*NS          ,MPI_INTEGER8        ,MPI_SUM,0,MPI_COMM_WORLD,ierr_mpi)
      else
        call MPI_REDUCE(histogram,histogram,histogram_size,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr_mpi)
        call MPI_REDUCE(npg      ,npg      ,NG*NS             ,MPI_INTEGER8        ,MPI_SUM,0,MPI_COMM_WORLD,ierr_mpi)
      endif
#endif

      do ivx = 1, VX_H
        histogram(ivx,1,:,:) = histogram(ivx,1,:,:)/(PI* 0.25d0 * dv**3)
      end do

      do ivy = 2, VY_H
        do ivx = 1, VX_H
          histogram(ivx,ivy,:,:) = histogram(ivx,ivy,:,:)/(2.0d0 * PI * dble(ivy - 1) * dv**3)
        end do
      end do

      do ig = 1, NG
        do is = 1, NS
          histogram(:,:,ig,is) = histogram(:,:,ig,is)/dble(npg(ig,is))/dble(istep)
          call outputBoxDiag(itgroup,is,ig,1,histogram(1:VX_H,:,ig,is),'velocity_histogram')
        end do
      end do    
    end if
  end subroutine vdistribution
end module CalcVelocityDistribution
