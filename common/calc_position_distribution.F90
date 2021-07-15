module CalcPositionDistribution
  use SimulationParameters, only: NX, cv, NS, NP, omega_p, omega_c,  &
                                  EXT_MIRROR, ext_mirror_a, ext_omega_b0, &
                                  ISKIP
  use SupplementalParameters, only: q, mass, XLEN, csq, IT_INTEGRAL_STEPS
  use InitFieldParticle, only: bx0
  use Slps, only: PI
  use CalcForwardBackwardWaves, only:  ey_fwd, ez_fwd, ey_bwd, ez_bwd, &
                                       by_fwd, bz_fwd, by_bwd, bz_bwd
  use OutputHDF5, only: outputSpaceDiag
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

  double precision, dimension(NX, NS, 2) :: j_ep_array, j_bp_array
  !DEC$ATTRIBUTES ALIGN: 64:: j_ep_array, j_bp_array

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! Access permission
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
  private
  public fdistribution, j_r_distribution_fb, anisotropy_distribution
contains
  subroutine fdistribution(it,x,vx,vy,vz,iq)
    integer         ,intent(in)                 :: it
    double precision,dimension(:),intent(in)    :: x
    double precision,dimension(:),intent(in)    :: vx, vy, vz
    integer(1)      ,dimension(:),intent(in)    :: iq

#ifndef _DEBUG
    integer :: f_size
#endif

    double precision,dimension(NX+8) :: fv0, fv1_para, fv1_perp, fv2_para, fv2_perp, fv2
    !DEC$ATTRIBUTES ALIGN: 64:: fv0, fv1_para, fv1_perp, fv2_para, fv2_perp, fv2

    ! Position in Field grid
    double precision :: x_p
    ! Indexs surrounding particle
    integer :: ixp
    !Weight for integration
    double precision :: sf1, sf2
    integer :: ifirst, ilast

    double precision :: gamma
    double precision :: tmp 
    double precision :: v_para_eq_2, v_perp_eq_2, v_para_eq, v_perp_eq
    double precision :: gamma_para, ek, ek_para, ek_perp, inv_sum_fx

    double precision :: NXC  

    integer :: is, ip, ix

#ifndef _DEBUG
    f_size = NX+8
#endif

    NXC = dble(NX)/2d0

    ifirst = 0
    ilast = 0

    do is = 1, NS
      ifirst = ilast + 1
      ilast  = ilast + NP(is)

      fv0 = 0d0
      fv1_para = 0d0
      fv1_perp = 0d0
      fv2_para = 0d0
      fv2_perp = 0d0
      fv2 = 0d0

      !$omp parallel do private(x_p, ixp, sf1, sf2), &
      !$omp private(gamma, tmp, v_para_eq_2, v_perp_eq_2, v_para_eq, v_perp_eq), &
      !$omp private(gamma_para, ek, ek_para, ek_perp), &
      !$omp reduction(+:fv0), &
      !$omp reduction(+:fv1_para), reduction(+:fv1_perp), &
      !$omp reduction(+:fv2_para), reduction(+:fv2_perp), &
      !$omp reduction(+:fv2)
      do ip = ifirst, ilast
        if ( iq(ip) == 0 ) cycle
        x_p = x(ip) + 2.0d0
        ixp = floor(x_p)
         
        sf2 = x_p - dble(ixp)
        sf1 = 1d0 - sf2

        fv0(ixp  ) = fv0(ixp  ) + sf1
        fv0(ixp+1) = fv0(ixp+1) + sf2

        gamma = cv/dsqrt(cv**2 - vx(ip)**2 - vy(ip)**2 - vz(ip)**2)
        
        v_perp_eq_2 = vy(ip)**2 + vz(ip)**2
        tmp =   ext_mirror_a * (x(ip)-NXC)**2 / (1d0 + ext_mirror_a * (x(ip)-NXC)**2)
        v_para_eq_2 = vx(ip)**2 + v_perp_eq_2 * tmp
        v_para_eq = max(dsqrt(v_para_eq_2), 0d0) * dsign(1d0, vx(ip))
        v_perp_eq = max(dsqrt(v_perp_eq_2), 0d0)

        fv1_para(ixp  ) = fv1_para(ixp  ) + sf1 * gamma * v_para_eq
        fv1_para(ixp+1) = fv1_para(ixp+1) + sf2 * gamma * v_para_eq
        fv1_perp(ixp  ) = fv1_perp(ixp  ) + sf1 * gamma * v_perp_eq
        fv1_perp(ixp+1) = fv1_perp(ixp+1) + sf2 * gamma * v_perp_eq

        gamma_para = cv/dsqrt(cv**2 - v_para_eq_2)
        ek      = (gamma - 1d0) * mass(is)*cv**2
        ek_para = (gamma_para - 1d0) * mass(is)*cv**2
        ek_perp = ek - ek_para

        fv2_para(ixp  ) = fv2_para(ixp  ) + sf1 * ek_para
        fv2_para(ixp+1) = fv2_para(ixp+1) + sf2 * ek_para
        fv2_perp(ixp  ) = fv2_perp(ixp  ) + sf1 * ek_perp
        fv2_perp(ixp+1) = fv2_perp(ixp+1) + sf2 * ek_perp
        fv2(ixp  ) = fv2(ixp  ) + sf1 * ek
        fv2(ixp+1) = fv2(ixp+1) + sf2 * ek
      end do
      !$omp end parallel do

#ifndef _DEBUG
      ! Particle decomposition =======================================================================
      if (irank_mpi == 0) then
        call MPI_REDUCE(MPI_IN_PLACE, fv0,     f_size,MPI_DOUBLE_PRECISION,MPI_SUM, 0, MPI_COMM_WORLD,ierr_mpi)
        call MPI_REDUCE(MPI_IN_PLACE, fv1_para,f_size,MPI_DOUBLE_PRECISION,MPI_SUM, 0, MPI_COMM_WORLD,ierr_mpi)
        call MPI_REDUCE(MPI_IN_PLACE, fv1_perp,f_size,MPI_DOUBLE_PRECISION,MPI_SUM, 0, MPI_COMM_WORLD,ierr_mpi)
        call MPI_REDUCE(MPI_IN_PLACE, fv2_para,f_size,MPI_DOUBLE_PRECISION,MPI_SUM, 0, MPI_COMM_WORLD,ierr_mpi)
        call MPI_REDUCE(MPI_IN_PLACE, fv2_perp,f_size,MPI_DOUBLE_PRECISION,MPI_SUM, 0, MPI_COMM_WORLD,ierr_mpi)
        call MPI_REDUCE(MPI_IN_PLACE, fv2,     f_size,MPI_DOUBLE_PRECISION,MPI_SUM, 0, MPI_COMM_WORLD,ierr_mpi)
      else
        call MPI_REDUCE(fv0,         fv0,     f_size,MPI_DOUBLE_PRECISION,MPI_SUM, 0, MPI_COMM_WORLD,ierr_mpi)
        call MPI_REDUCE(fv1_para,    fv1_para,f_size,MPI_DOUBLE_PRECISION,MPI_SUM, 0, MPI_COMM_WORLD,ierr_mpi)
        call MPI_REDUCE(fv1_perp,    fv1_perp,f_size,MPI_DOUBLE_PRECISION,MPI_SUM, 0, MPI_COMM_WORLD,ierr_mpi)
        call MPI_REDUCE(fv2_para,    fv2_para,f_size,MPI_DOUBLE_PRECISION,MPI_SUM, 0, MPI_COMM_WORLD,ierr_mpi)
        call MPI_REDUCE(fv2_perp,    fv2_perp,f_size,MPI_DOUBLE_PRECISION,MPI_SUM, 0, MPI_COMM_WORLD,ierr_mpi)
        call MPI_REDUCE(fv2,         fv2,     f_size,MPI_DOUBLE_PRECISION,MPI_SUM, 0, MPI_COMM_WORLD,ierr_mpi)
      endif
      !==============================================================================================
#endif

      fv0(2) = fv0(2) + fv0(NX+2)
      fv0(1) = fv0(NX+1)
      fv0(NX+2) = fv0(2)

      fv1_para(2) = fv1_para(2) + fv1_para(NX+2)
      fv1_para(1) = fv1_para(NX+1)
      fv1_para(NX+2) = fv1_para(2)

      fv1_perp(2) = fv1_perp(2) + fv1_perp(NX+2)
      fv1_perp(1) = fv1_perp(NX+1)
      fv1_perp(NX+2) = fv1_perp(2)

      fv2_para(2) = fv2_para(2) + fv2_para(NX+2)
      fv2_para(1) = fv2_para(NX+1)
      fv2_para(NX+2) = fv2_para(2)

      fv2_perp(2) = fv2_perp(2) + fv2_perp(NX+2)
      fv2_perp(1) = fv2_perp(NX+1)
      fv2_perp(NX+2) = fv2_perp(2)

      fv2(2) = fv2(2) + fv2(NX+2)
      fv2(1) = fv2(NX+1)
      fv2(NX+2) = fv2(2)

      inv_sum_fx = 1/sum(fv0(2:NX+1))

      ! Normalize
      do ix = 2, NX+2
         fv0(ix) = fv0(ix) * inv_sum_fx
         fv1_para(ix) = fv1_para(ix) * fv0(ix)
         fv1_perp(ix) = fv1_perp(ix) * fv0(ix)
         fv2_para(ix) = fv2_para(ix) * fv0(ix)
         fv2_perp(ix) = fv2_perp(ix) * fv0(ix)
         fv2(ix) = fv2(ix) * fv0(ix)
      end do

      call outputSpaceDiag(it,is,1,1,fv0(2:NX+1),'fv0')
      call outputSpaceDiag(it,is,1,1,fv1_perp(2:NX+1),'fv1_para')
      call outputSpaceDiag(it,is,1,1,fv1_para(2:NX+1),'fv1_perp')
      call outputSpaceDiag(it,is,1,1,fv2_para(2:NX+1),'fv2_para')
      call outputSpaceDiag(it,is,1,1,fv2_perp(2:NX+1),'fv2_perp')
      call outputSpaceDiag(it,is,1,1,fv2(2:NX+1),'fv2')

#ifndef _DEBUG
      call MPI_BARRIER(MPI_COMM_WORLD,ierr_mpi)
#endif
    end do
  end subroutine fdistribution

  subroutine j_r_distribution_fb(it,x,vy,vz,iq)
    integer         ,intent(in)               :: it
    double precision,dimension(:),intent(in)  :: x
    double precision,dimension(:),intent(in)  :: vy, vz
    integer(1)      ,dimension(:),intent(in)  :: iq

    integer :: is, iw, itgroup
    integer,save :: istep


    if (mod(it, ISKIP) == ISKIP - IT_INTEGRAL_STEPS/2 + 1 .or. it == 1) then
      istep = 1
      j_bp_array = 0d0
      j_ep_array = 0d0
    else
      istep = istep + 1
    end if

    call j_r_distribution    (it,1,x,vy,vz,iq,ey_fwd,ez_fwd,by_fwd,bz_fwd)
    call j_r_distribution    (it,2,x,vy,vz,iq,ey_bwd,ez_bwd,by_bwd,bz_bwd)
    
    if (mod(it, ISKIP) == IT_INTEGRAL_STEPS/2) then
      itgroup = it - IT_INTEGRAL_STEPS/2 + 1
      do iw = 1, 2
        do is = 1, NS
          j_ep_array(:,is,iw) = j_ep_array(:,is,iw)/dble(istep)
          j_bp_array(:,is,iw) = j_bp_array(:,is,iw)/dble(istep)

          call outputSpaceDiag(itgroup,is,1,iw,j_bp_array(2:NX+1,is,iw),'jb')
          call outputSpaceDiag(itgroup,is,1,iw,j_ep_array(2:NX+1,is,iw),'je')
        end do
      end do    
    end if
  end subroutine

  subroutine j_r_distribution(it,iw,x,vy,vz,iq,ey,ez,by,bz)
    integer         ,intent(in)               :: it,iw
    double precision,dimension(:),intent(in)  :: x
    double precision,dimension(:),intent(in)  :: vy, vz
    integer(1)      ,dimension(:),intent(in)  :: iq
    double precision,dimension(:),intent(in)  :: ey, ez
    double precision,dimension(:),intent(in)  :: by, bz

#ifndef _DEBUG
    integer :: f_size
#endif

    double precision, dimension(NX+8) :: j_ep, j_bp
    !DEC$ATTRIBUTES ALIGN: 64:: j_ep, j_bp

    integer :: is, ip, ix
    integer :: ifirst, ilast

    !Relocated field half grid for particles
    double precision,dimension(NX+8) :: by_avg
    !DEC$ATTRIBUTES ALIGN: 64:: by_avg

    ! Position in Field grid
    double precision :: x_p
    ! Indexs surrounding particle
    integer :: ixp
    !Weight for integration
    double precision :: sf1, sf2

    double precision :: p_ey, p_ez, p_by, p_bz
    double precision :: tmp_e, tmp_b
    double precision :: XLENC

#ifndef _DEBUG
    f_size = NX+8
#endif
    XLENC = XLEN/2d0

    !==============================================================================================
    ! Field cancelation to prevent from self-force oscillation
    !==============================================================================================
    
    by_avg = 0d0
    do ix = 2, NX+1
       by_avg(ix) = (by(ix+1) + by(ix))*0.5d0
    end do
    by_avg(1) = by_avg(NX+1)

    ifirst = 0
    ilast  = 0

    do is = 1, NS
      ifirst = ilast + 1
      ilast  = ilast + NP(is)

      j_ep = 0.0d0
      j_bp = 0.0d0

      !$omp parallel do private(x_p, ixp, sf1, sf2), &
      !$omp private(p_ey, p_ez, p_by, p_bz), &
      !$omp private(tmp_e, tmp_b), &
      !$omp reduction(+:j_ep), reduction(+:j_bp)
      do ip = ifirst, ilast
        if ( iq(ip) == 0 ) cycle
        x_p = x(ip) + 2.0d0
        ixp = floor(x_p)
        sf2 = x_p - dble(ixp)
        sf1 = 1 - sf2
        p_ey=sf1*ey(ixp)     + sf2*ey(ixp+1)
        p_ez=sf1*ez(ixp)     + sf2*ez(ixp+1)

        x_p = x(ip) + 1.5d0
        ixp = floor(x_p)
        sf2 = x_p - dble(ixp) 
        sf1 = 1d0 - sf2
        p_by=sf1*by_avg(ixp) + sf2*by_avg(ixp+1)
        p_bz=sf1*bz(ixp)     + sf2*bz(ixp+1)

        tmp_e = q(is) * (vy(ip) * p_ey + vz(ip) * p_ez)/dsqrt(p_ey **2 + p_ez ** 2)
        tmp_b = q(is) * (vy(ip) * p_by + vz(ip) * p_bz)/dsqrt(p_by **2 + p_bz ** 2)

        ! Landau resonance
        j_ep(ixp  ) = j_ep(ixp  ) + sf1 * tmp_e
        j_ep(ixp+1) = j_ep(ixp+1) + sf2 * tmp_e

        j_bp(ixp  ) = j_bp(ixp  ) + sf1 * tmp_b
        j_bp(ixp+1) = j_bp(ixp+1) + sf2 * tmp_b
      end do
      !$omp end parallel do

#ifndef _DEBUG
      ! Particle decomposition =======================================================================
      if (irank_mpi == 0) then
        call MPI_REDUCE(MPI_IN_PLACE, j_ep, f_size,MPI_DOUBLE_PRECISION,MPI_SUM,0, MPI_COMM_WORLD,ierr_mpi)
        call MPI_REDUCE(MPI_IN_PLACE, j_bp, f_size,MPI_DOUBLE_PRECISION,MPI_SUM,0, MPI_COMM_WORLD,ierr_mpi)
      else
        call MPI_REDUCE(j_ep,         j_ep, f_size,MPI_DOUBLE_PRECISION,MPI_SUM,0, MPI_COMM_WORLD,ierr_mpi)
        call MPI_REDUCE(j_bp,         j_bp, f_size,MPI_DOUBLE_PRECISION,MPI_SUM,0, MPI_COMM_WORLD,ierr_mpi)
      endif
      !==============================================================================================
#endif

! periodic boundary conditions-----------------------------------------------------------------

      j_ep(2) = j_ep(2) + j_ep(NX+2)
      j_ep(1) = j_ep(NX+1)
      j_ep(NX+2) = j_ep(2)

      j_bp(2) = j_bp(2) + j_bp(NX+2)
      j_bp(1) = j_bp(NX+1)
      j_bp(NX+2) = j_bp(2)

      j_ep_array(:,is,iw) = j_ep_array(:,is,iw) + j_ep(2:NX+1)
      j_bp_array(:,is,iw) = j_bp_array(:,is,iw) + j_bp(2:NX+1)

#ifndef _DEBUG
      call MPI_BARRIER(MPI_COMM_WORLD,ierr_mpi)
#endif
    end do
  end subroutine j_r_distribution

  subroutine anisotropy_distribution(it,x,vx,vy,vz,iq)
    integer         ,intent(in)                 :: it
    double precision,dimension(:),intent(in)    :: x
    double precision,dimension(:),intent(in)    :: vx, vy, vz
    integer(1)      ,dimension(:),intent(in)    :: iq

#ifndef _DEBUG
    integer :: f_size
#endif

    double precision,dimension(NX+8) :: T_para, T_perp, A
    !DEC$ATTRIBUTES ALIGN: 64:: T_para, T_perp, A

    ! Position in Field grid
    double precision :: x_p
    ! Indexs surrounding particle
    integer :: ixp
    !Weight for integration
    double precision :: sf1, sf2
    integer :: ifirst, ilast

    double precision :: p_T_para, p_T_perp, tmp
    double precision :: v_para_eq_2, v_perp_eq_2
    double precision :: NXC

    integer :: is, ip, ix

#ifndef _DEBUG
    f_size = NX+8 
#endif

    NXC = dble(NX)/2
    ifirst = 0
    ilast = 0

    do is = 1, NS
      ifirst = ilast + 1
      ilast  = ilast + NP(is)

      T_para = 0d0
      T_perp = 0d0
      A = 0d0

      !$omp parallel do private(x_p, ixp, sf1, sf2), &
      !$omp private(v_perp_eq_2, tmp, v_para_eq_2, p_T_para, p_T_perp), &
      !$omp reduction(+:T_para), reduction(+:T_perp)
      do ip = ifirst, ilast
        if ( iq(ip) == 0 ) cycle
        x_p = x(ip) + 2.0d0
        ixp = floor(x_p)
 
        sf2 = x_p - dble(ixp)
        sf1 = 1 - sf2

        v_perp_eq_2 = vy(ip)**2 + vz(ip)**2
        tmp =   ext_mirror_a * (x(ip)-NXC)**2 / (1d0 + ext_mirror_a * (x(ip)-NXC)**2)
        v_para_eq_2 = vx(ip)**2 + v_perp_eq_2 * tmp

        p_T_para = v_para_eq_2
        p_T_perp = v_perp_eq_2/2d0

        T_para(ixp  ) = T_para(ixp  ) + sf1 * p_T_para
        T_para(ixp+1) = T_para(ixp+1) + sf2 * p_T_para
        T_perp(ixp  ) = T_perp(ixp  ) + sf1 * p_T_perp
        T_perp(ixp+1) = T_perp(ixp+1) + sf2 * p_T_perp
      end do
      !$omp end parallel do

#ifndef _DEBUG
      ! Particle decomposition =======================================================================
      if (irank_mpi == 0) then
        call MPI_REDUCE(MPI_IN_PLACE, T_para,f_size,MPI_DOUBLE_PRECISION,MPI_SUM, 0, MPI_COMM_WORLD,ierr_mpi)
        call MPI_REDUCE(MPI_IN_PLACE, T_perp,f_size,MPI_DOUBLE_PRECISION,MPI_SUM, 0, MPI_COMM_WORLD,ierr_mpi)
      else
        call MPI_REDUCE(T_para,T_para,f_size,MPI_DOUBLE_PRECISION,MPI_SUM, 0, MPI_COMM_WORLD,ierr_mpi)
        call MPI_REDUCE(T_perp,T_perp,f_size,MPI_DOUBLE_PRECISION,MPI_SUM, 0, MPI_COMM_WORLD,ierr_mpi)
      endif
      !==============================================================================================
#endif

      T_para(2) = T_para(2) + T_para(NX+2)
      T_para(1) = T_para(NX+1)
      T_para(NX+2) = T_para(2)

      T_perp(2) = T_perp(2) + T_perp(NX+2)
      T_perp(1) = T_perp(NX+1)
      T_perp(NX+2) = T_perp(2)

      ! Calculate temperature anisotropy
      do ix = 2, NX+2
         tmp = T_perp(ix)/T_para(ix) - 1d0
         tmp = max(tmp, -1d0)
         tmp = min(tmp, 1d10)
         A(ix) = tmp
      end do

      call outputSpaceDiag(it,is,1,1,A,'A')
#ifndef _DEBUG
      call MPI_BARRIER(MPI_COMM_WORLD,ierr_mpi)
#endif
    end do
  end subroutine anisotropy_distribution
end module CalcPositionDistribution
