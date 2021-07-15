module CalcPhaseBunch
  use SimulationParameters, only: NTIME, NX, cv, NS, NP,  &
                                  WPIA_THETA_H => WPIA_NUM_ZETA, &
                                  WPIA_V_H => WPIA_NUM_V_PARA, &
                                  NG => iX_ZETA, &
                                  X_ZETA_MIN, X_ZETA_MAX, &
                                  omega_p, omega_c, EXT_MIRROR, ext_mirror_a, ext_omega_b0, &
                                  ISKIP, iXSKIP
  use SupplementalParameters, only: q, XLEN, csq, IT_INTEGRAL_STEPS
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

  integer,parameter :: NNX = min(NX, 1024)
  integer,parameter :: VY_H = WPIA_THETA_H+1

  double precision, save, dimension(NNX+8, VY_H, NG, NS, 4):: hist
  !DEC$ATTRIBUTES ALIGN: 64:: hist
 
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! Access permission
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
  private
  public x_zeta_distribution_fb
contains
  subroutine x_zeta_distribution_fb(it,x,vx,vy,vz,iq)
    integer         ,intent(in)                 :: it
    double precision,dimension(:),intent(in)    :: x
    double precision,dimension(:),intent(in)    :: vx, vy, vz
    integer(1)      ,dimension(:),intent(in)    :: iq
    
    integer :: iw, is, ig, ix, iz, itgroup, istep

#ifndef _DEBUG
    integer,parameter :: histogram_size = (NNX+8) * VY_H * NG * NS * 4
#endif

    if (mod(it, ISKIP) == ISKIP - IT_INTEGRAL_STEPS/2 + 1 .or. it == 1) then
      istep = 1
      hist = 0d0
    else 
      istep = istep + 1
    end if

    call x_zeta_distribution(1,x,vx,vy,vz,iq,by_fwd,bz_fwd)
    call x_zeta_distribution(2,x,vx,vy,vz,iq,by_bwd,bz_bwd)

    if (mod(it, ISKIP) == IT_INTEGRAL_STEPS/2) then

#ifndef _DEBUG
      if (irank_mpi == 0) then
        call MPI_REDUCE(MPI_IN_PLACE,hist       ,histogram_size,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr_mpi)
      else
        call MPI_REDUCE(hist       ,hist       ,histogram_size,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr_mpi)
      endif
#endif
      ! Boundary condition
      do iw = 1, 4
        do is = 1, NS
          do ig = 1, NG
            do iz = 1, VY_H-1
              hist(2 ,iz,ig,is,iw) = hist(2 ,iz,ig,is,iw) + hist(NNX+2,iz,ig,is,iw)
            end do

            do ix = 2, NNX+1
              hist(ix,1,ig,is,iw) = hist(ix,1,ig,is,iw) + hist(ix,VY_H,ig,is,iw)
            end do
            hist(2,1,ig,is,iw) = hist(2,1,ig,is,iw) + hist(NNX+2, VY_H,ig,is,iw)
          end do
        end do
      end do

#ifndef _DEBUG
      call MPI_BARRIER(MPI_COMM_WORLD,ierr_mpi)
#endif

      itgroup = it - IT_INTEGRAL_STEPS/2 + 1
      do iw = 1, 4
        do is = 1, NS
          do ig = 1, NG
            hist(:,:,ig,is,iw) = hist(:,:,ig,is,iw)/sum(hist(:,:,ig,is,iw))
          end do
        end do
      end do

      do iw = 1, 2
        do is = 1, NS
          do ig = 1, NG 
            call outputSpaceDiag(itgroup,is,ig,iw,hist(2:NNX+1,1:VY_H-1,ig,is,iw)/dble(istep),'x_zeta')
            call outputSpaceDiag(itgroup,is,ig,iw,hist(2:NNX+1,1:VY_H-1,ig,is,iw+2)/dble(istep),'x_zeta_vr')
          end do
        end do
      end do
    end if
   end subroutine

  subroutine x_zeta_distribution(iw,x,vx,vy,vz,iq,by,bz)
    integer,                      intent(in)  :: iw
    double precision,dimension(:),intent(in)  :: x
    double precision,dimension(:),intent(in)  :: vx, vy, vz
    integer(1)      ,dimension(:),intent(in)  :: iq
    double precision,dimension(:),intent(in)  :: by,bz

    double precision :: zeta_norm
    double precision :: v_abs, b_abs, v_dot_b, v_cross_b_x
    double precision :: cos_vb
    double precision :: zeta_rad
    double precision :: zetap

    integer :: is, ig, ip, ix, izetap
    integer :: ifirst, ilast

    !Relocated field half grid for particles
    double precision,dimension(NX+8) :: by_avg
    !DEC$ATTRIBUTES ALIGN: 64:: by_avg

    ! Position in Field grid
    double precision :: x_p
    ! Indexs surrounding particle
    integer :: ixp
    ! Width : x_field_grid - ixp
    double precision :: xw, zw
    !Weight for integration
    double precision :: sf1, sf2, sf3, sf4
    double precision :: p_by, p_bz

    ! wave direction for phase direction: Forward waves: -1d0, Backward waves: 1d0
    double precision :: wave_phase_sign
    double precision,parameter :: XLENC = XLEN/2d0
    double precision :: x_eq, v_phase, gamma, vr, tmp, omega_c_l, v_diff
    double precision :: vr_term_a, vr_term_b, vr_term_c, vr_zero, v_para
    logical :: is_skip_ip

    double precision,parameter :: X_NORM = dble(NNX)/dble(NX)

    ! Forward waves: -1d0, Backward waves: 1d0
    wave_phase_sign = -2d0 * dble(iw) + 3d0

    !==============================================================================================
    ! Field cancelation to prevent from self-force oscillation
    !==============================================================================================
    by_avg = 0.0d0
    do ix = 2, NX+1
       by_avg(ix) = (by(ix+1) + by(ix))*0.5d0
    end do
    by_avg(1) = by_avg(NX+1)

    zeta_norm = dble(WPIA_THETA_H) / (2d0*PI)
    is_skip_ip = .false.
    
    !$omp parallel private(is, ifirst, ilast) reduction(+:hist)
    ifirst = 0
    ilast = 0
    do is = 1, NS
      ifirst = ilast + 1
      ilast  = ilast + NP(is)

      !$omp do private(ip,ig,is_skip_ip, x_p, ixp, sf1, sf2), &
      !$omp private(v_abs, b_abs, v_dot_b, v_cross_b_x, cos_vb, zeta_rad), &
      !$omp private(zetap, izetap, xw, zw, sf3, sf4), &
      !$omp private(p_by, p_bz), &
      !$omp private(x_eq, omega_c_l, tmp, v_phase, gamma, vr, v_diff), &
      !$omp private(vr_term_a, vr_term_b, vr_term_c, vr_zero, v_para)
      do ip = ifirst, ilast
        if ( iq(ip) == 0 ) cycle
        do ig = 1, NG
          if( vx(ip) >= X_ZETA_MIN(ig) * cv .and. vx(ip) < X_ZETA_MAX(ig) *cv ) exit
          if(ig == NG) is_skip_ip = .true. ! the velocity of particle is out of range 
        end do         
        if(is_skip_ip) then 
          is_skip_ip = .false.
          cycle ! the velocity of particle is out of range
        endif

        x_p = x(ip) + 1.5d0
        ixp = floor(x_p)
            
        sf2 = x_p - dble(ixp) ! 0 <= xw <= 1.0
        sf1 = 1 - sf2

        ! Obtain electromagnetic field at the particle position by linear interporation
        p_by=sf1*by_avg(ixp) + sf2*by_avg(ixp+1)
        p_bz=sf1*bz(ixp)     + sf2*bz(ixp+1)

        ! Angle between v and B is calculated from cos theta = v*b/(|v||b|)
        ! Rotation derection is determined from v^b
        v_abs = max(dsqrt(vy(ip)**2 + vz(ip) **2),1d-30)
        b_abs = max(dsqrt(p_by**2 + p_bz **2)      ,1d-30)

        v_dot_b = vy(ip)*p_by + vz(ip)*p_bz
        v_cross_b_x = vy(ip) * p_bz - vz(ip) * p_by

        cos_vb = v_dot_b/(v_abs * b_abs)
        cos_vb = max(cos_vb, -1d0)
        cos_vb = min(cos_vb,  1d0)
        zeta_rad = dacos(cos_vb) * wave_phase_sign * dsign(1d0, q(is) * v_cross_b_x)
         ! Domain of definition is [0 2PI]
        zetap = zeta_rad * zeta_norm
        zetap = dmod(zetap+dble(WPIA_THETA_H), dble(WPIA_THETA_H)) + 1d0
        izetap = floor(zetap)
        zw = zetap  - dble(izetap)

        x_p = x(ip)*X_NORM+1.5d0
        ixp = floor(x_p)
            
        sf2 = x_p - dble(ixp) ! 0 <= xw <= 1.0
        sf1 = 1 - sf2

        xw = sf2
        sf3 = xw*zw
        sf2 = xw-sf3
        sf4 = zw-sf3
        sf1 = 1.0d0-xw-sf4

        hist(ixp  ,izetap  ,ig,is,iw) = hist(ixp  ,izetap  ,ig,is,iw) + sf1
        hist(ixp+1,izetap  ,ig,is,iw) = hist(ixp+1,izetap  ,ig,is,iw) + sf2
        hist(ixp+1,izetap+1,ig,is,iw) = hist(ixp+1,izetap+1,ig,is,iw) + sf3
        hist(ixp  ,izetap+1,ig,is,iw) = hist(ixp  ,izetap+1,ig,is,iw) + sf4

        x_eq = x(ip) - XLENC
        if (EXT_MIRROR == .true. ) then
          omega_c_l = dabs(omega_c) * (1d0 + ext_mirror_a * x_eq ** 2)
        else
          omega_c_l = dabs(omega_c)
        endif

        tmp = ext_omega_b0 *(omega_c_l - ext_omega_b0)
        v_phase = cv * dsqrt(tmp/(tmp + omega_p(1)**2))
        gamma = cv/dsqrt(csq-vx(ip)**2-vy(ip)**2-vz(ip)**2)
        vr = dabs(v_phase * (1d0 - omega_c_l/(gamma*ext_omega_b0)))
        vr = ((iw - 1) * 2d0 - 1d0) * vr ! Forward waves -> -, Backward waves -> +
        v_diff = vx(ip) - vr

        vr_term_a = csq + v_phase**2 * omega_c_l**2 /  ext_omega_b0**2
        vr_term_b = -cv * v_phase
        vr_term_c = v_phase ** 2*(1d0 - omega_c_l**2/ext_omega_b0**2)
        vr_zero = (-vr_term_b - sqrt(vr_term_b**2 - vr_term_a*vr_term_c)) / vr_term_a * cv
        vr_zero = ((1 - iw) * 2d0 + 1d0) * vr_zero ! Forward waves -> v_offset*1d0, Backward waves -> voffset*-1d0

        v_para = vr_zero + v_diff
        v_para = min(v_para, cv)
        v_para = max(v_para, -cv)

        do ig = 1, NG
          if( v_para >= X_ZETA_MIN(ig) * cv .and. v_para < X_ZETA_MAX(ig) *cv ) exit
          if(ig == NG) is_skip_ip = .true. ! the velocity of particle is out of range 
        end do         
        if(is_skip_ip) cycle ! the velocity of particle is out of range

        hist(ixp  ,izetap  ,ig,is,iw+2) = hist(ixp  ,izetap  ,ig,is,iw+2) + sf1
        hist(ixp+1,izetap  ,ig,is,iw+2) = hist(ixp+1,izetap  ,ig,is,iw+2) + sf2
        hist(ixp+1,izetap+1,ig,is,iw+2) = hist(ixp+1,izetap+1,ig,is,iw+2) + sf3
        hist(ixp  ,izetap+1,ig,is,iw+2) = hist(ixp  ,izetap+1,ig,is,iw+2) + sf4
       end do
       !$omp end do nowait
   end do
   !$omp end parallel
 end subroutine x_zeta_distribution
end module CalcPhaseBunch
