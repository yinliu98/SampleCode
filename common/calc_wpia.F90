module CalcWPIA
  use SimulationParameters, only: NX, cv, NS, NP,  &
                                  NG => iV_HISTO, &
                                  NUM_V => HISTOGRAM_NUM_V, &
                                  X_POS => HISTO_X_POS, &
                                  HXW => HISTO_XW, &
                                  WPIA_THETA_H => WPIA_NUM_ZETA, &
                                  WPIA_V_H => WPIA_NUM_V_PARA, &
                                  omega_p, omega_c, EXT_MIRROR, ext_mirror_a, ext_omega_b0 ,&
                                  ISKIP
  use SupplementalParameters, only: q, XLEN, csq, IT_INTEGRAL_STEPS
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

  integer,parameter :: VX_H = WPIA_THETA_H+1
  integer,parameter :: VY_H = 2*WPIA_V_H+1

  double precision, save, dimension(VX_H, VY_H, NG, NS, 4) :: hist
  !DEC$ATTRIBUTES ALIGN: 64:: hist
  integer(8), save, dimension(NG, NS, 2) :: npg
  
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! Access permission
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
  private
  public wave_particle_interaction_analyzer_fb
contains
  subroutine wave_particle_interaction_analyzer_fb(it,x,vx,vy,vz,iq)
    integer         ,intent(in)                 :: it
    double precision,dimension(:),intent(in)    :: x
    double precision,dimension(:),intent(in)    :: vx, vy, vz
    integer(1)      ,dimension(:),intent(in)    :: iq

    integer,save :: istep
    integer :: itgroup
    integer :: is , ig, iw

#ifndef _DEBUG
    integer,parameter :: histogram_size = VX_H*VY_H*NG*NS*4
#endif

    if (mod(it, ISKIP) == ISKIP - IT_INTEGRAL_STEPS/2 + 1 .or. it == 1) then
      istep = 1
      hist = 0d0
      npg = 0
    else
      istep = istep + 1
    end if

    call wave_particle_interaction_analyzer(it,1,x,vx,vy,vz,iq,by_fwd,bz_fwd)
    call wave_particle_interaction_analyzer(it,2,x,vx,vy,vz,iq,by_bwd,bz_fwd)

    if (mod(it, ISKIP) == IT_INTEGRAL_STEPS/2) then
#ifndef _DEBUG
      if (irank_mpi == 0) then
        call MPI_REDUCE(MPI_IN_PLACE,hist,histogram_size,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr_mpi)
        call MPI_REDUCE(MPI_IN_PLACE,npg,NG*NS*2        ,MPI_INTEGER8        ,MPI_SUM,0,MPI_COMM_WORLD,ierr_mpi)
      else
        call MPI_REDUCE(hist,hist,histogram_size,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr_mpi)
        call MPI_REDUCE(npg ,npg ,NG*NS*2       ,MPI_INTEGER8        ,MPI_SUM,0,MPI_COMM_WORLD,ierr_mpi)
      endif
#endif

      ! Boundary condition
      hist(1,:,:,:,:) = hist(1,:,:,:,:) + hist(VX_H,:,:,:,:)
      itgroup = it - IT_INTEGRAL_STEPS/2 + 1

      do is = 1, NS
        do ig = 1, NG
          do iw = 1, 2
            call outputBoxDiag(itgroup,is,ig,iw,hist(1:VX_H-1,1:VY_H, ig, is, iw)/dble(istep)/dble(npg(ig,is,iw)), 'zeta_v_para')
            call outputBoxDiag(itgroup,is,ig,iw,hist(1:VX_H-1,1:VY_H, ig, is, iw+2)/dble(istep)/dble(npg(ig,is,iw)), 'zeta_v_para_vr')
          end do
        end do
      end do    
    end if
  end subroutine

  subroutine wave_particle_interaction_analyzer(it,iw,x,vx,vy,vz,iq,by,bz)
    integer         ,intent(in)               :: it, iw
    double precision,dimension(:),intent(in)  :: x
    double precision,dimension(:),intent(in)  :: vx, vy, vz
    integer(1)      ,dimension(:),intent(in)  :: iq
    double precision,dimension(:),intent(in)  :: by, bz

    double precision :: v_norm, zeta_norm
    double precision :: v_abs, b_abs, v_dot_b, v_cross_b_x
    double precision :: cos_vb
    double precision :: zeta_rad
    double precision :: xmin, xmax
    double precision :: x_eq, v_phase, gamma, vr, tmp, omega_c_l, v_diff
    double precision :: vr_term_a, vr_term_b, vr_term_c, vr_zero

    integer :: is, ip, ig, ix
    integer :: ifirst, ilast
    !Relocated field half grid for particles
    double precision,dimension(NX+8) :: by_avg
    !DEC$ATTRIBUTES ALIGN: 64:: by_avg

    ! Position in Field grid
    double precision :: x_p
    ! Indexs surrounding particle
    integer :: ixp
    double precision :: sf1, sf2

    double precision :: ivpp, izetap
    double precision :: v_para, zetap
    double precision :: xw, yw
    double precision :: sf3, sf4

    double precision :: p_by, p_bz
    ! wave direction for phase direction: Forward waves: -1d0, Backward waves: 1d0
    double precision :: wave_phase_sign
    double precision,parameter :: XLENC = XLEN/2d0
    double precision :: cv_threshold 
        
    cv_threshold =  cv * dble(WPIA_V_H -1d0 - 1d-30)/dble(WPIA_V_H)
    ! Forward waves:  1d0, Backward waves: -1d0
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
    v_norm = dble(WPIA_V_H)/cv

    !$omp parallel private(is, ifirst, ilast) reduction(+:hist) reduction(+:npg)
    ifirst = 0
    ilast = 0
    do is = 1, NS
      ifirst = ilast + 1
      ilast  = ilast + NP(is)

      !$omp do private(ig, xmin, xmax, x_p, ixp, sf1, sf2), &
      !$omp private(v_abs, b_abs, v_dot_b, v_cross_b_x, cos_vb, zeta_rad), &
      !$omp private(zetap, v_para, izetap, ivpp, xw, yw, sf3, sf4), &
      !$omp private(x_eq,omega_c_l, tmp, v_phase, gamma, vr, v_diff), &
      !$omp private(vr_term_a, vr_term_b, vr_term_c, vr_zero), &
      !$omp private(p_by, p_bz)
      do ip = ifirst, ilast
        if ( iq(ip) == 0 ) cycle
        do ig = 1, NG
          xmin = dble(NX/2) + X_POS(ig) - HXW(ig)/2d0
          xmax = dble(NX/2) + X_POS(ig) + HXW(ig)/2d0
          if( x(ip) < xmin .or. x(ip) >= xmax ) cycle

          x_p = x(ip) + 1.5d0
          ixp = floor(x_p)

          sf2 = x_p - dble(ixp) ! 0 <= xw <= 1.0
          sf1 = 1 - sf2

          p_by=sf1*by_avg(ixp) + sf2*by_avg(ixp+1)
          p_bz=sf1*bz(ixp)     + sf2*bz(ixp+1)
          ! Angle between v and B is calculated from cos theta = v*b/(|v||b|)
          ! Rotation derection is determined from v^b
          v_abs = max(dsqrt(vy(ip)**2 + vz(ip) **2),1d-30)
          b_abs = max(dsqrt(p_by**2 + p_bz **2)    ,1d-30)

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

          v_para = vx(ip)  + cv
          v_para  = v_para * v_norm + 1d0
          ivpp = floor(v_para)

          xw = zetap  - dble(izetap)
          yw = v_para - dble(ivpp)

          sf3 = xw*yw
          sf2 = xw-sf3
          sf4 = yw-sf3
          sf1 = 1.0d0-xw-sf4

          hist(izetap  , ivpp  ,ig,is,iw) = hist(izetap  , ivpp  ,ig,is,iw) + sf1
          hist(izetap+1, ivpp  ,ig,is,iw) = hist(izetap+1, ivpp  ,ig,is,iw) + sf2
          hist(izetap+1, ivpp+1,ig,is,iw) = hist(izetap+1, ivpp+1,ig,is,iw) + sf3
          hist(izetap  , ivpp+1,ig,is,iw) = hist(izetap  , ivpp+1,ig,is,iw) + sf4
          npg(ig,is,iw) = npg(ig,is,iw) + iq(ip)

          ! For VR correction
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
          v_para = min(v_para, cv_threshold)
          v_para = max(v_para, -cv_threshold)

          v_para = v_para + cv
          v_para = v_para * v_norm + 1d0
          ivpp = floor(v_para)
          yw = v_para - dble(ivpp)

          sf3 = xw*yw
          sf2 = xw-sf3
          sf4 = yw-sf3
          sf1 = 1.0d0-xw-sf4

          hist(izetap  , ivpp  ,ig,is,iw+2) = hist(izetap  , ivpp  ,ig,is,iw+2) + sf1
          hist(izetap+1, ivpp  ,ig,is,iw+2) = hist(izetap+1, ivpp  ,ig,is,iw+2) + sf2
          hist(izetap+1, ivpp+1,ig,is,iw+2) = hist(izetap+1, ivpp+1,ig,is,iw+2) + sf3
          hist(izetap  , ivpp+1,ig,is,iw+2) = hist(izetap  , ivpp+1,ig,is,iw+2) + sf4
          exit
        end do
      end do
      !$omp end do nowait
    end do
    !$omp end parallel
  end subroutine wave_particle_interaction_analyzer
end module calcWPIA
