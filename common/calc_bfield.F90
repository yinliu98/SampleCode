module CalcBfield
  use SimulationParameters,   only: NX, FIELD_ABSORB_BOUNDARY
  use SupplementalParameters, only: damping_fm_rd, damping_fm_rr
  use InitFieldParticle,      only: bx0, by0
  implicit none
  private
  public bfield, bfieldBoundary
contains
  subroutine bfield(ey,ez,by,bz)
    double precision,intent(in),dimension(:)    :: ey, ez
    double precision,intent(inout),dimension(:) :: by, bz

    integer :: ix

    ! Solve Maxwell equation \nabla \cdot E = -\frac{\partial B}{\partial t}

    do ix = 2, NX+1
       by(ix) = by(ix) + (ez(ix) - ez(ix-1)) * damping_fm_rr(ix)
    end do

    do ix = 2, NX+1
       bz(ix) = bz(ix) - (ey(ix+1) - ey(ix)) * damping_fm_rr(ix)
    end do

    ! Absorption boundary
    if ( FIELD_ABSORB_BOUNDARY == .true.) then
       do ix = 2, NX+1
          by(ix) = (by(ix)-by0) * damping_fm_rd(ix) + by0
          bz(ix) = bz(ix) * damping_fm_rd(ix)
       end do
    endif

    ! Calculate boundary condition
    call bfieldBoundary(by,bz)
  end subroutine bfield

  subroutine bfieldBoundary(by,bz)
    double precision,intent(inout),dimension(:) :: by, bz

    if ( FIELD_ABSORB_BOUNDARY == .true.) then
      by(NX+1:NX+2) = 0d0
      bz(NX+1:NX+2) = 0d0
      by(1:2) = 0d0
      bz(1:2) = 0d0
    else
      by(NX+2) = by(2)
      bz(NX+2) = bz(2)
      by(1)    = by(NX+1)
      bz(1)    = bz(NX+1)
    endif
  end subroutine bfieldBoundary
end module CalcBfield
