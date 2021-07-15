module CalcEfield
  use SimulationParameters,   only: NX, FIELD_ABSORB_BOUNDARY, HAS_ELECTROSTATIC_COMPONENT
  use SupplementalParameters, only: tcs, damping_fm_rd, damping_fm_rr
  implicit none
  public efield, efieldBoundary
contains
  subroutine efield(by,bz,jx,jy,jz,ex,ey,ez)
    double precision,intent(in),dimension(:)    :: by, bz
    double precision,intent(in),dimension(:)    :: jx, jy, jz
    double precision,intent(inout),dimension(:) :: ex, ey, ez

    integer :: ix

    ! Solve Maxwell equation \nabla \cdot B = \mu_0 J + 1/(c^2) * \frac{\partial E}{\partial t}

    do ix = 2, NX+1
       ex(ix) = ex(ix) - jx(ix)*2.0d0
    end do

    do ix = 2, NX+1
       ey(ix) = ey(ix) - ((bz(ix)   - bz(ix-1)) * tcs + jy(ix)*2.0d0) * damping_fm_rr(ix)
    end do

    do ix = 2, NX+1
       ez(ix) = ez(ix) + ((by(ix+1) - by(ix  )) * tcs - jz(ix)*2.0d0) * damping_fm_rr(ix)
    end do

    ! Absorption boundary - Electromagnetic component only
    if ( FIELD_ABSORB_BOUNDARY == .true.) then
       do ix = 2, NX+1
          ey(ix) = ey(ix) * damping_fm_rd(ix)
          ez(ix) = ez(ix) * damping_fm_rd(ix)
       end do
    endif

    ! Calculate boundary condition
    ex(1) = ex(NX+1)
    ex(NX+2) = ex(2)
    call efieldBoundary(ey,ez)
    
    ! Reset ex 
    if (HAS_ELECTROSTATIC_COMPONENT == .false.) ex = 0d0
  end subroutine efield

  subroutine efieldBoundary(ey,ez)
    double precision,intent(inout),dimension(:) :: ey, ez

    if ( FIELD_ABSORB_BOUNDARY == .true.) then
      ey(NX+1:NX+2) = 0d0
      ez(NX+1:NX+2) = 0d0
      ey(1:2) = 0d0
      ez(1:2) = 0d0
    else
      ey(NX+2) = ey(2)
      ez(NX+2) = ez(2)
      ey(1)    = ey(NX+1)
      ez(1)    = ez(NX+1)
    endif
  end subroutine efieldBoundary
end module CalcEfield
