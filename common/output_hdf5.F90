module OutputHDF5
  use SimulationParameters, only: NTIME, DX, DT, NX, cv,&
                                  NS, NP, omega_p, QM, v_para, v_perp, v_d, PITCH_DEG,&
                                  omega_c, &
                                  THETA_DEG, &
                                  HAS_ELECTROSTATIC_COMPONENT, &
                                  SUBTRACT_BI_MAXWELLIAN, SUBTRACT_BETA, &
                                  SUBTRACT_RHO, &
                                  IS_PARABOLIC_DISTRIBUTION, &
                                  EXT_MIRROR, MIRROR_NTIME, ext_mirror_a, &
                                  PARTICLE_LOSSCORN_ELIMINATE, &
                                  LOSSCORN_ANGLE_DEG, &
                                  EXT_CURRENT, ctime, cwidth_begin, cwidth_end, &
                                  iANTENNA_OFFSET,  &
                                  ext_j_amp, ext_omega_a, ext_omega_b0,&
                                  iX_WAVE_NUMBER_DIV, &
                                  LOCALIZED_CHARGE_AVG, &
                                  FIELD_ABSORB_BOUNDARY, DAMPING_REGION_ND, &
                                  DAMPING_INTEGRAL_K, &
                                  ISKIP, IXSKIP, &
                                  iV_HISTO, &
                                  HISTO_X_POS, HISTO_XW, &
                                  HISTOGRAM_NUM_V, &
                                  WPIA_NUM_V_PARA, WPIA_NUM_ZETA, &
                                  iX_ZETA, X_ZETA_MIN, X_ZETA_MAX, &
                                  iZETA_V_PARA, TRACKER_MAX_PARTICLE_PER_PROCESS, &
                                  TRACKER_ITIME, TRACKER_SPECIES, TRACKER_POSITION_NUM, &
                                  TRACKER_ZETA_MIN, TRACKER_ZETA_MAX, &
                                  TRACKER_V_PARA_MIN, TRACKER_V_PARA_MAX, &
                                  TRACKER_LABEL, &
                                  VERSION, iSEED
  use SupplementalParameters, only: RN_X, RN_T, RN_V, RN_E_FIELD, RN_POTENTIAL,  &
                                    RN_B_FIELD, RN_CURRENT_DENS, RN_CHARGE_DENS, &
                                    RN_CHARGE, RN_MASS, &
                                    rho0, &
                                    damping_rr, damping_rd, damping_fm_rd, damping_fm_rr
  use CalcEnergy,           only: energy_kinetic, energy_kinetic_para, energy_kinetic_perp, &
                                  energy_e_field, energy_b_field, energy_total
  use CalcMagneticMomentum, only: magnetic_momentum_total
  use CalcForwardBackwardWaves, only:  ey_fwd, ez_fwd, ey_bwd, ez_bwd, &
                                       by_fwd, bz_fwd, by_bwd, bz_bwd
  use DealCommandArguments, only: is_output, file_name_out
  use InitFieldParticle,    only: bx0, by0
  use localTimeToUTC,       only: get_utc_date_and_time
  use h5HLLibrary
  use MPIParameters
  implicit none

  character (:),allocatable,save :: file_name
  character(8)                   :: date
  character(10)                  :: time

  interface outputBoxDiag
    module procedure box1d, box2d, box3d
  end interface outputBoxDiag

  interface outputSpaceDiag
    module procedure Space1d, Space2d, Space3d
  end interface outputSpaceDiag

  private
  public initFile
  public writeRenormalizeFactor, writeDampingFactor
  public outputRawData, output1dField, outputKx
  public outputBoxDiag, outputSpaceDiag, outputParticleTrack
  public WriteEnergyMomentum
  public outputAllData, outputFinalize, outputMirrorInit
  public outputParticleIndex
contains
  !===========================================================================
  ! integer to character for name of group
  !===========================================================================
  character(len=2) function is2char(is)
    integer,intent(in) :: is
    write(is2char,'(i2.2)') is
  end function

  character(len=3) function ig2char(ig)
    integer,intent(in) :: ig
    write(ig2char,'(i3.3)') ig
  end function

  character(len=10) function it2char(it)
    integer,intent(in) :: it
    write(it2char,'(i10.10)') it
  end function
  !===========================================================================
  ! File name
  !===========================================================================
  subroutine outputSetFileName
    if (is_output == .true.) then
      allocate (character(len(file_name_out)) :: file_name)
      file_name = file_name_out
    else
      call get_utc_date_and_time(date=date,time=time)
      allocate (character(18) :: file_name)
      file_name = date(1:8)//'_'//time(1:6)//'.h5'
    endif
  end subroutine outputSetFileName

  !===========================================================================
  ! Making groups and write parameters
  !===========================================================================
  subroutine initFile
    integer :: is
    call outputSetFileName

    if ( irank_mpi /= 0) return

    call initHDF5
    call h5File(file_name, 'x')
    call h5GroupCreate('renormalize_parameters')
    call h5GroupCreate('data')
    call h5GroupCreate('kx')
    call h5GroupCreate('1d_field')
    call h5GroupCreate('total')
    call h5GroupCreate('tracker_index')
    call h5GroupCreate('parameters', enable_grp_mv=.true.)
    call h5GroupCreate('species',   enable_grp_mv=.true.)
    do is =1,NS
      call h5GroupCreate(is2char(is))
    end do
    !===========================================================================
    ! WRITE PARAMETERS
    !===========================================================================
    call h5GroupOpen('parameters', is_abs=.true.)
    !>>>>>>>>>>>>>>>>>>>>>>>>>>s>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! simulation_parameters
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    call h5DatasetCreate('ntime',NTIME)
    call h5DatasetCreate('dx'   ,DX)
    call h5DatasetCreate('dt'   ,DT)
    call h5DatasetCreate('nx'   ,NX)
    call h5DatasetCreate('cv'   ,cv)
    call h5DatasetCreate('ns'   ,NS)

    do is =1,NS
      call h5GroupOpen('parameters/species/'//is2char(is), is_abs=.true.)
      call h5DatasetCreate('np'       ,NP(is))
      call h5DatasetCreate('omega_p'  ,omega_p(is))
      call h5DatasetCreate('qm'       ,QM(is))
      call h5DatasetCreate('v_para'   ,v_para(is))
      call h5DatasetCreate('v_perp'   ,v_perp(is))
      call h5DatasetCreate('v_d'      ,v_d(is))
      call h5DatasetCreate('pitch_deg',PITCH_DEG(is))
      call h5DatasetCreate('subtract_bi_maxwellian',SUBTRACT_BI_MAXWELLIAN(is))
      call h5DatasetCreate('subtract_beta',         SUBTRACT_BETA(is))
      call h5DatasetCreate('subtract_rho',          SUBTRACT_RHO(is))
      call h5DatasetCreate('is_parabolic_distribution',IS_PARABOLIC_DISTRIBUTION(is))
    end do

    call h5GroupOpen('parameters', is_abs=.true.)

    call h5DatasetCreate('omega_c'  ,omega_c  )
    call h5DatasetCreate('theta_deg',THETA_DEG)
    call h5DatasetCreate('ES_mode'  ,HAS_ELECTROSTATIC_COMPONENT)

    call h5GroupCreate('mirror', enable_grp_mv=.true.)
    call h5DatasetCreate('ext_mirror'                 ,EXT_MIRROR                 )
    call h5DatasetCreate('mirro_ntime'                ,MIRROR_NTIME               )
    call h5DatasetCreate('ext_mirror_a'               ,ext_mirror_a               )
    call h5DatasetCreate('particle_losscorn_eliminate',PARTICLE_LOSSCORN_ELIMINATE)
    call h5DatasetCreate('LOSSCORN_ANGLE_DEG'         ,LOSSCORN_ANGLE_DEG         )

    call h5GroupOpen('parameters', is_abs=.true.)
    call h5GroupCreate('antenna', enable_grp_mv=.true.)
    call h5DatasetCreate('EXT_CURRENT'       ,EXT_CURRENT      )
    call h5DatasetCreate('ctime'             ,ctime            )
    call h5DatasetCreate('cwidth_begin'      ,cwidth_begin     )
    call h5DatasetCreate('cwidth_end'        ,cwidth_end       )
    call h5DatasetCreate('iANTENNA_OFFSET'   ,iANTENNA_OFFSET              )
    call h5DatasetCreate('ext_j_amp'         ,ext_j_amp                    )
    call h5DatasetCreate('ext_omega_a'       ,ext_omega_a                  )
    call h5DatasetCreate('ext_omega_b'       ,ext_omega_b0                 )

    call h5GroupOpen('parameters', is_abs=.true.)
    call h5GroupCreate('dispersion', enable_grp_mv=.true.)
    call h5DatasetCreate('iX_WAVE_NUMBER_DIV' ,iX_WAVE_NUMBER_DIV)

    call h5GroupOpen('parameters', is_abs=.true.)
    call h5GroupCreate('poisson', enable_grp_mv=.true.)
    call h5DatasetCreate('LOCALIZED_CHARGE_AVG',LOCALIZED_CHARGE_AVG)

    call h5GroupOpen('parameters', is_abs=.true.)
    call h5GroupCreate('boundary', enable_grp_mv=.true.)
    call h5DatasetCreate('FIELD_ABSORB_BOUNDARY' ,FIELD_ABSORB_BOUNDARY)
    call h5DatasetCreate('DAMPING_REGION_ND'     ,DAMPING_REGION_ND    )
    call h5DatasetCreate('DAMPING_INTEGRAL_K'    ,DAMPING_INTEGRAL_K   )

    call h5GroupOpen('parameters', is_abs=.true.)
    call h5GroupCreate('skip', enable_grp_mv=.true.)
    call h5DatasetCreate('iSKIP'  ,iSKIP)
    call h5DatasetCreate('iXSKIP' ,iXSKIP)

    call h5GroupOpen('parameters', is_abs=.true.)
    call h5GroupCreate('box_diag_info', enable_grp_mv=.true.)
    call h5DatasetCreate('iv_histo'   , iV_HISTO)
    call h5DatasetCreate('histo_x_pos', HISTO_X_POS)
    call h5DatasetCreate('histo_xw'   , HISTO_XW)
    call h5DatasetCreate('histogram_num_v',HISTOGRAM_NUM_V)
    call h5DatasetCreate('wpia_num_v_para',WPIA_NUM_V_PARA)
    call h5DatasetCreate('wpia_num_zeta'  ,WPIA_NUM_ZETA)

    call h5GroupOpen('parameters', is_abs=.true.)
    call h5GroupCreate('space_diag_info', enable_grp_mv=.true.)
    call h5DatasetCreate('ix_zeta'    ,iX_ZETA)
    call h5DatasetCreate('x_zeta_min' ,X_ZETA_MIN)
    call h5DatasetCreate('x_zeta_max' ,X_ZETA_MAX)

    call h5GroupOpen('parameters', is_abs=.true.)
    call h5GroupCreate('tracker_info', enable_grp_mv=.true.)
    call h5DatasetCreate('izeta_v_para' ,iZETA_V_PARA)
    call h5DatasetCreate('tr_max_particle_per_process', TRACKER_MAX_PARTICLE_PER_PROCESS)
    call h5DatasetCreate('tr_itime'     ,TRACKER_ITIME)
    call h5DatasetCreate('tr_species'   ,TRACKER_SPECIES)
    call h5DatasetCreate('tr_pos_num'   ,TRACKER_POSITION_NUM)
    call h5DatasetCreate('tr_zeta_min'  ,TRACKER_ZETA_MIN)
    call h5DatasetCreate('tr_zeta_max'  ,TRACKER_ZETA_MAX)
    call h5DatasetCreate('tr_v_para_min',TRACKER_V_PARA_MIN)
    call h5DatasetCreate('tr_v_para_max',TRACKER_V_PARA_MAX)
    call h5DatasetCreate('tr_label'     ,TRACKER_LABEL)

    call h5GroupOpen('parameters', is_abs=.true.)
    call h5GroupCreate('info', enable_grp_mv=.true.)
    call h5DatasetCreate('version',version)
    call h5DatasetCreate('iseed',iseed)
    call h5DatasetCreate('iprocess',iprocess)
    call h5GroupClose
  end subroutine initFile

  subroutine writeRenormalizeFactor
    if ( irank_mpi /= 0) return
    ! Renormalize factor
    call h5GroupOpen('renormalize_parameters')
    call h5DatasetCreate('rn_x',rn_x)
    call h5DatasetCreate('rn_t',rn_t)
    call h5DatasetCreate('rn_v',rn_v)
    call h5DatasetCreate('rn_e_field',rn_e_field)
    call h5DatasetCreate('rn_potential',rn_potential)
    call h5DatasetCreate('rn_b_field',rn_b_field)
    call h5DatasetCreate('rn_current_dens',rn_current_dens)
    call h5DatasetCreate('rn_charge_dens',rn_charge_dens)
    call h5DatasetCreate('rn_charge',rn_charge)
    call h5DatasetCreate('rn_mass',rn_mass)
    call h5GroupClose
  end subroutine writeRenormalizeFactor

  subroutine writeDampingFactor
    if ( irank_mpi /= 0) return
    ! Renormalize factor
    call h5GroupOpen('parameters/boundary')
    call h5DatasetCreate('damping_rr', damping_rr      )
    call h5DatasetCreate('damping_rd', damping_rd      )
    call h5DatasetCreate('damping_fm_rr', damping_fm_rr)
    call h5DatasetCreate('damping_fm_rd', damping_fm_rd)
    call h5GroupClose
  end subroutine writeDampingFactor

  subroutine outputRawData(it,ex,ey,ez,by,bz,jx,jy,jz,rho,x,vx,vy,vz,iq)
    integer,intent(in) :: it
    double precision,dimension(:),intent(in) :: ex, ey, ez
    double precision,dimension(:),intent(in) :: by, bz
    double precision,dimension(:),intent(in) :: jx, jy, jz
    double precision,dimension(:),intent(in) :: rho
    double precision,dimension(:),intent(in) :: x
    double precision,dimension(:),intent(in) :: vx, vy, vz
    integer(1),      dimension(:),intent(in) :: iq

    integer :: ifirst, ilast
    integer :: is

    if ( irank_mpi /= 0) return

    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! Field
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    call h5GroupOpen('data')
    call h5GroupCreate(it2char(it), enable_grp_mv=.true.)
    call h5GroupCreate('particle')
    call h5GroupCreate('space_diag')
    call h5GroupCreate('box_diag')
    call h5GroupCreate('track_data')
    call h5GroupCreate('field', enable_grp_mv=.true.)
    call h5DatasetCreate('Ex',ex(2:nx+1))
    call h5DatasetCreate('Ey',ey(2:nx+1))
    call h5DatasetCreate('Ez',ez(2:nx+1))
    call h5DatasetCreate('By',by(2:nx+1))
    call h5DatasetCreate('Bz',bz(2:nx+1))
    call h5DatasetCreate('Jx',jx(2:nx+1))
    call h5DatasetCreate('Jy',jy(2:nx+1))
    call h5DatasetCreate('Jz',jz(2:nx+1))
    call h5DatasetCreate('rho',rho(2:nx+1))
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! particle
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ifirst = 0
    ilast  = 0

    do is=1, NS
      !	Initialize iterter
      ifirst = ilast + 1
      ilast  = ilast + NP(is)
      !Get 128 particles data
      ilast  = min(ifirst+128-1,ilast)
      call h5GroupOpen('data/'//it2char(it)//'/particle', is_abs=.true.)
      call h5GroupCreate(is2char(is), enable_grp_mv=.true.)
      call h5DatasetCreate('vx' ,vx(ifirst:ilast))
      call h5DatasetCreate('vy' ,vy(ifirst:ilast))
      call h5DatasetCreate('vz' ,vz(ifirst:ilast))
      call h5DatasetCreate('x'  ,x(ifirst:ilast))
      call h5DatasetCreate('iq' ,int(iq(ifirst:ilast)))
    end do
    call h5GroupClose
  end subroutine

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! KxKy Diag
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine outputKx(it, name, k, enable_create_dir)
    integer,intent(in) :: it
    character(*),intent(in) :: name
    double precision,dimension(:),intent(in) :: k
    logical :: enable_create_dir

    if ( irank_mpi /= 0) return

    call h5GroupOpen('kx')
    if (enable_create_dir) then
      call h5GroupCreate(it2char(it),enable_grp_mv=.true.)
    else
      call h5GroupOpen(it2char(it))
    endif

    call h5DatasetCreate(trim(name) ,k)
    call h5GroupClose
  end subroutine

  subroutine Box1d(it,is,ig,iw,var,name)
    integer,intent(in) :: it
    integer,intent(in) :: is,ig,iw
    double precision,dimension(:),intent(in) :: var
    character(*),intent(in)  :: name

    if (irank_mpi /= 0) return
    call h5GroupOpen('data/'//it2char(it)//'/box_diag')
    call h5DatasetCreate(is2char(is)//'_'//ig2char(ig)//'_'//ig2char(iw)//'_'//name, var)
    call h5GroupClose
  end subroutine

  subroutine Box2d(it,is,ig,iw,var,name)
    integer,intent(in) :: it
    integer,intent(in) :: is,ig,iw
    double precision,dimension(:,:),intent(in) :: var
    character(*),intent(in)  :: name

    if (irank_mpi /= 0) return
    call h5GroupOpen('data/'//it2char(it)//'/box_diag')
    call h5DatasetCreate(is2char(is)//'_'//ig2char(ig)//'_'//ig2char(iw)//'_'//name, var)
    call h5GroupClose
  end subroutine

  subroutine Box3d(it,is,ig,iw,var,name)
    integer,intent(in) :: it
    integer,intent(in) :: is,ig,iw
    double precision,dimension(:,:,:),intent(in) :: var
    character(*),intent(in)  :: name

    if (irank_mpi /= 0) return
    call h5GroupOpen('data/'//it2char(it)//'/box_diag')
    call h5DatasetCreate(is2char(is)//'_'//ig2char(ig)//'_'//ig2char(iw)//'_'//name, var)
    call h5GroupClose
  end subroutine

  subroutine Space1d(it,is,ig,iw,var,name)
    integer,intent(in) :: it
    integer,intent(in) :: is,ig,iw
    double precision,dimension(:),intent(in) :: var
    character(*),intent(in)  :: name

    if (irank_mpi /= 0) return
    call h5GroupOpen('data/'//it2char(it)//'/space_diag')
    call h5DatasetCreate(is2char(is)//'_'//ig2char(ig)//'_'//ig2char(iw)//'_'//name, var)
    call h5GroupClose
  end subroutine

  subroutine Space2d(it,is,ig,iw,var,name)
    integer,intent(in) :: it
    integer,intent(in) :: is,ig,iw
    double precision,dimension(:,:),intent(in) :: var
    character(*),intent(in)  :: name

    if (irank_mpi /= 0) return
    call h5GroupOpen('data/'//it2char(it)//'/space_diag')
    call h5DatasetCreate(is2char(is)//'_'//ig2char(ig)//'_'//ig2char(iw)//'_'//name, var)
    call h5GroupClose
  end subroutine

  subroutine Space3d(it,is,ig,iw,var,name)
    integer,intent(in) :: it
    integer,intent(in) :: is,ig,iw
    double precision,dimension(:,:,:),intent(in) :: var
    character(*),intent(in)  :: name

    if (irank_mpi /= 0) return
    call h5GroupOpen('data/'//it2char(it)//'/space_diag')
    call h5DatasetCreate(is2char(is)//'_'//ig2char(ig)//'_'//ig2char(iw)//'_'//name, var)
    call h5GroupClose
  end subroutine

  subroutine outputParticleTrack(it,ig,var,name)
    integer,intent(in) :: it
    integer,intent(in) :: ig
    double precision,dimension(:),intent(in) :: var
    character(*),intent(in)  :: name

    if (irank_mpi /= 0) return
    call h5GroupOpen('data/'//it2char(it)//'/track_data')
    call h5DatasetCreate(ig2char(ig)//'_'//name, var)
    call h5GroupClose
  end subroutine

  subroutine WriteEnergyMomentum
    if (irank_mpi /= 0) return
    call h5DatasetCreate('total/E_kinetic'     ,energy_kinetic)
    call h5DatasetCreate('total/E_kinetic_para',energy_kinetic_para)
    call h5DatasetCreate('total/E_kinetic_perp',energy_kinetic_perp)
    call h5DatasetCreate('total/E_e_field',energy_e_field)
    call h5DatasetCreate('total/E_b_field',energy_b_field)
    call h5DatasetCreate('total/E_total'  ,energy_total)
    call h5DatasetCreate('total/mu_total' ,magnetic_momentum_total)
  end subroutine WriteEnergyMomentum

  subroutine output1dField(it,ex)
    integer,intent(in) :: it
    double precision,dimension(:),intent(in) :: ex
    if (irank_mpi /= 0) return

    call h5GroupOpen('1d_field')
    call h5GroupCreate(it2char(it), enable_grp_mv=.true.)
    call h5DatasetCreate('Ex'     ,ex(2:nx+1:IXSKIP))
    call h5DatasetCreate('Ey_fwd' ,ey_fwd(2:nx+1:IXSKIP))
    call h5DatasetCreate('Ey_bwd' ,ey_bwd(2:nx+1:IXSKIP))
    call h5DatasetCreate('Ez_fwd' ,ez_fwd(2:nx+1:IXSKIP))
    call h5DatasetCreate('Ez_bwd' ,ez_bwd(2:nx+1:IXSKIP))
    call h5DatasetCreate('By_fwd' ,by_fwd(2:nx+1:IXSKIP)-by0)
    call h5DatasetCreate('By_bwd' ,by_bwd(2:nx+1:IXSKIP)-by0)
    call h5DatasetCreate('Bz_fwd' ,bz_fwd(2:nx+1:IXSKIP))
    call h5DatasetCreate('Bz_bwd' ,bz_bwd(2:nx+1:IXSKIP))
    call h5GroupClose
  end subroutine

  subroutine outputAllData(it,ex,ey,ez,by,bz,jx,jy,jz,rho,x,vx,vy,vz,iq)
    integer                      ,intent(in) :: it
    double precision,dimension(:),intent(in) :: ex, ey, ez
    double precision,dimension(:),intent(in) :: by, bz
    double precision,dimension(:),intent(in) :: jx, jy, jz
    double precision,dimension(:),intent(in) :: rho

    double precision,dimension(:),intent(in) :: x
    double precision,dimension(:),intent(in) :: vx, vy, vz
    integer(1),      dimension(:),intent(in) :: iq

    integer :: ifile_name

    character(len=5) :: crank
    character (:),allocatable,save :: datfile_name

#ifndef _DEBUG
    write(crank,'(i5.5)') irank_mpi
#else
    crank = '00000'
#endif
    ifile_name = len(file_name)-3 ! without extension
    allocate(character(ifile_name+10) :: datfile_name)
    datfile_name =  file_name(1:ifile_name) //'_'// crank // '.dat'

    open(99,file=datfile_name, form='unformatted', access='stream', status='replace')
    write(99) NX, cv, NS
    write(99) NP, omega_p, QM, v_para, v_perp, v_d, PITCH_DEG
    write(99) omega_c, THETA_DEG
    write(99) it
    write(99) ex, ey, ez, by, bz, jx, jy, jz, rho, rho0, x, vx, vy, vz, iq
    close(99)

    deallocate(datfile_name)
  end subroutine outputAllData

  subroutine outputFinalize
    deallocate(file_name)
    
  end subroutine outputFinalize

  subroutine outputMirrorInit(x,vx,vy,vz,iq)
    double precision,dimension(:),intent(in) :: x
    double precision,dimension(:),intent(in) :: vx, vy, vz
    integer(1),      dimension(:),intent(in) :: iq

    character(len=5) :: crank
    character (:),allocatable,save :: datfile_name

#ifndef _DEBUG
    write(crank,'(i5.5)') irank_mpi
#else
    crank = '00000'
#endif

    if (is_output == .true.) then
      allocate (character(len(file_name_out)-3+10) :: datfile_name)
      datfile_name = file_name_out(1:len(file_name_out)-3)//'_'// crank // '.dat'
    else
      allocate(character(11+10) :: datfile_name)
      datfile_name =  'mirror_init'//'_'// crank // '.dat'
    endif

    open(99,file=datfile_name, form='unformatted', access='stream', status='replace')
    write(99) NX, cv, NS
    write(99) NP, omega_p, QM, v_para, v_perp, v_d, PITCH_DEG
    write(99) omega_c, THETA_DEG
    write(99) x, vx, vy, vz, iq
    close(99)

    deallocate(datfile_name)
  end subroutine outputMirrorInit

  subroutine outputParticleIndex(izvp, npg, recvcounts, displs, index_list)
    integer,              intent(in) :: izvp, npg
    integer,dimension(:), intent(in) :: recvcounts, displs, index_list

    character(len=3 ) :: grp_z_v
    if ( irank_mpi /= 0) return

    write(grp_z_v,'(i3.3)') izvp
    call h5GroupOpen('tracker_index')
    call h5GroupCreate(grp_z_v, enable_grp_mv=.true.)
    call h5DatasetCreate('np'        , npg       )
    call h5DatasetCreate('recvcounts', recvcounts)
    call h5DatasetCreate('displs'    , displs    )
    call h5DatasetCreate('index_list', index_list)
    call h5GroupClose
  end subroutine outputParticleIndex
end module OutputHDF5
