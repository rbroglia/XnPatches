!> This module contains the definition of procedures and variables useful for post-process the Xnavis data.
module Data_Type_PostProcess
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision                     ! Integers and reals precision definition.
USE Data_Type_Command_Line_Interface ! Definition of Type_Command_Line_Interface.
USE Data_Type_OS                     ! Definition of Type_OS.
USE Data_Type_Vector                 ! Definition of Type_Vector.
USE Lib_IO_Misc                      ! Library for miscellanea IO procedures.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> Derived type containing the post-processing options.
!> @ingroup Data_Type_PostProcessDerivedType
type, public:: Type_PostProcess
  integer(I4P)::              myrank         = 0_I4P     !< Rank of current processed file.
  character(100)::            File_grd       = 'unset'   !< GRD file name.
  character(100)::            File_icc       = 'unset'   !< ICC file name.
  character(100)::            File_sol       = 'unset'   !< Solution file name.
  character(100)::            File_out       = 'unset'   !< Output file name.
  integer(I4P)::              unit_out       = 0_I4P     !< Unit of Output file.
  integer(I4P)::              unit_grd       = 0_I4P     !< Unit of GRD file.
  integer(I4P)::              unit_icc       = 0_I4P     !< Unit of ICC file.
  integer(I4P)::              unit_sol       = 0_I4P     !< Unit of Solution file.
  integer(I4P)::              unit_for       = 0_I4P     !< Unit of forces.dat file.
  integer(I4P)::              Ng             = 0_I4P     !< Number of groups.
  integer(I4P), allocatable:: unit_g_for(:)              !< Unit of forces.dat file for each group.
  integer(I4P)::              patch          = 1_I4P     !< Boundary condition value of post-processed patches.
  integer(I4P)::              s_offset       = 0_I4P     !< Streamline patch offset.
  real(R8P)::                 Re             = -1._R_P   !< Reynolds number.
  real(R8P)::                 Fr             = -1._R_P   !< Froude number.
  real(R8P)::                 rFr2           = 0._R_P    !< 1/(Froude number)^2.
  real(R8P)::                 zfs            = 1.D100    !< Z quote of free surface.
  type(Type_Vector)::         fsum_p                     !< Global sum of forces vector, pressure part.
  type(Type_Vector)::         fsum_v                     !< Global sum of forces vector, viscous part.
  type(Type_Vector)::         msum_p                     !< Global sum of moments vector, pressure part.
  type(Type_Vector)::         msum_v                     !< Global sum of moments vector, viscous part.
  real(R8P)::                 Ssum           = 0._R8P    !< Global sum of "wet" surface.
  logical::                   sol            = .false.   !< Inquiring flag for solution file.
  logical::                   cell           = .false.   !< Inquiring flag for interpolation (or not) variables at nodes.
  logical::                   level_set      = .false.   !< Inquiring flag for level set variable.
  logical::                   laminar        = .false.   !< Inquiring flag for no turbulent model.
  logical::                   zeroeq         = .false.   !< Inquiring flag for zero equations turbulent variables.
  logical::                   oneeq          = .true.    !< Inquiring flag for one  equations turbulent variables.
  logical::                   twoeq          = .false.   !< Inquiring flag for two  equations turbulent variables.
  logical::                   binary         = .true.    !< Inquiring flag for binary output file.
  logical::                   tec            = .true.    !< Inquiring flag for tecplot file format.
  logical::                   vtk            = .false.   !< Inquiring flag for vtk file format.
  logical::                   global_blk_num = .false.   !< Flag for inquiring the activation of global blocks numeration.
  logical::                   forces         = .false.   !< Inquiring flag for forces computing.
  logical::                   yp             = .false.   !< Inquiring flag for yplus computing.
  logical::                   metrics        = .false.   !< Inquiring flag for metrics saving.
  logical::                   verbose        = .false.   !< Verbose output.
  type(Type_OS)::             OS                         !< Running architecture.
  integer(I4P), allocatable:: blockmap(:,:),procmap(:,:) !< Blocks maps.
  contains
    procedure:: set_from_cli       ! Procedure for setting post process data from CLI.
    procedure:: set_from_mbpar     ! Procedure for setting post process data from mb.par if found.
    procedure:: set_from_procinput ! Procedure for setting global block numeration if proc.input is found.
    procedure:: set_blockmap       ! Procedure for setting block maps if local numeration is used.
    procedure:: input_files_init   ! Procedure for initializing input files.
    procedure:: save_forces        ! Procedure for saving forces global integral.
    procedure:: finalize           ! Procedure for finalizing post-processor.
endtype Type_PostProcess
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @ingroup Data_Type_PostProcessPrivateProcedure
  !> @{
  !> @brief Procedure for setting post process data from CLA.
  subroutine set_from_cli(pp,cli)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_PostProcess),           intent(INOUT):: pp       !< Post-processing options.
  type(Type_Command_Line_Interface), intent(INOUT):: cli      !< Command Line Interface (CLI).
  integer(I4P)::                                     turb_eq  !< Equations of Turbulent model.
  logical::                                          ascii    !< Flag for setting ascii output files.
  character(3)::                                     tec_f    !< Flag for setting Tecplot output files.
  character(3)::                                     vtk_f    !< Flag for setting VTK output files.
  character(3)::                                     os_type  !< Type of OS.
  integer(I4P)::                                     error    !< Error trapping flag.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! getting CLA values
  call cli%get(switch='-g',       val=pp%File_grd, error=error,pref='|-->'); if (error/=0) stop
  call cli%get(switch='-o',       val=pp%File_out, error=error,pref='|-->'); if (error/=0) stop
  call cli%get(switch='-i',       val=pp%File_icc, error=error,pref='|-->'); if (error/=0) stop
  call cli%get(switch='-s',       val=pp%File_sol, error=error,pref='|-->'); if (error/=0) stop
  call cli%get(switch='-p',       val=pp%patch,    error=error,pref='|-->'); if (error/=0) stop
  call cli%get(switch='-Fr',      val=pp%Fr,       error=error,pref='|-->'); if (error/=0) stop
  call cli%get(switch='-Re',      val=pp%Re,       error=error,pref='|-->'); if (error/=0) stop
  call cli%get(switch='-zfs',     val=pp%zfs,      error=error,pref='|-->'); if (error/=0) stop
  call cli%get(switch='-cell',    val=pp%cell,     error=error,pref='|-->'); if (error/=0) stop
  call cli%get(switch='-forces',  val=pp%forces,   error=error,pref='|-->'); if (error/=0) stop
  call cli%get(switch='-metrics', val=pp%metrics,  error=error,pref='|-->'); if (error/=0) stop
  call cli%get(switch='-yplus',   val=pp%yp,       error=error,pref='|-->'); if (error/=0) stop
  call cli%get(switch='-ls',      val=pp%level_set,error=error,pref='|-->'); if (error/=0) stop
  call cli%get(switch='-nt',      val=pp%laminar,  error=error,pref='|-->'); if (error/=0) stop
  call cli%get(switch='-eq',      val=turb_eq,     error=error,pref='|-->'); if (error/=0) stop
  call cli%get(switch='-stream',  val=pp%s_offset, error=error,pref='|-->'); if (error/=0) stop
  call cli%get(switch='-ascii',   val=ascii,       error=error,pref='|-->'); if (error/=0) stop
  call cli%get(switch='-tec',     val=tec_f,       error=error,pref='|-->'); if (error/=0) stop
  call cli%get(switch='-vtk',     val=vtk_f,       error=error,pref='|-->'); if (error/=0) stop
  call cli%get(switch='-proc',    val=pp%myrank,   error=error,pref='|-->'); if (error/=0) stop
  call cli%get(switch='-os',      val=os_type,     error=error,pref='|-->'); if (error/=0) stop
  call cli%get(switch='-v',       val=pp%verbose,  error=error,pref='|-->'); if (error/=0) stop
  ! using CLA values for driving XnPlot
  if (adjustl(trim(pp%File_sol))/='unset') pp%sol = .true.
  if (adjustl(trim(pp%File_out))=='unset') pp%File_out = adjustl(trim(pp%File_grd))
  if (pp%laminar) then
    pp%zeroeq = .false.
     pp%oneeq = .false.
     pp%twoeq = .false.
  endif
  select case(turb_eq)
  case(0)
    pp%zeroeq = .true.
    pp%oneeq = .false.
  case(1)
    pp%oneeq = .true.
  case(2)
    pp%twoeq = .true.
    pp%oneeq = .false.
  endselect
  if (pp%Fr>0) pp%rFr2 = 1._R8P/(pp%Fr*pp%Fr)
  if (ascii) pp%binary = .false.
  pp%tec = (Upper_Case(adjustl(trim(tec_f)))=='YES')
  pp%vtk = (Upper_Case(adjustl(trim(vtk_f)))=='YES')
  if (pp%myrank/=-1) then
    pp%global_blk_num = .true.
  else
    pp%myrank = 0
  endif
  call pp%OS%init(c_id=Upper_Case(os_type))
  if ( pp%laminar .and. (pp%zeroeq .or. pp%oneeq .or. pp%twoeq) ) then
     write(stderr,'(A)') '+--> Incompatible switches, laminar disables turbulent model'
     stop
  end if
  ! converting the directory separators of files names according to the OS using
  pp%File_grd = adjustl(trim(pp%File_grd)) ; pp%File_grd = pp%OS%string_separator_fix(string=pp%File_grd)
  pp%File_icc = adjustl(trim(pp%File_icc)) ; pp%File_icc = pp%OS%string_separator_fix(string=pp%File_icc)
  pp%File_sol = adjustl(trim(pp%File_sol)) ; pp%File_sol = pp%OS%string_separator_fix(string=pp%File_sol)
  pp%File_out = adjustl(trim(pp%File_out)) ; pp%File_out = pp%OS%string_separator_fix(string=pp%File_out)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine set_from_cli

  !> @brief Procedure for setting post process data from mb.par if found.
  subroutine set_from_mbpar(pp)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_PostProcess), intent(INOUT):: pp       !< Post-processing options.
  logical::                                is_file  !< Flag for inquiring the presence of file.
  integer(I4P)::                           error    !< Error trapping flag.
  integer(I4P)::                           unitfree !< Free logical unit.
  integer(I4P)::                           e        !< Counter.
  character(100)::                         string   !< Dummy string.
  logical::                                balom
  logical::                                sgs
  logical::                                spall
  logical::                                des
  logical::                                ddes
  logical::                                lambr
  logical::                                chang
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  balom = .false.
  sgs   = .false.
  spall = .false.
  des   = .false.
  ddes  = .false.
  lambr = .false.
  chang = .false.
  inquire(file="mb.par",exist=is_file,iostat=error)
  if (is_file) then
    open(unit=Get_Unit(unitfree),file="mb.par",action='READ')
    do e=1,8 ! skip the first 8 records
      read(unitfree,*)
    enddo
    read(unitfree,*) e      ! number grid levels
    do e=1,e+9 ! skip the e+9 records
      read(unitfree,*)
    enddo
    read(unitfree,*) pp%Re  ! Reynolds number
    read(unitfree,*) pp%Fr  ! Foude number
    read(unitfree,*)
    read(unitfree,*) string ! turbulence model
    read(unitfree,*)
    read(unitfree,*) pp%zfs ! quote of free surface
    close(unitfree)
    if (pp%Fr>0._R_P) then
      pp%level_set = .true.
      pp%rFr2 = 1._R_P/(pp%Fr*pp%Fr)
    endif
    ! checking turbulence model
    if (string(1:3)=='bal'.or.string(1:3)=='BAL') balom = .true.
    if (string(1:3)=='sgs'.or.string(1:3)=='SGS') sgs   = .true.
    if (string(1:3)=='spa'.or.string(1:3)=='SPA') spall = .true.
    if (string(1:3)=='lam'.or.string(1:3)=='LAM') lambr = .true.
    if (string(1:3)=='cha'.or.string(1:3)=='CHA') chang = .true.
    if (string(1:3)=='des'.or.string(1:3)=='DES') then
       spall = .true.
       des   = .true.
    endif
    if (string(1:3)=='dde'.or.string(1:3)=='DDE') then
       spall = .true.
       des   = .true.
       ddes  = .true.
    endif
    pp%zeroeq = (sgs.OR.balom)
    pp%oneeq  = spall
    pp%twoeq  = (lambr.OR.chang)
    write(stdout,'(A)') '+--> Found mb.par'
    write(stdout,'(A)') '|--> Reynolds number '//str(n=pp%Re)
    write(stdout,'(A)') '|--> Froude   number '//str(n=pp%Fr)
    if (pp%level_set) then
      write(stdout,'(A)') '|--> Free sruface present'
    else
      write(stdout,'(A)') '|--> Free sruface not present'
    endif
    if (pp%zeroeq) then
      write(stdout,'(A)') '|--> Zero equation turbulence model'
    elseif (pp%oneeq) then
      write(stdout,'(A)') '|--> One equation turbulence model'
    elseif(pp%twoeq) then
      write(stdout,'(A)') '|--> Two equations turbulence model'
    endif
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine set_from_mbpar

  !> @brief Procedure for setting global block numeration if proc.input is found.
  subroutine set_from_procinput(pp)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_PostProcess), intent(INOUT):: pp       !< Post-processing options.
  logical::                                is_file  !< Flag for inquiring the presence of file.
  integer(I4P)::                           error    !< Error trapping flag.
  integer(I4P)::                           unitfree !< Free logical unit.
  integer(I4P)::                           b,c      !< Counters.
  integer(I4P)::                           nprocs   !< number of processes
  integer(I4P)::                           Nb       !< Number of local blocks.
  integer(I4P)::                           Nbtot    !< Number of total blocks.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (pp%global_blk_num) then
    inquire(file='proc.input',exist=is_file)
    if (.not.is_file) then
      error = File_Not_Found(filename='proc.input',cpn='parse_command_arguments')
    else
      write(stdout,'(A)') '+--> Found proc.input'
      open(unit=Get_Unit(unitfree),file="proc.input",action='READ')
      read(unitfree,*) ! record skipped
      read(unitfree,*) ! record skipped
      read(unitfree,*) ! record skipped
      read(unitfree,*) Nbtot ! number of total blocks
      if (allocated(pp%procmap)) deallocate(pp%procmap) ; allocate(pp%procmap(1:2,1:Nbtot)) ; pp%procmap = 0_I_P
      read(unitfree,*) ! record skipped
      read(unitfree,*) ! record skipped
      read(unitfree,*) ! record skipped
      do b=1,Nbtot
        read(unitfree,*) c,pp%procmap(1,b),c,pp%procmap(2,b)
      end do
      close(unitfree)
      nprocs = maxval(pp%procmap(2,:))
      pp%Ng  = maxval(pp%procmap(1,:))
      ! computing the local (of myrank) number of blocks
      Nb = 0
      do b=1,Nbtot
        if (pp%procmap(2,b)==pp%myrank) Nb = Nb + 1
      enddo
      if (allocated(pp%blockmap)) deallocate(pp%blockmap) ; allocate(pp%blockmap(1:2,1:Nb)) ; pp%blockmap = 0_I_P
      c = 0
      do b=1,Nbtot
        if (pp%procmap(2,b)==pp%myrank) then
          c = c + 1
          pp%blockmap(1,c) = pp%procmap(1,b)
          pp%blockmap(2,c) = b
        endif
      enddo
    endif
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine set_from_procinput

  !> @brief Procedure for setting block maps if local numeration is used.
  subroutine set_blockmap(pp,Nb)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_PostProcess), intent(INOUT):: pp !< Post-processing options.
  integer(I4P),            intent(IN)::    Nb !< Number of local blocks.
  integer(I4P)::                           b  !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(pp%blockmap)) deallocate(pp%blockmap) ; allocate(pp%blockmap(1:2,1:Nb)) ; pp%blockmap = 0_I_P
  do b=1,Nb
    pp%blockmap(2,b) = b
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine set_blockmap

  !> @brief Procedure for initializing input files.
  subroutine input_files_init(pp)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_PostProcess), intent(INOUT):: pp      !< Post-processing options.
  logical::                                is_file !< Flag for inquiring the presence of file.
  integer(I4P)::                           error   !< Error trapping flag.
  integer(I4P)::                           g       !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  inquire(file=adjustl(trim(pp%File_grd)),exist=is_file,iostat=error)
  if (.NOT.is_file) then
    error = File_Not_Found(filename=pp%File_grd,cpn='input_files_init')
  else
    open(unit=Get_Unit(pp%unit_grd),file=adjustl(trim(pp%File_grd)),form='UNFORMATTED',action='READ')
  endif
  inquire(file=adjustl(trim(pp%File_icc)),exist=is_file,iostat=error)
  if (.NOT.is_file) then
    error = File_Not_Found(filename=pp%File_icc,cpn='input_files_init')
  else
    open(unit=Get_Unit(pp%unit_icc),file=adjustl(trim(pp%File_icc)),form='UNFORMATTED',action='READ')
  endif
  if (pp%sol) then
    inquire(file=adjustl(trim(pp%File_sol)),exist=is_file,iostat=error)
    if (.NOT.is_file) then
      error = File_Not_Found(filename=pp%File_sol,cpn='input_files_init')
    else
      open(unit=Get_Unit(pp%unit_sol),file=adjustl(trim(pp%File_sol)),form='UNFORMATTED',action='READ')
    endif
  endif
  if (pp%forces) then
    !if (pp%forcesRB) then
    !  open(unit=Get_Unit(unit_for_RB),file=adjustl(trim(File_out))//"-forces.RB",form='UNFORMATTED')
    !  open(unit=Get_Unit(unit_for_RB_scr),form='UNFORMATTED',status='SCRATCH')
    !endif
    if (pp%Ng>0) then
      if (allocated(pp%unit_g_for)) deallocate(pp%unit_g_for) ; allocate(pp%unit_g_for(0:pp%Ng))
      do g=0,pp%Ng
        open(unit=Get_Unit(pp%unit_g_for(g)),file=adjustl(trim(pp%File_out))//"-forces-grp_"//trim(strz(3,g))//".dat")
      enddo
    endif
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine input_files_init

  !> @brief Procedure for saving forces global integral.
  function save_forces(pp) result(error)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_PostProcess), intent(INOUT):: pp    !< Post-processing options.
  integer(I4P)::                           error !< Error trapping flag.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  open(unit=Get_Unit(pp%unit_for),file=adjustl(trim(pp%File_out))//"-forces.dat")
  write(pp%unit_for,'(A)',iostat=error)trim(str(n=(pp%fsum_v%x+pp%fsum_p%x)))//' '// &
                                       trim(str(n=(pp%fsum_v%y+pp%fsum_p%y)))//' '// &
                                       trim(str(n=(pp%fsum_v%z+pp%fsum_p%z)))//' '// &
                                       trim(str(n=(            pp%fsum_p%x)))//' '// &
                                       trim(str(n=(            pp%fsum_p%y)))//' '// &
                                       trim(str(n=(            pp%fsum_p%z)))//' '// &
                                       trim(str(n=(pp%fsum_v%x            )))//' '// &
                                       trim(str(n=(pp%fsum_v%y            )))//' '// &
                                       trim(str(n=(pp%fsum_v%z            )))//' '// &
                                       trim(str(n=(pp%msum_v%x+pp%msum_p%x)))//' '// &
                                       trim(str(n=(pp%msum_v%y+pp%msum_p%y)))//' '// &
                                       trim(str(n=(pp%msum_v%z+pp%msum_p%z)))//' '// &
                                       trim(str(n=(            pp%msum_p%x)))//' '// &
                                       trim(str(n=(            pp%msum_p%y)))//' '// &
                                       trim(str(n=(            pp%msum_p%z)))//' '// &
                                       trim(str(n=(pp%msum_v%x            )))//' '// &
                                       trim(str(n=(pp%msum_v%y            )))//' '// &
                                       trim(str(n=(pp%msum_v%z            )))//' '// &
                                       trim(str(n=(pp%Ssum                )))//' '// &
                                       'Fx,Fy,Fz,Fx_p,Fy_p,Fz_p,Fx_v,Fy_v,Fz_v,Mx,My,Mz,Mx_p,My_p,Mz_p,Mx_v,My_v,Mz_v,S of file "'&
                                       //adjustl(trim(pp%File_sol))//'"'
  close(pp%unit_for)
!  if (forcesRB) then
    !if (Np>0) then
      !if (allocated(Ni)) deallocate(Ni) ; allocate(Ni(1:1))
      !if (allocated(Nj)) deallocate(Nj) ; allocate(Nj(1:1))
      !if (allocated(Nk)) deallocate(Nk) ; allocate(Nk(1:1))
      !write(unit_for_RB)Np
      !rewind(unit_for_RB_scr)
      !do
        !read(unit_for_RB_scr,err=10,end=10)b,v,Ni(1),Nj(1),Nk(1)
        !write(unit_for_RB)b,v,Ni(1),Nj(1),Nk(1)
        !if (allocated(dummy)) deallocate(dummy) ; allocate(dummy(1:Ni(1),1:Nj(1),1:Nk(1)))
        !do i=1,12
          !read(unit_for_RB_scr,err=10)dummy
        !enddo
      !enddo
      !10 continue
      !rewind(unit_for_RB_scr)
      !do
        !read(unit_for_RB_scr,err=20,end=20)b,v,Ni(1),Nj(1),Nk(1)
        !if (allocated(dummy)) deallocate(dummy) ; allocate(dummy(1:Ni(1),1:Nj(1),1:Nk(1)))
        !do i=1,12
          !read(unit_for_RB_scr,err=20,end=20)dummy
          !write(unit_for_RB)dummy
        !enddo
      !enddo
      !20 continue
    !endif
    !close(unit_for_RB)     ! Riccardo Broglia forces.RB
    !close(unit_for_RB_scr) ! Riccardo Broglia forces.RB scratch file
  !endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction save_forces

  !> @brief Procedure for finalizing post-processor.
  subroutine finalize(pp)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_PostProcess), intent(INOUT):: pp !< Post-processing options.
  integer(I4P)::                           g  !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  close(pp%unit_grd)
  close(pp%unit_icc)
  if (pp%sol) then
    close(pp%unit_sol)
  endif
  if (pp%Ng>0.and.pp%forces) then
    do g=0,pp%Ng
      close(pp%unit_g_for(g))
    enddo
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine finalize
  !> @}
endmodule Data_Type_PostProcess
