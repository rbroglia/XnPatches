!> This module contains the definition of procedures and variables useful for post-process the Xnavis data.
module Data_Type_PostProcess
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision                     ! Integers and reals precision definition.
USE Data_Type_Command_Line_Interface ! Definition of Type_Command_Line_Interface.
USE Data_Type_OS                     ! Definition of Type_OS.
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
  logical::                   fcc            = .false.   !< Inquiring flag for icc file.
  logical::                   ngc            = .false.   !< Inquiring flag for grd without ghosts cells.
  logical::                   sol            = .false.   !< Inquiring flag for solution file.
  logical::                   cell           = .false.   !< Inquiring flag for interpolation (or not) variables at nodes.
  logical::                   level_set      = .false.   !< Inquiring flag for level set variable.
  logical::                   laminar        = .false.   !< Inquiring flag for no turbulent model.
  logical::                   zeroeq         = .false.   !< Inquiring flag for zero equations turbulent variables.
  logical::                   oneeq          = .true.    !< Inquiring flag for one  equations turbulent variables.
  logical::                   twoeq          = .false.   !< Inquiring flag for two  equations turbulent variables.
  logical::                   vordet         = .false.   !< Inquiring flag for vordet variable computing.
  logical::                   binary         = .true.    !< Inquiring flag for binary output file.
  logical::                   tec            = .true.    !< Inquiring flag for tecplot file format.
  logical::                   vtk            = .false.   !< Inquiring flag for vtk file format.
  logical::                   global_blk_num = .false.   !< Flag for inquiring the activation of global blocks numeration.
  logical::                   verbose        = .false.   !< Verbose output.
  type(Type_OS)::             OS                         !< Running architecture.
  integer(I4P), allocatable:: blockmap(:,:),procmap(:,:) !< Blocks maps.
  contains
    procedure:: set_from_cli       ! Procedure for setting post process data from CLI.
    procedure:: set_from_mbpar     ! Procedure for setting post process data from mb.par if found.
    procedure:: set_from_procinput ! Procedure for setting global block numeration if proc.input is found.
    procedure:: set_blockmap       ! Procedure for setting block maps if local numeration is used.
    procedure:: input_files_init   ! Procedure for initializing input files.
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
  call cli%get(switch='-g',     val=pp%File_grd, error=error,pref='|-->'); if (error/=0) stop
  call cli%get(switch='-o',     val=pp%File_out, error=error,pref='|-->'); if (error/=0) stop
  call cli%get(switch='-i',     val=pp%File_icc, error=error,pref='|-->'); if (error/=0) stop
  call cli%get(switch='-s',     val=pp%File_sol, error=error,pref='|-->'); if (error/=0) stop
  call cli%get(switch='-ngc',   val=pp%ngc,      error=error,pref='|-->'); if (error/=0) stop
  call cli%get(switch='-cell',  val=pp%cell,     error=error,pref='|-->'); if (error/=0) stop
  call cli%get(switch='-ls',    val=pp%level_set,error=error,pref='|-->'); if (error/=0) stop
  call cli%get(switch='-nt',    val=pp%laminar,  error=error,pref='|-->'); if (error/=0) stop
  call cli%get(switch='-eq',    val=turb_eq,     error=error,pref='|-->'); if (error/=0) stop
  call cli%get(switch='-vordet',val=pp%vordet,   error=error,pref='|-->'); if (error/=0) stop
  call cli%get(switch='-ascii', val=ascii,       error=error,pref='|-->'); if (error/=0) stop
  call cli%get(switch='-tec',   val=tec_f,       error=error,pref='|-->'); if (error/=0) stop
  call cli%get(switch='-vtk',   val=vtk_f,       error=error,pref='|-->'); if (error/=0) stop
  call cli%get(switch='-proc',  val=pp%myrank,   error=error,pref='|-->'); if (error/=0) stop
  call cli%get(switch='-os',    val=os_type,     error=error,pref='|-->'); if (error/=0) stop
  call cli%get(switch='-v',     val=pp%verbose,  error=error,pref='|-->'); if (error/=0) stop
  ! using CLA values for driving XnPlot
  if (adjustl(trim(pp%File_icc))/='unset') pp%fcc = .true.
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
  if ((pp%vordet.and.(.not.pp%fcc)).or.(pp%vordet.and.(.not.pp%sol))) then
    write(stderr,'(A)') '+--> In order to compute "vordet" variables the icc and sol files must be provided.'
    stop
  endif
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
  real(R8P)::                              Fr       !< Froude number.
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
    read(unitfree,*)        ! Reynolds number
    read(unitfree,*) Fr     ! Foude number
    if (Fr>0._R8P) pp%level_set = .true.
    read(unitfree,*)
    read(unitfree,*) string ! turbulence model
    close(unitfree)
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
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  inquire(file=adjustl(trim(pp%File_grd)),exist=is_file,iostat=error)
  if (.NOT.is_file) then
    error = File_Not_Found(filename=pp%File_grd,cpn='input_files_init')
  else
    open(unit=Get_Unit(pp%unit_grd),file=adjustl(trim(pp%File_grd)),form='UNFORMATTED',action='READ')
  endif
  if (pp%fcc) then
    inquire(file=adjustl(trim(pp%File_icc)),exist=is_file,iostat=error)
    if (.NOT.is_file) then
      error = File_Not_Found(filename=pp%File_icc,cpn='input_files_init')
    else
      open(unit=Get_Unit(pp%unit_icc),file=adjustl(trim(pp%File_icc)),form='UNFORMATTED',action='READ')
    endif
  endif
  if (pp%sol) then
    inquire(file=adjustl(trim(pp%File_sol)),exist=is_file,iostat=error)
    if (.NOT.is_file) then
      error = File_Not_Found(filename=pp%File_sol,cpn='input_files_init')
    else
      open(unit=Get_Unit(pp%unit_sol),file=adjustl(trim(pp%File_sol)),form='UNFORMATTED',action='READ')
    endif
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine input_files_init

  !> @brief Procedure for finalizing post-processor.
  subroutine finalize(pp)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_PostProcess), intent(INOUT):: pp !< Post-processing options.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  close(pp%unit_grd)
  if (pp%fcc) then
    close(pp%unit_icc)
  endif
  if (pp%sol) then
    close(pp%unit_sol)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine finalize
  !> @}
endmodule Data_Type_PostProcess
