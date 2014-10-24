!> @brief This module contains procedures for saving Tecplot output.
module Lib_TEC
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision          ! Integers and reals precision definition.
USE Block_Variables       ! Block variables definition.
USE Data_Type_PostProcess ! Definition of Type_PostProcess.
USE Data_Type_Vector      ! Definition of Type_Vector.
USE Lib_IO_Misc           ! Library for miscellanea IO procedures.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
save
public:: file_tec
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
type:: Type_File_Tec
  character(1)::                  tecendrec = char(0) !< End-character for binary-record end.
  character(len=:), allocatable:: tecvarname          !< Variables name for tecplot header file.
  integer(I4P),     allocatable:: tecvarloc(:)        !< Tecplot array of variables location.
  character(len=:), allocatable:: tecvarlocstr        !< Tecplot string of variables location.
  integer(I4P),     allocatable:: tecnull(:)          !< Tecplot null array.
  integer(I4P)::                  nvar = 3            !< Number of variables saved.
  contains
    procedure::         init       ! Procedure for initializing file writing.
    procedure::         save_block ! Procedure for saving one block.
    procedure, nopass:: finalize   ! Procedure for finalizing file writing.
endtype Type_File_Tec
type(Type_File_Tec):: file_tec
! tecio functions
integer(I4P), external::  tecini112,    &
                          tecauxstr112, &
                          teczne112,    &
                          tecdat112,    &
                          tecend112
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @brief Procedure for initializing file writing.
  subroutine init(file_d,pp)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_File_Tec),   intent(INOUT):: file_d !< File data.
  type(Type_PostProcess), intent(INOUT):: pp     !< Post-processor data.
  integer(I4P)::                          error  !< Error trapping flag.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (pp%binary) then
    file_d%tecvarname = 'x y z'
  else
    file_d%tecvarname = ' VARIABLES ="x" "y" "z"'
  endif
  if (pp%fcc) then
    file_d%nvar = file_d%nvar + 1
    if (pp%binary) then
      file_d%tecvarname = trim(file_d%tecvarname)//' icc'
    else
      file_d%tecvarname = trim(file_d%tecvarname)//' "icc"'
    endif
  endif
  if (pp%sol) then
    file_d%nvar = file_d%nvar + 4
    if (pp%binary) then
      file_d%tecvarname = trim(file_d%tecvarname)//' u v w p'
    else
      file_d%tecvarname = trim(file_d%tecvarname)//' "u" "v" "w" "p"'
    endif
    if (pp%zeroeq) then
      file_d%nvar = file_d%nvar + 1 ! visc must be saved
      if (pp%binary) then
        file_d%tecvarname = trim(file_d%tecvarname)//' visc'
      else
        file_d%tecvarname = trim(file_d%tecvarname)//' "visc"'
      endif
    elseif (pp%oneeq) then
      file_d%nvar = file_d%nvar + 2 ! visc and vitl must be saved
      if (pp%binary) then
        file_d%tecvarname = trim(file_d%tecvarname)//' visc vitl'
      else
        file_d%tecvarname = trim(file_d%tecvarname)//' "visc" "vitl"'
      endif
    elseif (pp%twoeq) then
      file_d%nvar = file_d%nvar + 3 ! visc, ken and eps must be saved
      if (pp%binary) then
        file_d%tecvarname = trim(file_d%tecvarname)//' visc ken eps'
      else
        file_d%tecvarname = trim(file_d%tecvarname)//' "visc" "ken" "eps"'
      endif
    endif
    if (pp%level_set) then
      file_d%nvar = file_d%nvar + 2 ! f and f0 must be saved
      if (pp%binary) then
        file_d%tecvarname = trim(file_d%tecvarname)//' f f0'
      else
        file_d%tecvarname = trim(file_d%tecvarname)//' "f" "f0"'
      endif
    endif
    if (pp%vordet) then
      file_d%nvar = file_d%nvar + 6 ! vorticity variables must be saved
      if (pp%binary) then
        file_d%tecvarname = trim(file_d%tecvarname)//' vordet qfactor helicity vorticity-x vorticity-y vorticity-z'
      else
        file_d%tecvarname = trim(file_d%tecvarname)//' "vordet" "qfactor" "helicity" "vorticity-x" "vorticity-y" "vorticity-z"'
      endif
    endif
  endif
  if (allocated(file_d%tecnull  )) deallocate(file_d%tecnull  ) ; allocate(file_d%tecnull(  1:file_d%nvar))
  if (allocated(file_d%tecvarloc)) deallocate(file_d%tecvarloc) ; allocate(file_d%tecvarloc(1:file_d%nvar))
  if (pp%binary) then
    file_d%tecnull = 0
    file_d%tecvarloc = 1
    if (file_d%nvar>3.and.pp%cell) file_d%tecvarloc(4:file_d%nvar)= 0
  endif
  if (pp%binary) then
    error = tecini112(file_d%tecendrec,                          &
                      trim(file_d%tecvarname)//file_d%tecendrec, &
                      adjustl(trim(pp%File_out))//".plt"//file_d%tecendrec,'.'//file_d%tecendrec,0,0,1)
  else
    open(unit=Get_Unit(pp%unit_out),file=adjustl(trim(pp%File_out))//".dat")
    write(pp%unit_out,'(A)',iostat=error)trim(file_d%tecvarname)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine init

  !> @brief Procedure for saving block.
  function save_block(file_d,pp,b,Ni,Nj,Nk) result(error)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_File_Tec),   intent(INOUT):: file_d                               !< File data.
  type(Type_PostProcess), intent(IN)::    pp                                   !< Post-processor data.
  integer(I4P),           intent(IN)::    b                                    !< Actual block number.
  integer(I4P),           intent(IN)::    Ni,Nj,Nk                             !< Number of cells.
  integer(I4P)::                          error                                !< Error trapping flag.
  integer(I4P)::                          nnode,ncell                          !< Number of nodes and cells.
  real(R8P), allocatable::                vari(:,:,:)                          !< Interpolated generic variables.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.pp%cell) then ! variables must be interpolated at nodes, allocating dummy variable
    if (allocated(vari)) deallocate(vari) ; allocate(vari(0:Ni,0:Nj,0:Nk))
  endif
  nnode = (Ni+1)*(Nj+1)*(Nk+1) ! computing number of nodes
  ncell =  Ni   * Nj   * Nk    ! computing number of cells
  associate(binary=>pp%binary,cell=>pp%cell,blockmap=>pp%blockmap,sol=>pp%sol,fcc=>pp%fcc,vordet=>pp%vordet,&
            zeroeq=>pp%zeroeq,oneeq=>pp%oneeq,twoeq=>pp%twoeq,level_set=>pp%level_set,unit_out=>pp%unit_out,nvar=>file_d%nvar)
  if (binary) then
    error = teczne112('blk_'//trim(strz(4,blockmap(2,b)))//'_grp_'//trim(strz(3,blockmap(1,b)))//file_d%tecendrec, &
                      0,                                                                                           &
                      Ni+1,                                                                                        &
                      Nj+1,                                                                                        &
                      Nk+1,                                                                                        &
                      0,                                                                                           &
                      0,                                                                                           &
                      0,                                                                                           &
                      0.0,                                                                                         &
                      0,                                                                                           &
                      0,                                                                                           &
                      1,                                                                                           &
                      0,                                                                                           &
                      0,                                                                                           &
                      0,                                                                                           &
                      0,                                                                                           &
                      0,                                                                                           &
                      file_d%tecnull,                                                                              &
                      file_d%tecvarloc,                                                                            &
                      file_d%tecnull,                                                                              &
                      0)
  else
    if (nvar>3.and.cell) then
      write(unit_out,'(A)',iostat=error)' ZONE  T="blk_'//trim(strz(4,blockmap(2,b)))//'_grp_'//trim(strz(3,blockmap(1,b)))//'"'//&
                                        ', I='//trim(str(no_sign=.true.,n=Ni+1))//                                                &
                                        ', J='//trim(str(no_sign=.true.,n=Nj+1))//                                                &
                                        ', K='//trim(str(no_sign=.true.,n=Nk+1))//                                                &
                                        ', DATAPACKING=BLOCK'//                                                                   &
                                        ', VARLOCATION=([1-3]=NODAL,[4-'//trim(str(.true.,nvar))//']=CELLCENTERED)'
    else
      write(unit_out,'(A)',iostat=error)' ZONE  T="blk_'//trim(strz(4,blockmap(2,b)))//'_grp_'//trim(strz(3,blockmap(1,b)))//'"'//&
                                        ', I='//trim(str(no_sign=.true.,n=Ni+1))//                                                &
                                        ', J='//trim(str(no_sign=.true.,n=Nj+1))//                                                &
                                        ', K='//trim(str(no_sign=.true.,n=Nk+1))//                                                &
                                        ', DATAPACKING=BLOCK'//                                                                   &
                                        ', VARLOCATION=([1-'//trim(str(.true.,nvar))//']=NODAL)'
    endif
  endif
  error = tec_var(pp=pp,n=nnode,var=node(0:Ni,0:Nj,0:Nk)%x,d=1)
  error = tec_var(pp=pp,n=nnode,var=node(0:Ni,0:Nj,0:Nk)%y,d=1)
  error = tec_var(pp=pp,n=nnode,var=node(0:Ni,0:Nj,0:Nk)%z,d=1)
  if (fcc) then
    if (cell) then
      error = tec_var(pp=pp,n=ncell,var=real(ricc(1:Ni,1:Nj,1:Nk),R8P),d=1)
    else
      error = tec_var(pp=pp,n=nnode,var=real(vicc(0:Ni,0:Nj,0:Nk),R8P),d=1)
    endif
  endif
  if (sol) then
    if (cell) then
      error = tec_var(pp=pp,n=ncell,var=momentum(1:Ni,1:Nj,1:Nk)%x,d=1)
      error = tec_var(pp=pp,n=ncell,var=momentum(1:Ni,1:Nj,1:Nk)%y,d=1)
      error = tec_var(pp=pp,n=ncell,var=momentum(1:Ni,1:Nj,1:Nk)%z,d=1)
      error = tec_var(pp=pp,n=ncell,var=pressure(1:Ni,1:Nj,1:Nk)  ,d=1)
      if (zeroeq) then
        error = tec_var(pp=pp,n=ncell,var=visc(1:Ni,1:Nj,1:Nk),d=1)
      elseif (oneeq) then
        error = tec_var(pp=pp,n=ncell,var=visc(1:Ni,1:Nj,1:Nk),d=1)
        error = tec_var(pp=pp,n=ncell,var=vitl(1:Ni,1:Nj,1:Nk),d=1)
      elseif (twoeq) then
        error = tec_var(pp=pp,n=ncell,var=visc(1:Ni,1:Nj,1:Nk),d=1)
        error = tec_var(pp=pp,n=ncell,var=ken (1:Ni,1:Nj,1:Nk),d=1)
        error = tec_var(pp=pp,n=ncell,var=eps (1:Ni,1:Nj,1:Nk),d=1)
      endif
      if (level_set) then
        error = tec_var(pp=pp,n=ncell,var=f (1:Ni,1:Nj,1:Nk),d=1)
        error = tec_var(pp=pp,n=ncell,var=f0(1:Ni,1:Nj,1:Nk),d=1)
      endif
      if (vordet) then
        error = tec_var(pp=pp,n=ncell,var=vord(1:Ni,1:Nj,1:Nk),d=1)
        error = tec_var(pp=pp,n=ncell,var=qfac(1:Ni,1:Nj,1:Nk),d=1)
        error = tec_var(pp=pp,n=ncell,var=heli(1:Ni,1:Nj,1:Nk),d=1)
        error = tec_var(pp=pp,n=ncell,var=vorticity(1:Ni,1:Nj,1:Nk)%x,d=1)
        error = tec_var(pp=pp,n=ncell,var=vorticity(1:Ni,1:Nj,1:Nk)%y,d=1)
        error = tec_var(pp=pp,n=ncell,var=vorticity(1:Ni,1:Nj,1:Nk)%z,d=1)
      endif
    else
      call varinterpolation(Ni=Ni,Nj=Nj,Nk=Nk,var=momentum(0:Ni+1,0:Nj+1,0:Nk+1)%x,vari=vari)
      error = tec_var(pp=pp,n=nnode,var=vari,d=1)
      call varinterpolation(Ni=Ni,Nj=Nj,Nk=Nk,var=momentum(0:Ni+1,0:Nj+1,0:Nk+1)%y,vari=vari)
      error = tec_var(pp=pp,n=nnode,var=vari,d=1)
      call varinterpolation(Ni=Ni,Nj=Nj,Nk=Nk,var=momentum(0:Ni+1,0:Nj+1,0:Nk+1)%z,vari=vari)
      error = tec_var(pp=pp,n=nnode,var=vari,d=1)
      call varinterpolation(Ni=Ni,Nj=Nj,Nk=Nk,var=pressure(0:Ni+1,0:Nj+1,0:Nk+1)  ,vari=vari)
      error = tec_var(pp=pp,n=nnode,var=vari,d=1)
      if (zeroeq) then
        call varinterpolation(Ni=Ni,Nj=Nj,Nk=Nk,var=visc(0:Ni+1,0:Nj+1,0:Nk+1),vari=vari)
        error = tec_var(pp=pp,n=nnode,var=vari,d=1)
      elseif (oneeq) then
        call varinterpolation(Ni=Ni,Nj=Nj,Nk=Nk,var=visc(0:Ni+1,0:Nj+1,0:Nk+1),vari=vari)
        error = tec_var(pp=pp,n=nnode,var=vari,d=1)
        call varinterpolation(Ni=Ni,Nj=Nj,Nk=Nk,var=vitl(0:Ni+1,0:Nj+1,0:Nk+1),vari=vari)
        error = tec_var(pp=pp,n=nnode,var=vari,d=1)
      elseif (twoeq) then
        call varinterpolation(Ni=Ni,Nj=Nj,Nk=Nk,var=visc(0:Ni+1,0:Nj+1,0:Nk+1),vari=vari)
        error = tec_var(pp=pp,n=nnode,var=vari,d=1)
        call varinterpolation(Ni=Ni,Nj=Nj,Nk=Nk,var=ken (0:Ni+1,0:Nj+1,0:Nk+1),vari=vari)
        error = tec_var(pp=pp,n=nnode,var=vari,d=1)
        call varinterpolation(Ni=Ni,Nj=Nj,Nk=Nk,var=eps (0:Ni+1,0:Nj+1,0:Nk+1),vari=vari)
        error = tec_var(pp=pp,n=nnode,var=vari,d=1)
      endif
      if (level_set) then
        call varinterpolation(Ni=Ni,Nj=Nj,Nk=Nk,var=f (0:Ni+1,0:Nj+1,0:Nk+1),vari=vari)
        error = tec_var(pp=pp,n=nnode,var=vari,d=1)
        call varinterpolation(Ni=Ni,Nj=Nj,Nk=Nk,var=f0(0:Ni+1,0:Nj+1,0:Nk+1),vari=vari)
        error = tec_var(pp=pp,n=nnode,var=vari,d=1)
      endif
      if (vordet) then
        call varinterpolation(Ni=Ni,Nj=Nj,Nk=Nk,var=vord(     0:Ni+1,0:Nj+1,0:Nk+1),  vari=vari)
        error = tec_var(pp=pp,n=nnode,var=vari,d=1)
        call varinterpolation(Ni=Ni,Nj=Nj,Nk=Nk,var=qfac(     0:Ni+1,0:Nj+1,0:Nk+1),  vari=vari)
        error = tec_var(pp=pp,n=nnode,var=vari,d=1)
        call varinterpolation(Ni=Ni,Nj=Nj,Nk=Nk,var=heli(     0:Ni+1,0:Nj+1,0:Nk+1),  vari=vari)
        error = tec_var(pp=pp,n=nnode,var=vari,d=1)
        call varinterpolation(Ni=Ni,Nj=Nj,Nk=Nk,var=vorticity(0:Ni+1,0:Nj+1,0:Nk+1)%x,vari=vari)
        error = tec_var(pp=pp,n=nnode,var=vari,d=1)
        call varinterpolation(Ni=Ni,Nj=Nj,Nk=Nk,var=vorticity(0:Ni+1,0:Nj+1,0:Nk+1)%y,vari=vari)
        error = tec_var(pp=pp,n=nnode,var=vari,d=1)
        call varinterpolation(Ni=Ni,Nj=Nj,Nk=Nk,var=vorticity(0:Ni+1,0:Nj+1,0:Nk+1)%z,vari=vari)
        error = tec_var(pp=pp,n=nnode,var=vari,d=1)
      endif
    endif
  endif
  endassociate
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction save_block

  !> @brief Interface function for saving variables into Tecplot file.
  function tec_var(pp,n,var,d) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_PostProcess), intent(IN):: pp       !< Post-processor data.
  integer(I_P),           intent(IN):: n        !< Number of var elements.
  real(R_P),              intent(IN):: var(1:n) !< Variable to be saved.
  integer(I_P),           intent(IN):: d        !< Double precision output (1 yes, 0 no).
  integer(I_P)::                       err      !< Error traping flag: 0 no errors, >0 error occours.
  integer(I_P)::                       e        !< Element counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (pp%binary) then
    err = tecdat112(n,var,d)
  else
    write(pp%unit_out,FR8P,iostat=err)(var(e),e=1,n)
  endif
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction tec_var

  !> @brief Procedure for finalizing file writing.
  subroutine finalize(pp)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_PostProcess), intent(IN)::    pp     !< Post-processor data.
  integer(I4P)::                          error  !< Error trapping flag.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (pp%binary) then
    error = tecend112()
  else
    close(pp%unit_out)
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine finalize
endmodule Lib_TEC
