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
  integer(I4P)::                  nvar = 4            !< Number of variables saved.
  contains
    procedure::         init       ! Procedure for initializing file writing.
    procedure::         save_patch ! Procedure for saving one patch.
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
    file_d%tecvarname = 'x y z icc'
  else
    file_d%tecvarname = ' VARIABLES ="x" "y" "z" "icc"'
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
    if (pp%forces) then
      file_d%nvar = file_d%nvar + 18 ! fx, fy, fz, mx, my and mz (total, pressure and viscous parts)
      if (pp%binary) then
        file_d%tecvarname = trim(file_d%tecvarname)//&
                            ' fx fy fz fx_p fy_p fz_p fx_v fy_v fz_v mx my mz mx_p my_p mz_p mx_v my_v mz_v'
      else
        file_d%tecvarname = trim(file_d%tecvarname)//' "fx" "fy" "fz" "fx_p" "fy_p" "fz_p" "fx_v" "fy_v" "fz_v"'//&
                                                     ' "mx" "my" "mz" "mx_p" "my_p" "mz_p" "mx_v" "my_v" "mz_v"'
      endif
    endif
    if (pp%yp) then
      file_d%nvar = file_d%nvar + 1 ! yplus
      if (pp%binary) then
        file_d%tecvarname = trim(file_d%tecvarname)//' yplus'
      else
        file_d%tecvarname = trim(file_d%tecvarname)//' "yplus"'
      endif
    endif
    if (pp%tau) then
      file_d%nvar = file_d%nvar + 3 ! taux,tauy,tauz
      if (pp%binary) then
        file_d%tecvarname = trim(file_d%tecvarname)//' taux tauy tauz'
      else
        file_d%tecvarname = trim(file_d%tecvarname)//' "taux" "tauy" "tauz"'
      endif
    endif
  endif
  if (pp%metrics) then
    file_d%nvar = file_d%nvar + 4 ! Nx,Ny,Nz and S
    if (pp%binary) then
      file_d%tecvarname = trim(file_d%tecvarname)//' Nx Ny Nz S'
    else
      file_d%tecvarname = trim(file_d%tecvarname)//' "Nx" "Ny" "Nz" "S"'
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
  function save_patch(file_d,Np,pp,b,face,ni1,ni2,nj1,nj2,nk1,nk2,ci1,ci2,cj1,cj2,ck1,ck2) &
           result(error)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_File_Tec),   intent(INOUT):: file_d      !< File data.
  integer(I4P),           intent(IN)::    Np          !< Number of patches.
  type(Type_PostProcess), intent(IN)::    pp          !< Post-processor data.
  integer(I4P),           intent(IN)::    b           !< Actual block number.
  integer(I4P),           intent(IN)::    face        !< Face where patch is defined: 1,2,3,4,5,6.
  integer(I4P),           intent(IN)::    ni1,ni2     !< First and last node i indexes.
  integer(I4P),           intent(IN)::    nj1,nj2     !< First and last node j indexes.
  integer(I4P),           intent(IN)::    nk1,nk2     !< First and last node k indexes.
  integer(I4P),           intent(IN)::    ci1,ci2     !< First and last cell i indexes.
  integer(I4P),           intent(IN)::    cj1,cj2     !< First and last cell j indexes.
  integer(I4P),           intent(IN)::    ck1,ck2     !< First and last cell k indexes.
  integer(I4P)::                          error       !< Error trapping flag.
  integer(I4P)::                          nnode,ncell !< Number of nodes and cells.
  real(R8P), allocatable::                vari(:,:,:) !< Interpolated generic variables.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.pp%cell) then ! variables must be interpolated at nodes, allocating dummy variable
    if (allocated(vari)) deallocate(vari) ; allocate(vari(ni1:ni2,nj1:nj2,nk1:nk2))
  endif
  nnode = (ni2-ni1+1)*(nj2-nj1+1)*(nk2-nk1+1) ! computing number of nodes
  ncell = (ci2-ci1+1)*(cj2-cj1+1)*(ck2-ck1+1) ! computing number of cells

  associate(binary=>pp%binary,cell=>pp%cell,blockmap=>pp%blockmap,sol=>pp%sol,yp=>pp%yp,forces=>pp%forces,patch=>pp%patch,      &
            metrics=>pp%metrics,zeroeq=>pp%zeroeq,oneeq=>pp%oneeq,twoeq=>pp%twoeq,level_set=>pp%level_set,unit_out=>pp%unit_out,&
            nvar=>file_d%nvar)
  if (binary) then
    error = teczne112('ptc_'//trim(strz(3,patch))//'_n_'//trim(strz(5,Np))// &
                      '_blk_'//trim(strz(4,blockmap(2,b)))//                 &
                      '_fac_'//trim(strz(1,face))//                          &
                      '_grp_'//trim(strz(3,blockmap(1,b)))//file_d%tecendrec,&
                      0,                                                     &
                      ni2-ni1+1,                                             &
                      nj2-nj1+1,                                             &
                      nk2-nk1+1,                                             &
                      0,                                                     &
                      0,                                                     &
                      0,                                                     &
                      0.0,                                                   &
                      0,                                                     &
                      0,                                                     &
                      1,                                                     &
                      0,                                                     &
                      0,                                                     &
                      0,                                                     &
                      0,                                                     &
                      0,                                                     &
                      file_d%tecnull,                                        &
                      file_d%tecvarloc,                                      &
                      file_d%tecnull,                                        &
                      0)
  else
    if (nvar>3.and.cell) then
      write(unit_out,'(A)',iostat=error)' ZONE  T="ptc_'//trim(strz(3,patch))//'_n_'//trim(strz(5,Np))// &
                                        '_blk_'//trim(strz(4,blockmap(2,b)))//                           &
                                        '_fac_'//trim(strz(1,face))//                                    &
                                        '_grp_'//trim(strz(3,blockmap(1,b)))//'"'//                      &
                                        ', I='//trim(str(no_sign=.true.,n=ni2-ni1+1))//                  &
                                        ', J='//trim(str(no_sign=.true.,n=nj2-nj1+1))//                  &
                                        ', K='//trim(str(no_sign=.true.,n=nk2-nk1+1))//                  &
                                        ', DATAPACKING=BLOCK'//                                          &
                                        ', VARLOCATION=([1-3]=NODAL,[4-'//trim(str(.true.,Nvar))//']=CELLCENTERED)'
    else
      write(unit_out,'(A)',iostat=error)' ZONE  T="ptc_'//trim(strz(3,patch))//'_n_'//trim(strz(5,Np))// &
                                        '_blk_'//trim(strz(4,blockmap(2,b)))//                           &
                                        '_fac_'//trim(strz(1,face))//                                    &
                                        '_grp_'//trim(strz(3,blockmap(1,b)))//'"'//                      &
                                        ', I='//trim(str(no_sign=.true.,n=ni2-ni1+1))//                  &
                                        ', J='//trim(str(no_sign=.true.,n=nj2-nj1+1))//                  &
                                        ', K='//trim(str(no_sign=.true.,n=nk2-nk1+1))//                  &
                                        ', DATAPACKING=BLOCK'//                                          &
                                        ', VARLOCATION=([1-'//trim(str(.true.,Nvar))//']=NODAL)'
    endif
  endif
  error = tec_var(pp=pp,n=nnode,var=node(ni1:ni2,nj1:nj2,nk1:nk2)%x,d=1)
  error = tec_var(pp=pp,n=nnode,var=node(ni1:ni2,nj1:nj2,nk1:nk2)%y,d=1)
  error = tec_var(pp=pp,n=nnode,var=node(ni1:ni2,nj1:nj2,nk1:nk2)%z,d=1)
  if (cell) then
    error = tec_var(pp=pp,n=ncell,var=real(ricc(ci1:ci2,cj1:cj2,ck1:ck2),R8P),d=1)
  else
    error = tec_var(pp=pp,n=nnode,var=real(vicc(ni1:ni2,nj1:nj2,nk1:nk2),R8P),d=1)
  endif
  if (sol) then
    if (cell) then
      error = tec_var(pp=pp,n=ncell,var=momentum(ci1:ci2,cj1:cj2,ck1:ck2)%x,d=1)
      error = tec_var(pp=pp,n=ncell,var=momentum(ci1:ci2,cj1:cj2,ck1:ck2)%y,d=1)
      error = tec_var(pp=pp,n=ncell,var=momentum(ci1:ci2,cj1:cj2,ck1:ck2)%z,d=1)
      error = tec_var(pp=pp,n=ncell,var=pressure(ci1:ci2,cj1:cj2,ck1:ck2)  ,d=1)
      if (zeroeq) then
        error = tec_var(pp=pp,n=ncell,var=visc(ci1:ci2,cj1:cj2,ck1:ck2),d=1)
      elseif (oneeq) then
        error = tec_var(pp=pp,n=ncell,var=visc(ci1:ci2,cj1:cj2,ck1:ck2),d=1)
        error = tec_var(pp=pp,n=ncell,var=vitl(ci1:ci2,cj1:cj2,ck1:ck2),d=1)
      elseif (twoeq) then
        error = tec_var(pp=pp,n=ncell,var=visc(ci1:ci2,cj1:cj2,ck1:ck2),d=1)
        error = tec_var(pp=pp,n=ncell,var=ken (ci1:ci2,cj1:cj2,ck1:ck2),d=1)
        error = tec_var(pp=pp,n=ncell,var=eps (ci1:ci2,cj1:cj2,ck1:ck2),d=1)
      endif
      if (level_set) then
        error = tec_var(pp=pp,n=ncell,var=f (ci1:ci2,cj1:cj2,ck1:ck2),d=1)
        error = tec_var(pp=pp,n=ncell,var=f0(ci1:ci2,cj1:cj2,ck1:ck2),d=1)
      endif
      if (forces) then
        error = tec_var(pp=pp,n=ncell,var=(f_p(ci1:ci2,cj1:cj2,ck1:ck2)%x + &
                                           f_v(ci1:ci2,cj1:cj2,ck1:ck2)%x),d=1)
        error = tec_var(pp=pp,n=ncell,var=(f_p(ci1:ci2,cj1:cj2,ck1:ck2)%y + &
                                           f_v(ci1:ci2,cj1:cj2,ck1:ck2)%y),d=1)
        error = tec_var(pp=pp,n=ncell,var=(f_p(ci1:ci2,cj1:cj2,ck1:ck2)%z + &
                                           f_v(ci1:ci2,cj1:cj2,ck1:ck2)%z),d=1)
        error = tec_var(pp=pp,n=ncell,var= f_p(ci1:ci2,cj1:cj2,ck1:ck2)%x,d=1)
        error = tec_var(pp=pp,n=ncell,var= f_p(ci1:ci2,cj1:cj2,ck1:ck2)%y,d=1)
        error = tec_var(pp=pp,n=ncell,var= f_p(ci1:ci2,cj1:cj2,ck1:ck2)%z,d=1)
        error = tec_var(pp=pp,n=ncell,var= f_v(ci1:ci2,cj1:cj2,ck1:ck2)%x,d=1)
        error = tec_var(pp=pp,n=ncell,var= f_v(ci1:ci2,cj1:cj2,ck1:ck2)%y,d=1)
        error = tec_var(pp=pp,n=ncell,var= f_v(ci1:ci2,cj1:cj2,ck1:ck2)%z,d=1)
        error = tec_var(pp=pp,n=ncell,var=(m_p(ci1:ci2,cj1:cj2,ck1:ck2)%x + &
                                           m_v(ci1:ci2,cj1:cj2,ck1:ck2)%x),d=1)
        error = tec_var(pp=pp,n=ncell,var=(m_p(ci1:ci2,cj1:cj2,ck1:ck2)%y + &
                                           m_v(ci1:ci2,cj1:cj2,ck1:ck2)%y),d=1)
        error = tec_var(pp=pp,n=ncell,var=(m_p(ci1:ci2,cj1:cj2,ck1:ck2)%z + &
                                           m_v(ci1:ci2,cj1:cj2,ck1:ck2)%z),d=1)
        error = tec_var(pp=pp,n=ncell,var= m_p(ci1:ci2,cj1:cj2,ck1:ck2)%x,d=1)
        error = tec_var(pp=pp,n=ncell,var= m_p(ci1:ci2,cj1:cj2,ck1:ck2)%y,d=1)
        error = tec_var(pp=pp,n=ncell,var= m_p(ci1:ci2,cj1:cj2,ck1:ck2)%z,d=1)
        error = tec_var(pp=pp,n=ncell,var= m_v(ci1:ci2,cj1:cj2,ck1:ck2)%x,d=1)
        error = tec_var(pp=pp,n=ncell,var= m_v(ci1:ci2,cj1:cj2,ck1:ck2)%y,d=1)
        error = tec_var(pp=pp,n=ncell,var= m_v(ci1:ci2,cj1:cj2,ck1:ck2)%z,d=1)
      endif
      if (yp) then
        error = tec_var(pp=pp,n=ncell,var=yplus(ni1:ni2,cj1:cj2,ck1:ck2),d=1)
      endif
      if (pp%tau) then
        error = tec_var(pp=pp,n=ncell,var=tau(ni1:ni2,cj1:cj2,ck1:ck2)%x,d=1)
        error = tec_var(pp=pp,n=ncell,var=tau(ni1:ni2,cj1:cj2,ck1:ck2)%y,d=1)
        error = tec_var(pp=pp,n=ncell,var=tau(ni1:ni2,cj1:cj2,ck1:ck2)%z,d=1)
      endif
    else
      call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                            var  = momentum(ni1:ni2+1,nj1:nj2+1,nk1:nk2+1)%x,         &
                            vari = vari(    ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
      error = tec_var(pp=pp,n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
      call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                            var  = momentum(ni1:ni2+1,nj1:nj2+1,nk1:nk2+1)%y,         &
                            vari = vari(    ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
      error = tec_var(pp=pp,n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
      call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                            var  = momentum(ni1:ni2+1,nj1:nj2+1,nk1:nk2+1)%z,         &
                            vari = vari(    ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
      error = tec_var(pp=pp,n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
      call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                            var  = pressure(ni1:ni2+1,nj1:nj2+1,nk1:nk2+1),           &
                            vari = vari(    ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
      error = tec_var(pp=pp,n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
      if (zeroeq) then
        call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                              var  = visc(ni1:ni2+1,nj1:nj2+1,nk1:nk2+1),               &
                              vari = vari(ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
        error = tec_var(pp=pp,n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
      elseif (oneeq) then
        call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                              var  = visc(ni1:ni2+1,nj1:nj2+1,nk1:nk2+1),               &
                              vari = vari(ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
        error = tec_var(pp=pp,n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
        call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                              var  = vitl(ni1:ni2+1,nj1:nj2+1,nk1:nk2+1),               &
                              vari = vari(ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
        error = tec_var(pp=pp,n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
      elseif (twoeq) then
        call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                              var  = visc(ni1:ni2+1,nj1:nj2+1,nk1:nk2+1),               &
                              vari = vari(ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
        error = tec_var(pp=pp,n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
        call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                              var  = ken( ni1:ni2+1,nj1:nj2+1,nk1:nk2+1),               &
                              vari = vari(ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
        error = tec_var(pp=pp,n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
        call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                              var  = eps( ni1:ni2+1,nj1:nj2+1,nk1:nk2+1),               &
                              vari = vari(ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
        error = tec_var(pp=pp,n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
      endif
      if (level_set) then
        call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                              var  = f(   ni1:ni2+1,nj1:nj2+1,nk1:nk2+1),               &
                              vari = vari(ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
        error = tec_var(pp=pp,n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
        call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                              var  = f0(  ni1:ni2+1,nj1:nj2+1,nk1:nk2+1),               &
                              vari = vari(ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
        error = tec_var(pp=pp,n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
      endif
      if (forces) then
        call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                              var  = f_p (ni1:ni2+1,nj1:nj2+1,nk1:nk2+1)%x +            &
                                     f_v (ni1:ni2+1,nj1:nj2+1,nk1:nk2+1)%x,             &
                              vari = vari(ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
        error = tec_var(pp=pp,n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
        call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                              var  = f_p (ni1:ni2+1,nj1:nj2+1,nk1:nk2+1)%y +            &
                                     f_v (ni1:ni2+1,nj1:nj2+1,nk1:nk2+1)%y,             &
                              vari = vari(ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
        error = tec_var(pp=pp,n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
        call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                              var  = f_p (ni1:ni2+1,nj1:nj2+1,nk1:nk2+1)%z +            &
                                     f_v (ni1:ni2+1,nj1:nj2+1,nk1:nk2+1)%z,             &
                              vari = vari(ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
        error = tec_var(pp=pp,n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
        call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                              var  = f_p (ni1:ni2+1,nj1:nj2+1,nk1:nk2+1)%x,             &
                              vari = vari(ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
        error = tec_var(pp=pp,n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
        call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                              var  = f_p (ni1:ni2+1,nj1:nj2+1,nk1:nk2+1)%y,             &
                              vari = vari(ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
        error = tec_var(pp=pp,n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
        call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                              var  = f_p (ni1:ni2+1,nj1:nj2+1,nk1:nk2+1)%z,             &
                              vari = vari(ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
        error = tec_var(pp=pp,n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
        call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                              var  = f_v (ni1:ni2+1,nj1:nj2+1,nk1:nk2+1)%x,             &
                              vari = vari(ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
        error = tec_var(pp=pp,n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
        call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                              var  = f_v (ni1:ni2+1,nj1:nj2+1,nk1:nk2+1)%y,             &
                              vari = vari(ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
        error = tec_var(pp=pp,n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
        call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                              var  = f_v (ni1:ni2+1,nj1:nj2+1,nk1:nk2+1)%z,             &
                              vari = vari(ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
        error = tec_var(pp=pp,n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
        call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                              var  = m_p (ni1:ni2+1,nj1:nj2+1,nk1:nk2+1)%x +            &
                                     m_v (ni1:ni2+1,nj1:nj2+1,nk1:nk2+1)%x,             &
                              vari = vari(ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
        error = tec_var(pp=pp,n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
        call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                              var  = m_p (ni1:ni2+1,nj1:nj2+1,nk1:nk2+1)%y +            &
                                     m_v (ni1:ni2+1,nj1:nj2+1,nk1:nk2+1)%y,             &
                              vari = vari(ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
        error = tec_var(pp=pp,n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
        call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                              var  = m_p (ni1:ni2+1,nj1:nj2+1,nk1:nk2+1)%z +            &
                                     m_v (ni1:ni2+1,nj1:nj2+1,nk1:nk2+1)%z,             &
                              vari = vari(ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
        error = tec_var(pp=pp,n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
        call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                              var  = m_p (ni1:ni2+1,nj1:nj2+1,nk1:nk2+1)%x,             &
                              vari = vari(ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
        error = tec_var(pp=pp,n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
        call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                              var  = m_p (ni1:ni2+1,nj1:nj2+1,nk1:nk2+1)%y,             &
                              vari = vari(ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
        error = tec_var(pp=pp,n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
        call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                              var  = m_p (ni1:ni2+1,nj1:nj2+1,nk1:nk2+1)%z,             &
                              vari = vari(ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
        error = tec_var(pp=pp,n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
        call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                              var  = m_v (ni1:ni2+1,nj1:nj2+1,nk1:nk2+1)%x,             &
                              vari = vari(ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
        error = tec_var(pp=pp,n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
        call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                              var  = m_v (ni1:ni2+1,nj1:nj2+1,nk1:nk2+1)%y,             &
                              vari = vari(ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
        error = tec_var(pp=pp,n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
        call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                              var  = m_v (ni1:ni2+1,nj1:nj2+1,nk1:nk2+1)%z,             &
                              vari = vari(ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
        error = tec_var(pp=pp,n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
      endif
      if (yp) then
        call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                              var  = yplus(ni1:ni2+1,nj1:nj2+1,nk1:nk2+1),              &
                              vari = vari( ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
        error = tec_var(pp=pp,n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
      endif
      if (pp%tau) then
        call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                              var  = tau( ni1:ni2+1,nj1:nj2+1,nk1:nk2+1)%x,             &
                              vari = vari(ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
        error = tec_var(pp=pp,n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
        call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                              var  = tau( ni1:ni2+1,nj1:nj2+1,nk1:nk2+1)%y,             &
                              vari = vari(ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
        error = tec_var(pp=pp,n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
        call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                              var  = tau( ni1:ni2+1,nj1:nj2+1,nk1:nk2+1)%z,             &
                              vari = vari(ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
        error = tec_var(pp=pp,n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
      endif
    endif
  endif
  if (metrics) then
    if (cell) then
      select case(face)
      case(1,2)
        error = tec_var(pp=pp,n=ncell,var=NFi(ni1,cj1:cj2,ck1:ck2)%x,d=1)
        error = tec_var(pp=pp,n=ncell,var=NFi(ni1,cj1:cj2,ck1:ck2)%y,d=1)
        error = tec_var(pp=pp,n=ncell,var=NFi(ni1,cj1:cj2,ck1:ck2)%z,d=1)
        error = tec_var(pp=pp,n=ncell,var= Si(ni1,cj1:cj2,ck1:ck2)  ,d=1)
      case(3,4)
        error = tec_var(pp=pp,n=ncell,var=NFj(ci1:ci2,nj1,ck1:ck2)%x,d=1)
        error = tec_var(pp=pp,n=ncell,var=NFj(ci1:ci2,nj1,ck1:ck2)%y,d=1)
        error = tec_var(pp=pp,n=ncell,var=NFj(ci1:ci2,nj1,ck1:ck2)%z,d=1)
        error = tec_var(pp=pp,n=ncell,var= Sj(ci1:ci2,nj1,ck1:ck2)  ,d=1)
      case(5,6)
        error = tec_var(pp=pp,n=ncell,var=NFk(ci1:ci2,cj1:cj2,nk1)%x,d=1)
        error = tec_var(pp=pp,n=ncell,var=NFk(ci1:ci2,cj1:cj2,nk1)%y,d=1)
        error = tec_var(pp=pp,n=ncell,var=NFk(ci1:ci2,cj1:cj2,nk1)%z,d=1)
        error = tec_var(pp=pp,n=ncell,var= Sk(ci1:ci2,cj1:cj2,nk1)  ,d=1)
      endselect
    else
      select case(face)
      case(1,2)
        call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                              var  =  NFi(ni1:ni2+1,nj1:nj2+1,nk1:nk2+1)%x,             &
                              vari = vari(ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
        error = tec_var(pp=pp,n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
        call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                              var  =  NFi(ni1:ni2+1,nj1:nj2+1,nk1:nk2+1)%y,             &
                              vari = vari(ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
        error = tec_var(pp=pp,n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
        call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                              var  =  NFi(ni1:ni2+1,nj1:nj2+1,nk1:nk2+1)%z,             &
                              vari = vari(ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
        error = tec_var(pp=pp,n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
        call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                              var  =   Si(ni1:ni2+1,nj1:nj2+1,nk1:nk2+1),               &
                              vari = vari(ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
        error = tec_var(pp=pp,n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
      case(3,4)
        call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                              var  =  NFj(ni1:ni2+1,nj1:nj2+1,nk1:nk2+1)%x,             &
                              vari = vari(ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
        error = tec_var(pp=pp,n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
        call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                              var  =  NFj(ni1:ni2+1,nj1:nj2+1,nk1:nk2+1)%y,             &
                              vari = vari(ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
        error = tec_var(pp=pp,n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
        call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                              var  =  NFj(ni1:ni2+1,nj1:nj2+1,nk1:nk2+1)%z,             &
                              vari = vari(ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
        error = tec_var(pp=pp,n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
        call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                              var  =   Sj(ni1:ni2+1,nj1:nj2+1,nk1:nk2+1),               &
                              vari = vari(ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
        error = tec_var(pp=pp,n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
      case(5,6)
        call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                              var  =  NFk(ni1:ni2+1,nj1:nj2+1,nk1:nk2+1)%x,             &
                              vari = vari(ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
        error = tec_var(pp=pp,n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
        call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                              var  =  NFk(ni1:ni2+1,nj1:nj2+1,nk1:nk2+1)%y,             &
                              vari = vari(ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
        error = tec_var(pp=pp,n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
        call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                              var  =  NFk(ni1:ni2+1,nj1:nj2+1,nk1:nk2+1)%z,             &
                              vari = vari(ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
        error = tec_var(pp=pp,n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
        call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                              var  =   Sk(ni1:ni2+1,nj1:nj2+1,nk1:nk2+1),               &
                              vari = vari(ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
        error = tec_var(pp=pp,n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
      endselect
    endif
  endif
  endassociate
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction save_patch

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
