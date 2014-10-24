program XnPatches
!-----------------------------------------------------------------------------------------------------------------------------------
! Code Name XnPatches
! Author    Stefano Zaghi
! Version   0.0.1
! Date      14-02-2012
!
! XnPatches loads Xnavis mesh and solution files and produces post-processed Tecplot file of selected patches.
!
! TODO: Graphical User Interface
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision                                                                    ! Integers and reals precision definition.
USE Data_Type_Vector                                                                ! Definition of Type_Vector.
USE Data_Type_OS                                                                    ! Definition of Type_OS.
USE Lib_IO_Misc                                                                     ! Library for miscellanea IO procedures.
USE Lib_VTK_IO                                                                      ! Library for I/O VTK files.
USE, intrinsic:: ISO_FORTRAN_ENV, only: stdout => OUTPUT_UNIT, stderr => ERROR_UNIT ! Standard output/error logical units.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
character(100)::            File_grd ='unset'          ! GRD file name.
character(100)::            File_icc ='unset'          ! ICC file name.
character(100)::            File_sol ='unset'          ! Solution file name.
character(100)::            File_out ='unset'          ! Output file name.
integer(I_P)::              unit_grd                   ! Unit of GRD file.
integer(I_P)::              unit_icc                   ! Unit of ICC file.
integer(I_P)::              unit_sol                   ! Unit of Solution file.
integer(I_P)::              unit_out                   ! Unit of Output file.
integer(I_P)::              unit_for                   ! Unit of forces.dat file.
integer(I_P), allocatable:: unit_g_for(:)              ! Unit of forces.dat file for each group.
logical::                   is_file                    ! Flag for inquiring the presence of file.
integer(I_P)::              patch = 1_I_P              ! Boundary condition value of post-processed patches.
integer(I_P)::              ni1,ni2                    ! First and last node i indexes.
integer(I_P)::              nj1,nj2                    ! First and last node j indexes.
integer(I_P)::              nk1,nk2                    ! First and last node k indexes.
integer(I_P)::              ci1,ci2                    ! First and last cell i indexes.
integer(I_P)::              cj1,cj2                    ! First and last cell j indexes.
integer(I_P)::              ck1,ck2                    ! First and last cell k indexes.

integer(I_P)::              istrml                     ! streamlines

logical::                   patch_found                ! Inquiring flag for patches searching.
integer(I_P)::              i,j,k,b,v                  ! Counters.
integer(I_P)::              err                        ! Error traping flag: 0 no errors, >0 error occours.
integer(I_P)::              Np = 0_I_P                 ! Number of patches.
integer(I_P)::              nprocs = 0_I_P             ! number of processes
integer(I_P)::              myrank = 0_I_P             ! Rank of current processed file.
integer(I_P), allocatable:: blockmap(:,:),procmap(:,:) ! Blocks maps.
logical::                   global_blk_num = .false.   ! Flag for inquiring the activation of global blocks numeration.
integer(I_P)::              Nbtot = 0_I_P              ! Number of total blocks.
! mesh variables
integer(I_P)::                   Nb = 0_I_P        ! Number of blocks.
integer(I_P)::                   Ng = 0_I_P        ! Number of groups.
integer(I_P),      allocatable:: gc(:,:)           ! Number of ghost cells.
integer(I_P),      allocatable:: Ni(:),Nj(:),Nk(:) ! Number of cells of each block.
type(Type_Vector), allocatable:: node(:,:,:)       ! Nodes coordinates.
! icc variables
integer(I_P),      allocatable:: icc(:,:,:)   ! Cell centered icc values.
integer(I_P),      allocatable:: ricc(:,:,:)  ! Cell centered rcc values.
integer(I_P),      allocatable:: ticc(:,:,:)  ! Cell centered rcc values (only type of cc).
integer(I_P),      allocatable:: vicc(:,:,:)  ! Node centered rcc values.
real(R4P),         allocatable:: rcc(:)       ! rcc array.
integer(I_P)::                   Nrcc = 0_I_P ! rcc array dimension.
! metrics variables
type(Type_Vector), allocatable:: NFi(:,:,:)        ! Face i normal versor.
type(Type_Vector), allocatable:: NFj(:,:,:)        ! Face j normal versor.
type(Type_Vector), allocatable:: NFk(:,:,:)        ! Face k normal versor.
type(Type_Vector), allocatable:: NFiS(:,:,:)       ! Face i normal versor with surface area module.
type(Type_Vector), allocatable:: NFjS(:,:,:)       ! Face j normal versor with surface area module.
type(Type_Vector), allocatable:: NFkS(:,:,:)       ! Face k normal versor with surface area module.
real(R_P),         allocatable:: Si(:,:,:)         ! Face i area.
real(R_P),         allocatable:: Sj(:,:,:)         ! Face j area.
real(R_P),         allocatable:: Sk(:,:,:)         ! Face k area.
real(R_P),         allocatable:: volume(:,:,:)     ! Volumes of cells.
! solution variables
integer(I_P)::                   Nvar = 0_I_P    ! Number of variables
type(Type_Vector), allocatable:: momentum(:,:,:) ! Momentum.
real(R_P),         allocatable:: pressure(:,:,:) ! Pressure.
real(R_P),         allocatable:: f(:,:,:)        ! Level set function.
real(R_P),         allocatable:: f0(:,:,:)       ! Level 0 (level set).
real(R_P),         allocatable:: visc(:,:,:)     ! Viscosity.
real(R_P),         allocatable:: vitl(:,:,:)     ! Turbulent viscosity.
real(R_P),         allocatable:: ken(:,:,:)      ! Turbulent kinetic energy.
real(R_P),         allocatable:: eps(:,:,:)      ! Turbulent kinetic energy dissipation.
real(R_P),         allocatable:: yplus(:,:,:)    ! Estimation of y+.
type(Type_Vector), allocatable:: f_p(:,:,:)      ! Forces vector, pressure part.
type(Type_Vector), allocatable:: f_v(:,:,:)      ! Forces vector, viscous  part.
type(Type_Vector), allocatable:: m_p(:,:,:)      ! Moments vector, pressure part.
type(Type_Vector), allocatable:: m_v(:,:,:)      ! Moments vector, viscous  part.
type(Type_Vector)::              fsum_p          ! Global sum of forces vector, pressure part.
type(Type_Vector)::              fsum_v          ! Global sum of forces vector, viscous part.
type(Type_Vector)::              msum_p          ! Global sum of moments vector, pressure part.
type(Type_Vector)::              msum_v          ! Global sum of moments vector, viscous part.
real(R_P)::                      Ssum = 0._R_P   ! Global sum of "wet" surface.
real(R_P)::                      Re = -1._R_P    ! Reynolds number.
real(R_P)::                      Fr = -1._R_P    ! Froude number.
real(R_P)::                      rFr2 = 0._R_P   ! 1/(Froude number)^2.
real(R_P)::                      zfs = MaxR_P    ! Z quote of free surface.
real(R_P)::                      t               ! Solution time.
! misc variables
logical::                 sol       = .false. ! Inquiring flag for patches with solution (x,y,z,icc,u,v...).
logical::                 cell      = .false. ! Inquiring flag for interpolation (or not) variables at nodes.
logical::                 level_set = .false. ! Inquiring flag for level set variable.
logical::                 laminar   = .false. ! Inquiring flag for no turbulent model.
logical::                 zeroeq    = .false. ! Inquiring flag for zero equations turbulent variables.
logical::                 oneeq     = .true.  ! Inquiring flag for one  equations turbulent variables.
logical::                 twoeq     = .false. ! Inquiring flag for two  equations turbulent variables.
logical::                 forces    = .false. ! Inquiring flag for forces computing.
logical::                 yp        = .false. ! Inquiring flag for yplus computing.
logical::                 metrics   = .false. ! Inquiring flag for metrics saving.
logical::                 binary    = .true.  ! Inquiring flag for binary output file.
logical::                 tec       = .true.  ! Inquiring flag for tecplot file format.
logical::                 vtk       = .false. ! Inquiring flag for vtk file format.
integer(I_P), external::  tecini112,    &     ! |
                          tecauxstr112, &     ! |
                          teczne112,    &     ! | Tecplot external binary functions.
                          tecdat112,    &     ! |
                          tecend112           ! |
character(1), parameter:: tecendrec = char(0) ! End-character for binary-record finalize.
character(800)::          tecvarname          ! Variables name for tecplot header file.
character(100)::          vtkbdir             ! VTK base name of output directory.
character(100)::          vtkbfile            ! VTK base name of output files.
! Riccardo Broglia vars
integer(I_P)::           unit_for_RB     ! Unit of forces.RB (Riccardo Broglia).
integer(I_P)::           unit_for_RB_scr ! Unit of forces.RB (Riccardo Broglia) scratch file.
real(R_P), allocatable:: dummy(:,:,:)    ! Dummy var for reading scratch file forces.RB Riccardo Broglia.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
! parsing command line arguments
call parse_command_arguments

! opening input files
inquire(file=adjustl(trim(File_grd)),exist=is_file,iostat=err)
if (.NOT.is_file) then
  call File_Not_Found(0,File_grd,'XnPatches')
else
  unit_grd = Get_Unit() ; open(unit=unit_grd,file=adjustl(trim(File_grd)),form='UNFORMATTED',action='READ')
endif
inquire(file=adjustl(trim(File_icc)),exist=is_file,iostat=err)
if (.NOT.is_file) then
  call File_Not_Found(0,File_icc,'XnPatches')
else
  unit_icc = Get_Unit() ; open(unit=unit_icc,file=adjustl(trim(File_icc)),form='UNFORMATTED',action='READ')
endif
if (sol) then
  inquire(file=adjustl(trim(File_sol)),exist=is_file,iostat=err)
  if (.NOT.is_file) then
    call File_Not_Found(0,File_sol,'XnPatches')
  else
    unit_sol = Get_Unit() ; open(unit=unit_sol,file=adjustl(trim(File_sol)),form='UNFORMATTED',action='READ')
  endif
endif

! opening output file
if (tec) then
  if (binary) then
    err = tecini112('PATCH_'//trim(strz(2,patch))//tecendrec, &
                    trim(tecvarname)//tecendrec,              &
                    adjustl(trim(File_out))//"."//trim(strz(2,patch))//".plt"//tecendrec,'.'//tecendrec,0,0,1)
  else
    unit_out = Get_Unit() ; open(unit=unit_out,file=adjustl(trim(File_out))//"."//trim(strz(2,patch))//".dat")
    write(unit_out,'(A)',iostat=err)' TITLES ="PATCH_'//trim(strz(2,patch))//'"'
    write(unit_out,'(A)',iostat=err)trim(tecvarname)
  endif
endif
if (vtk) then
  ! making VTS_files directory
  err = make_dir(directory=vtkbdir)
endif
if (forces) then
  unit_for = Get_Unit() ; open(unit=unit_for,file=adjustl(trim(File_out))//"-forces.dat")
  unit_for_RB = Get_Unit() ; open(unit=unit_for_RB,file=adjustl(trim(File_out))//"-forces.RB",form='UNFORMATTED')
  unit_for_RB_scr = Get_Unit() ; open(unit=unit_for_RB_scr,form='UNFORMATTED',status='SCRATCH')
  if (Ng>0) then
    if (allocated(unit_g_for)) deallocate(unit_g_for) ; allocate(unit_g_for(0:Ng))
    do i=0,Ng
      unit_g_for(i) = Get_Unit() ; open(unit=unit_g_for(i),file=adjustl(trim(File_out))//"-forces-grp_"//trim(strz(3,i))//".dat")
    enddo
  endif
  ! initializing global sum of forces and moments
  fsum_p = 0._R_P
  fsum_v = 0._R_P
  msum_p = 0._R_P
  msum_v = 0._R_P
endif

! loading mesh dimensions
read(unit_grd)Nb
if (.not.global_blk_num) then
  ! there is no global blocks numerations: the blocks map is identity
  if (allocated(blockmap)) deallocate(blockmap) ; allocate(blockmap(1:2,1:Nb)) ; blockmap = 0_I_P
  do b=1,Nb
    blockmap(2,b) = b
  enddo
endif
read(unit_icc)b
if (b/=Nb) then
  write(stderr,'(A)')' Inconsistent ICC file with GRD one:'
  write(stderr,'(A)')' Nb GRD file: '//trim(str(.true.,Nb))
  write(stderr,'(A)')' Nb ICC file: '//trim(str(.true.,b))
  stop
endif
if (sol) then
  read(unit_sol)t
endif
write(stdout,'(A)')' Number of block of input files: '//trim(str(.true.,Nb))
if (allocated(gc)) deallocate(gc) ; allocate(gc(1:6,1:Nb)) ; gc = 2
if (allocated(Ni)) deallocate(Ni) ; allocate(Ni(    1:Nb))
if (allocated(Nj)) deallocate(Nj) ; allocate(Nj(    1:Nb))
if (allocated(Nk)) deallocate(Nk) ; allocate(Nk(    1:Nb))
write(stdout,'(A)')' Blocks dimensions'
do b=1,Nb
  read(unit_grd,err=1)Ni(b),Nj(b),Nk(b),gc(1,b),gc(2,b),gc(3,b),gc(4,b),gc(5,b),gc(6,b)
  1 continue
  write(stdout,'(A)')'   b,Ni,Nk,Nk,gc(6): '//trim(str(.true.,b))//', '//       &
                                              trim(str(.true.,Ni(b)))//', '//   &
                                              trim(str(.true.,Nj(b)))//', '//   &
                                              trim(str(.true.,Nk(b)))//', '//   &
                                              trim(str(.true.,gc(1,b)))//', '// &
                                              trim(str(.true.,gc(2,b)))//', '// &
                                              trim(str(.true.,gc(3,b)))//', '// &
                                              trim(str(.true.,gc(4,b)))//', '// &
                                              trim(str(.true.,gc(5,b)))//', '// &
                                              trim(str(.true.,gc(6,b)))
  read(unit_icc,err=2)i,j,k,v,v,v,v,v,v
  2 continue
  if (i/=Ni(b).OR.j/=Nj(b).OR.k/=Nk(b)) then
    write(stderr,'(A)')' Inconsistent ICC file with GRD one:'
    write(stderr,'(A)')' Block: '//trim(str(.true.,b))
    write(stderr,'(A)')' Ni,Nj,Nk GRD file: '//trim(str(.true.,Ni(b)))//','//trim(str(.true.,Nj(b)))//','//trim(str(.true.,Nk(b)))
    write(stderr,'(A)')' Ni,Nj,Nk ICC file: '//trim(str(.true.,i))//','//trim(str(.true.,j))//','//trim(str(.true.,k))
    stop
  endif
  if (sol) then
    read(unit_sol)i,j,k
    if (i/=Ni(b).OR.j/=Nj(b).OR.k/=Nk(b)) then
      write(stderr,'(A)')' Inconsistent SOL file with GRD one:'
      write(stderr,'(A)')' Block: '//trim(str(.true.,b))
      write(stderr,'(A)')' Ni,Nj,Nk GRD file: '//trim(str(.true.,Ni(b)))//','//trim(str(.true.,Nj(b)))//','//trim(str(.true.,Nk(b)))
      write(stderr,'(A)')' Ni,Nj,Nk SOL file: '//trim(str(.true.,i))//','//trim(str(.true.,j))//','//trim(str(.true.,k))
      stop
    endif
  endif
enddo
do b=1,Nb
  read(unit_icc)(((v,i=1-gc(1,b),Ni(b)+gc(2,b)),j=1-gc(3,b),Nj(b)+gc(4,b)),k=1-gc(5,b),Nk(b)+gc(6,b))
enddo
read(unit_icc)Nrcc
if (allocated(rcc)) deallocate(rcc) ; allocate(rcc(1:Nrcc))
read(unit_icc)(rcc(v),v=1,Nrcc)
! rewind icc file at icc pointer record
rewind(unit_icc)
read(unit_icc)b
do b=1,Nb
  read(unit_icc)i,j,k
enddo

Np = 0 ! initialize patch number
do b=1,Nb ! post processing patches
  ! allocating mesh variables
  if (allocated(node))   deallocate(node)
  allocate(node(0-gc(1,b):Ni(b)+gc(2,b),0-gc(3,b):Nj(b)+gc(4,b),0-gc(5,b):Nk(b)+gc(6,b)))
  ! loading block nodes coordinates
  read(unit_grd)(((node(i,j,k)%x,i=0-gc(1,b),Ni(b)+gc(2,b)),j=0-gc(3,b),Nj(b)+gc(4,b)),k=0-gc(5,b),Nk(b)+gc(6,b))
  read(unit_grd)(((node(i,j,k)%y,i=0-gc(1,b),Ni(b)+gc(2,b)),j=0-gc(3,b),Nj(b)+gc(4,b)),k=0-gc(5,b),Nk(b)+gc(6,b))
  read(unit_grd)(((node(i,j,k)%z,i=0-gc(1,b),Ni(b)+gc(2,b)),j=0-gc(3,b),Nj(b)+gc(4,b)),k=0-gc(5,b),Nk(b)+gc(6,b))
  ! allocating icc variables
  if (allocated(icc))  deallocate(icc)
  if (allocated(ricc)) deallocate(ricc)
  if (allocated(ticc)) deallocate(ticc)
  if (allocated(vicc)) deallocate(vicc)
                 allocate(icc (1-gc(1,b):Ni(b)+gc(2,b),1-gc(3,b):Nj(b)+gc(4,b),1-gc(5,b):Nk(b)+gc(6,b)))
                 allocate(ricc(0        :Ni(b)+1      ,0        :Nj(b)+1      ,0        :Nk(b)+1      ))
                 allocate(ticc(0        :Ni(b)+1      ,0        :Nj(b)+1      ,0        :Nk(b)+1      ))
  if (.not.cell) allocate(vicc(0        :Ni(b)+1      ,0        :Nj(b)+1      ,0        :Nk(b)+1      ))
  ! loading icc variables
  read(unit_icc)(((icc(i,j,k),i=1-gc(1,b),Ni(b)+gc(2,b)),j=1-gc(3,b),Nj(b)+gc(4,b)),k=1-gc(5,b),Nk(b)+gc(6,b))
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(i,j,k)         &
  !$OMP SHARED(b,Ni,Nj,Nk,icc,rcc,ticc,ricc)
  ! extracting type of bc
  !$OMP DO
  do k=0,Nk(b)+1
    do j=0,Nj(b)+1
      do i=0,Ni(b)+1
        ticc(i,j,k) = 0
        if (icc(i,j,k)>0) ticc(i,j,k) = abs(nint(rcc(icc(i,j,k))))
      enddo
    enddo
  enddo
  ! computing ricc for blanking
  !$OMP DO
  do k=1,Nk(b)
    do j=1,Nj(b)
      do i=1,Ni(b)
        if (icc(i,j,k)>0) then
          ricc(i,j,k) = nint(rcc(icc(i,j,k)))
        else
          ricc(i,j,k) = 0
        end if
        if (ricc(i,j,k)>28) ricc(i,j,k) = 0
      enddo
    enddo
  enddo
  !$OMP END PARALLEL
  if (.not.cell) then ! ricc must be interpolated at nodes
    ricc(0      ,:      ,:      ) = ricc(1    ,:    ,:    )
    ricc(Ni(b)+1,:      ,:      ) = ricc(Ni(b),:    ,:    )
    ricc(:      ,0      ,:      ) = ricc(:    ,1    ,:    )
    ricc(:      ,Nj(b)+1,:      ) = ricc(:    ,Nj(b),:    )
    ricc(:      ,:      ,0      ) = ricc(:    ,:    ,1    )
    ricc(:      ,:      ,Nk(b)+1) = ricc(:    ,:    ,Nk(b))
    do k=0,Nk(b)
      do j=0,Nj(b)
        do i=0,Ni(b)
          vicc(i,j,k) = max(0,maxval(ricc(i:i+1,j:j+1,k:k+1)))
        enddo
      enddo
    enddo
  endif
  if (metrics.or.forces.or.yp) then ! allocating and computing metrics variables
    ! allocating metrics variables
    if (allocated(NFi))    deallocate(NFi)
    if (allocated(NFj))    deallocate(NFj)
    if (allocated(NFk))    deallocate(NFk)
    if (allocated(NFiS))   deallocate(NFiS)
    if (allocated(NFjS))   deallocate(NFjS)
    if (allocated(NFkS))   deallocate(NFkS)
    if (allocated(Si))     deallocate(Si)
    if (allocated(Sj))     deallocate(Sj)
    if (allocated(Sk))     deallocate(Sk)
    if (allocated(volume)) deallocate(volume)
    allocate(NFi   (0-gc(1,b):Ni(b)+gc(2,b),0-gc(3,b):Nj(b)+gc(4,b),0-gc(5,b):Nk(b)+gc(6,b)))
    allocate(NFj   (0-gc(1,b):Ni(b)+gc(2,b),0-gc(3,b):Nj(b)+gc(4,b),0-gc(5,b):Nk(b)+gc(6,b)))
    allocate(NFk   (0-gc(1,b):Ni(b)+gc(2,b),0-gc(3,b):Nj(b)+gc(4,b),0-gc(5,b):Nk(b)+gc(6,b)))
    allocate(NFiS  (0-gc(1,b):Ni(b)+gc(2,b),0-gc(3,b):Nj(b)+gc(4,b),0-gc(5,b):Nk(b)+gc(6,b)))
    allocate(NFjS  (0-gc(1,b):Ni(b)+gc(2,b),0-gc(3,b):Nj(b)+gc(4,b),0-gc(5,b):Nk(b)+gc(6,b)))
    allocate(NFkS  (0-gc(1,b):Ni(b)+gc(2,b),0-gc(3,b):Nj(b)+gc(4,b),0-gc(5,b):Nk(b)+gc(6,b)))
    allocate(Si    (0-gc(1,b):Ni(b)+gc(2,b),0-gc(3,b):Nj(b)+gc(4,b),0-gc(5,b):Nk(b)+gc(6,b)))
    allocate(Sj    (0-gc(1,b):Ni(b)+gc(2,b),0-gc(3,b):Nj(b)+gc(4,b),0-gc(5,b):Nk(b)+gc(6,b)))
    allocate(Sk    (0-gc(1,b):Ni(b)+gc(2,b),0-gc(3,b):Nj(b)+gc(4,b),0-gc(5,b):Nk(b)+gc(6,b)))
    allocate(volume(1-gc(1,b):Ni(b)+gc(2,b),1-gc(3,b):Nj(b)+gc(4,b),1-gc(5,b):Nk(b)+gc(6,b)))
    ! computing the metrics of cells
    call compute_metrics(b=b)
    ! correcting the metrics of boundary conditions cells
    call bc_metrics_correction(b=b)
  endif
  if (sol) then ! allocating and loading solution variables
    if (allocated(momentum)) deallocate(momentum)
    if (allocated(pressure)) deallocate(pressure)
    allocate(momentum(1-gc(1,b):Ni(b)+gc(2,b),1-gc(3,b):Nj(b)+gc(4,b),1-gc(5,b):Nk(b)+gc(6,b)))
    allocate(pressure(1-gc(1,b):Ni(b)+gc(2,b),1-gc(3,b):Nj(b)+gc(4,b),1-gc(5,b):Nk(b)+gc(6,b)))
    read(unit_sol)(((momentum(i,j,k)%x,i=0,Ni(b)+1),j=0,Nj(b)+1),k=0,Nk(b)+1)
    read(unit_sol)(((momentum(i,j,k)%y,i=0,Ni(b)+1),j=0,Nj(b)+1),k=0,Nk(b)+1)
    read(unit_sol)(((momentum(i,j,k)%z,i=0,Ni(b)+1),j=0,Nj(b)+1),k=0,Nk(b)+1)
    read(unit_sol)(((pressure(i,j,k)  ,i=0,Ni(b)+1),j=0,Nj(b)+1),k=0,Nk(b)+1)
    if (level_set) then
      if (allocated(f))  deallocate(f)  ; allocate(f (1-gc(1,b):Ni(b)+gc(2,b),1-gc(3,b):Nj(b)+gc(4,b),1-gc(5,b):Nk(b)+gc(6,b)))
      if (allocated(f0)) deallocate(f0) ; allocate(f0(1-gc(1,b):Ni(b)+gc(2,b),1-gc(3,b):Nj(b)+gc(4,b),1-gc(5,b):Nk(b)+gc(6,b)))
      read(unit_sol)(((f (i,j,k),i=0,Ni(b)+1),j=0,Nj(b)+1),k=0,Nk(b)+1)
      read(unit_sol)(((f0(i,j,k),i=0,Ni(b)+1),j=0,Nj(b)+1),k=0,Nk(b)+1)
    endif
    if (zeroeq.or.oneeq.or.twoeq) then
      if (allocated(visc)) deallocate(visc)
      allocate(visc(1-gc(1,b):Ni(b)+gc(2,b),1-gc(3,b):Nj(b)+gc(4,b),1-gc(5,b):Nk(b)+gc(6,b)))
      read(unit_sol)(((visc(i,j,k),i=0,Ni(b)+1),j=0,Nj(b)+1),k=0,Nk(b)+1)
    endif
    if (oneeq) then
      if (allocated(vitl)) deallocate(vitl)
      allocate(vitl(1-gc(1,b):Ni(b)+gc(2,b),1-gc(3,b):Nj(b)+gc(4,b),1-gc(5,b):Nk(b)+gc(6,b)))
      read(unit_sol)(((vitl(i,j,k),i=0,Ni(b)+1),j=0,Nj(b)+1),k=0,Nk(b)+1)
    endif
    if (twoeq) then
      if (allocated(ken)) deallocate(ken) ; allocate(ken(1-gc(1,b):Ni(b)+gc(2,b),1-gc(3,b):Nj(b)+gc(4,b),1-gc(5,b):Nk(b)+gc(6,b)))
      if (allocated(eps)) deallocate(eps) ; allocate(eps(1-gc(1,b):Ni(b)+gc(2,b),1-gc(3,b):Nj(b)+gc(4,b),1-gc(5,b):Nk(b)+gc(6,b)))
      read(unit_sol)(((ken(i,j,k),i=0,Ni(b)+1),j=0,Nj(b)+1),k=0,Nk(b)+1)
      read(unit_sol)(((eps(i,j,k),i=0,Ni(b)+1),j=0,Nj(b)+1),k=0,Nk(b)+1)
    endif
  endif
  if (forces) then ! allocating and initializing forces variables
    if (allocated(f_p)) deallocate(f_p)
    if (allocated(f_v)) deallocate(f_v)
    if (allocated(m_p)) deallocate(m_p)
    if (allocated(m_v)) deallocate(m_v)
    allocate(f_p(0-gc(1,b):Ni(b)+gc(2,b),0-gc(3,b):Nj(b)+gc(4,b),0-gc(5,b):Nk(b)+gc(6,b)))
    allocate(f_v(0-gc(1,b):Ni(b)+gc(2,b),0-gc(3,b):Nj(b)+gc(4,b),0-gc(5,b):Nk(b)+gc(6,b)))
    allocate(m_p(0-gc(1,b):Ni(b)+gc(2,b),0-gc(3,b):Nj(b)+gc(4,b),0-gc(5,b):Nk(b)+gc(6,b)))
    allocate(m_v(0-gc(1,b):Ni(b)+gc(2,b),0-gc(3,b):Nj(b)+gc(4,b),0-gc(5,b):Nk(b)+gc(6,b)))
    f_p = 0._R_P
    f_v = 0._R_P
    m_p = 0._R_P
    m_v = 0._R_P
  endif
  if (yp) then ! allocating y+ variable
    if (allocated(yplus))  deallocate(yplus)
    allocate(yplus(0-gc(1,b):Ni(b)+gc(2,b),0-gc(3,b):Nj(b)+gc(4,b),0-gc(5,b):Nk(b)+gc(6,b)))
  endif
  ! checking the presence of the patch
  patch_found = (any(ticc(0,:,:)==patch)).OR.(any(ticc(Ni(b)+1,:      ,:      )==patch)).OR. &
                (any(ticc(:,0,:)==patch)).OR.(any(ticc(:      ,Nj(b)+1,:      )==patch)).OR. &
                (any(ticc(:,:,0)==patch)).OR.(any(ticc(:      ,:      ,Nk(b)+1)==patch))
  patch_found = (patch_found.OR.((any(ticc(0,:,:)==patch+10)).OR.(any(ticc(Ni(b)+1,:      ,:      )==patch+10)).OR. &
                                 (any(ticc(:,0,:)==patch+10)).OR.(any(ticc(:      ,Nj(b)+1,:      )==patch+10)).OR. &
                                 (any(ticc(:,:,0)==patch+10)).OR.(any(ticc(:      ,:      ,Nk(b)+1)==patch+10))))
  if (patch_found) then
    istrml=0!3
    if (any(ticc(0,:,:)==patch).OR.any(ticc(0,:,:)==patch+10)) then
      Np = Np + 1

      ci1 = 1 ; ci2 = 1
      cj1 = 1 ; cj2 = Nj(b)
      ck1 = 1 ; ck2 = Nk(b)

      ni1 = 0+istrml ; ni2 = 0+istrml
      nj1 = cj1-1    ; nj2 = cj2
      nk1 = ck1-1    ; nk2 = ck2

      err = patch_save(ni1 = ni1, ni2 = ni2, &
                       nj1 = nj1, nj2 = nj2, &
                       nk1 = nk1, nk2 = nk2, &
                       ci1 = ci1, ci2 = ci2, &
                       cj1 = cj1, cj2 = cj2, &
                       ck1 = ck1, ck2 = ck2, &
                       face = 1, b = b)
    endif
    if (any(ticc(Ni(b)+1,:,:)==patch).OR.any(ticc(Ni(b)+1,:,:)==patch+10)) then
      Np = Np + 1

      ci1 = Ni(b) ; ci2 = Ni(b)
      cj1 = 1     ; cj2 = Nj(b)
      ck1 = 1     ; ck2 = Nk(b)

      ni1 = Ni(b)-istrml ; ni2 = Ni(b)-istrml
      nj1 = cj1-1 ; nj2 = cj2
      nk1 = ck1-1 ; nk2 = ck2

      err = patch_save(ni1 = ni1, ni2 = ni2, &
                       nj1 = nj1, nj2 = nj2, &
                       nk1 = nk1, nk2 = nk2, &
                       ci1 = ci1, ci2 = ci2, &
                       cj1 = cj1, cj2 = cj2, &
                       ck1 = ck1, ck2 = ck2, &
                       face = 2, b = b)
    endif
    if (any(ticc(:,0,:)==patch).OR.any(ticc(:,0,:)==patch+10)) then
      Np = Np + 1

      cj1 = 1 ; cj2 = 1
      ci1 = 1 ; ci2 = Ni(b)
      ck1 = 1 ; ck2 = Nk(b)

      nj1 = 0+istrml     ; nj2 = 0+istrml
      ni1 = ci1-1 ; ni2 = ci2
      nk1 = ck1-1 ; nk2 = ck2

      err = patch_save(ni1 = ni1, ni2 = ni2, &
                       nj1 = nj1, nj2 = nj2, &
                       nk1 = nk1, nk2 = nk2, &
                       ci1 = ci1, ci2 = ci2, &
                       cj1 = cj1, cj2 = cj2, &
                       ck1 = ck1, ck2 = ck2, &
                       face = 3, b = b)
    endif
    if (any(ticc(:,Nj(b)+1,:)==patch).OR.any(ticc(:,Nj(b)+1,:)==patch+10)) then
      Np = Np + 1

      cj1 = Nj(b) ; cj2 = Nj(b)
      ci1 = 1     ; ci2 = Ni(b)
      ck1 = 1     ; ck2 = Nk(b)

      nj1 = Nj(b)-istrml ; nj2 = Nj(b)-istrml
      ni1 = ci1-1 ; ni2 = ci2
      nk1 = ck1-1 ; nk2 = ck2

      err = patch_save(ni1 = ni1, ni2 = ni2, &
                       nj1 = nj1, nj2 = nj2, &
                       nk1 = nk1, nk2 = nk2, &
                       ci1 = ci1, ci2 = ci2, &
                       cj1 = cj1, cj2 = cj2, &
                       ck1 = ck1, ck2 = ck2, &
                       face = 4, b = b)
    endif
    if (any(ticc(:,:,0)==patch).OR.any(ticc(:,:,0)==patch+10)) then
      Np = Np + 1

      ck1 = 1 ; ck2 = 1
      ci1 = 1 ; ci2 = Ni(b)
      cj1 = 1 ; cj2 = Nj(b)

      nk1 = 0+istrml ; nk2 = 0+istrml
      ni1 = ci1-1 ; ni2 = ci2
      nj1 = cj1-1 ; nj2 = cj2

      err = patch_save(ni1 = ni1, ni2 = ni2, &
                       nj1 = nj1, nj2 = nj2, &
                       nk1 = nk1, nk2 = nk2, &
                       ci1 = ci1, ci2 = ci2, &
                       cj1 = cj1, cj2 = cj2, &
                       ck1 = ck1, ck2 = ck2, &
                       face = 5, b = b)
    endif
    if (any(ticc(:,:,Nk(b)+1)==patch).OR.any(ticc(:,:,Nk(b)+1)==patch+10)) then
      Np = Np + 1

      ck1 = Nk(b) ; ck2 = Nk(b)
      ci1 = 1     ; ci2 = Ni(b)
      cj1 = 1     ; cj2 = Nj(b)

      nk1 = Nk(b)-istrml ; nk2 = Nk(b)-istrml
      ni1 = ci1-1 ; ni2 = ci2
      nj1 = cj1-1 ; nj2 = cj2

      err = patch_save(ni1 = ni1, ni2 = ni2, &
                       nj1 = nj1, nj2 = nj2, &
                       nk1 = nk1, nk2 = nk2, &
                       ci1 = ci1, ci2 = ci2, &
                       cj1 = cj1, cj2 = cj2, &
                       ck1 = ck1, ck2 = ck2, &
                       face = 6, b = b)
    endif
  else
    cycle
  endif
enddo

close(unit_grd)
close(unit_icc)
if (allocated(rcc)) deallocate(rcc)
if (sol) then
  close(unit_sol)
endif
if (tec) then
  if (binary) then
    err = tecend112()
  else
    close(unit_out)
  endif
endif
if (vtk) then
  ! saving the multiblock VTM wrapper
  err = VTM_INI_XML(filename=adjustl(trim(File_out))//"."//trim(strz(2,patch))//'.vtm')
  err = VTM_BLK_XML(block_action='open')
  err = VTM_WRF_XML(flist=(/(adjustl(trim(basename(File_out)))//'-vts_files'//OS%sep//     &
                             adjustl(trim(basename(File_out)))//"."//trim(strz(2,patch))// &
                            '_n_'//trim(strz(5,b))//'.vts',                                &
                             b=1,Np)/))
  err = VTM_BLK_XML(block_action='close')
  err = VTM_END_XML()
endif
if (forces) then
  ! updating global file of forces
  write(unit_for,'(A)',iostat=err)trim(str(n=(fsum_v%x+fsum_p%x)))//' '// &
                                  trim(str(n=(fsum_v%y+fsum_p%y)))//' '// &
                                  trim(str(n=(fsum_v%z+fsum_p%z)))//' '// &
                                  trim(str(n=(         fsum_p%x)))//' '// &
                                  trim(str(n=(         fsum_p%y)))//' '// &
                                  trim(str(n=(         fsum_p%z)))//' '// &
                                  trim(str(n=(fsum_v%x         )))//' '// &
                                  trim(str(n=(fsum_v%y         )))//' '// &
                                  trim(str(n=(fsum_v%z         )))//' '// &
                                  trim(str(n=(msum_v%x+msum_p%x)))//' '// &
                                  trim(str(n=(msum_v%y+msum_p%y)))//' '// &
                                  trim(str(n=(msum_v%z+msum_p%z)))//' '// &
                                  trim(str(n=(         msum_p%x)))//' '// &
                                  trim(str(n=(         msum_p%y)))//' '// &
                                  trim(str(n=(         msum_p%z)))//' '// &
                                  trim(str(n=(msum_v%x         )))//' '// &
                                  trim(str(n=(msum_v%y         )))//' '// &
                                  trim(str(n=(msum_v%z         )))//' '// &
                                  trim(str(n=(Ssum             )))//' '// &
                                  'Fx,Fy,Fz,Fx_p,Fy_p,Fz_p,Fx_v,Fy_v,Fz_v,Mx,My,Mz,Mx_p,My_p,Mz_p,Mx_v,My_v,Mz_v,S of file "'&
                                  //adjustl(trim(File_sol))//'"'
  close(unit_for)
  if (Np>0) then
    if (allocated(Ni)) deallocate(Ni) ; allocate(Ni(1:1))
    if (allocated(Nj)) deallocate(Nj) ; allocate(Nj(1:1))
    if (allocated(Nk)) deallocate(Nk) ; allocate(Nk(1:1))
    write(unit_for_RB)Np
    rewind(unit_for_RB_scr)
    do
      read(unit_for_RB_scr,err=10,end=10)b,v,Ni(1),Nj(1),Nk(1)
      write(unit_for_RB)b,v,Ni(1),Nj(1),Nk(1)
      if (allocated(dummy)) deallocate(dummy) ; allocate(dummy(1:Ni(1),1:Nj(1),1:Nk(1)))
      do i=1,12
        read(unit_for_RB_scr,err=10)dummy
      enddo
    enddo
    10 continue
    rewind(unit_for_RB_scr)
    do
      read(unit_for_RB_scr,err=20,end=20)b,v,Ni(1),Nj(1),Nk(1)
      if (allocated(dummy)) deallocate(dummy) ; allocate(dummy(1:Ni(1),1:Nj(1),1:Nk(1)))
      do i=1,12
        read(unit_for_RB_scr,err=20,end=20)dummy
        write(unit_for_RB)dummy
      enddo
    enddo
    20 continue
  endif
  close(unit_for_RB)     ! Riccardo Broglia forces.RB
  close(unit_for_RB_scr) ! Riccardo Broglia forces.RB scratch file
  if (Ng>0) then
    do i=0,Ng
      close(unit_g_for(i)) ! forces-grp(#) file
    enddo
  endif
endif

if (Np==0) then
  ! no patches have been found
  write(stdout,'(A)')' No patches have been found!'
  write(stdout,'(A)')' Removing output file containing no patch data'
  if (tec) then
    if (binary) then
      err = remove_file(adjustl(trim(File_out))//"."//trim(strz(2,patch))//".plt")
    else
      err = remove_file(adjustl(trim(File_out))//"."//trim(strz(2,patch))//".dat")
    endif
  endif
  if (vtk) err = remove_file(adjustl(trim(File_out))//"."//trim(strz(2,patch))//'.vtm')
  if (forces) then
    err = remove_file(adjustl(trim(File_out))//"-forces.dat")
    err = remove_file(adjustl(trim(File_out))//"-forces.RB")  ! Riccardo Broglia forces.RB
  endif
endif
stop
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  subroutine print_usage()
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Printing XnPatches usage.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(stdout,*)
  write(stdout,'(A)') ' XnPatches'
  write(stdout,'(A)') ' Post processing code for Xnavis code'
  write(stdout,'(A)') ' Usage:'
  write(stdout,'(A)') '   XnPatches -g file_geo -i file_icc'
  write(stdout,'(A)') '            [-o file_output'
  write(stdout,'(A)') '             -s file_solution'
  write(stdout,'(A)') '             -p #bc_of_processed_patch'
  write(stdout,'(A)') '             -cell'
  write(stdout,'(A)') '             -ls'
  write(stdout,'(A)') '             -nt => no turbulent model, laminar flow'
  write(stdout,'(A)') '             -eq #turbulent_eq_model(0,1,2)'
  write(stdout,'(A)') '             -Re #Re'
  write(stdout,'(A)') '             -Fr #Fr'
  write(stdout,'(A)') '             -zfs #zfs'
  write(stdout,'(A)') '             -forces'
  write(stdout,'(A)') '             -yplus'
  write(stdout,'(A)') '             -metrics'
  write(stdout,'(A)') '             -ascii'
  write(stdout,'(A)') '             -tec yes/no'
  write(stdout,'(A)') '             -vtk yes/no'
  write(stdout,'(A)') '             -proc #proc'
  write(stdout,'(A)') '             -os UIX/WIN]'
  write(stdout,*)
  write(stdout,'(A)') ' Optional arguments meaning and default values:'
  write(stdout,'(A)') '  -o file_output   => output file name; default is basename of grd file with the proper extension'
  write(stdout,'(A)') '  -s file_solution => solution file name; if passed the solution variables are saved'
  write(stdout,'(A)') '  -p 1             => walls are default processed patches'
  write(stdout,'(A)') '  -cell            => all variables other than nodes coord. are cell centered (default no, node centered)'
  write(stdout,'(A)') '  -ls              => solution with level set, Froude number and free surface quote are necessary if'
  write(stdout,'(A)') '                      forces are computed (default no)'
  write(stdout,'(A)') '  -eq 1'
  write(stdout,'(A)') '  -Re #Re          => Reynolds number'
  write(stdout,'(A)') '  -Fr #Fr          => Froude number'
  write(stdout,'(A)') '  -zfs #zfs        => Quote of free surface'
  write(stdout,'(A)') '  -forces          => compute forces, Reynolds number is necessary (default no, do not compute them)'
  write(stdout,'(A)') '  -yplus           => compute y+, Reynolds number ans solution are necessary (default no, do not compute)'
  write(stdout,'(A)') '  -metrics         => save metrics variables (default no, do not save metrics)'
  write(stdout,'(A)') '  -ascii           => write ascii output file (default no, write binary one)'
  write(stdout,'(A)') '  -tec yes/no      => write (or not) Tecplot file format (default yes)'
  write(stdout,'(A)') '  -vtk yes/no      => write (or not) VTK file format (default no)'
  write(stdout,'(A)') '  -proc #proc      => if file "proc.input" if found global blocks numeration is used; #proc is the process'
  write(stdout,'(A)') '                      number of the current processed file'
  write(stdout,'(A)') '  -os UIX/WIN      => type of Operating System write (default *UIX OS type)'
  write(stdout,*)
  write(stdout,'(A)') ' Examples:'
  write(stdout,'(A)') '  XnPatches -g cc.01.grd -i cc.01                              -o mesh.01 (process only mesh of wall)'
  write(stdout,'(A)') '  XnPatches -g cc.01.grd -i cc.01 -s sol.00.01 -forces -Re 1d6 -o sol.01  (solution and forces are saved)'
  write(stdout,*)
  write(stdout,'(A)') ' Note:'
  write(stdout,'(A)') '   1) the output file name extension is not necessary because it assigned according to the type of output:'
  write(stdout,'(A)') '      binary       Tecplot => .plt'
  write(stdout,'(A)') '      ascii        Tecplot => .dat'
  write(stdout,'(A)') '      binary/ascii VTK     => .vtm'
  write(stdout,'(A)') '   2) if a file name "mb.par" is present into the directory where XnPatches is executed the values of Re,'
  write(stdout,'(A)') '      Fr, zfs and the turbulence model are loaded from this file thus they can be omitted from the command'
  write(stdout,'(A)') '      line arguments list.'
  write(stdout,'(A)') '   3) if forces are computed output file with suffix "-forces.dat" is also saved: this file contains the'
  write(stdout,'(A)') '      integral of the forces over all the patches and the "wet" surface.'
  write(stdout,'(A)') '   4) all the variables other than the nodes coordinates are saved at cell center if "-cell" option is used;'
  write(stdout,'(A)') '      if blanking is used the blanking mode must be "any corners" or "primary values".'
  write(stdout,*)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine print_usage

  subroutine parse_command_arguments()
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Subroutine for parsing command line arguments.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P)::   Nca = 0         ! Number of command line arguments.
  character(8)::   switch          ! Switch string.
  character(100):: string          ! Dummy string.
  character(3)::   os_type = 'UIX' ! OS type 'UIX' or 'WIN'.
  integer(I_P)::   c,e             ! Counter.
  integer(I_P)::   unitfree        ! Free logical unit.
  ! turbulence models inquiring flags
  logical:: balom = .false.
  logical:: sgs   = .false.
  logical:: spall = .false.
  logical:: des   = .false.
  logical:: ddes  = .false.
  logical:: lambr = .false.
  logical:: chang = .false.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  Nca = command_argument_count()
  if (Nca==0) then
    call print_usage
    stop
  endif
  c = 0
  do while (c<Nca)
    c = c + 1
    call get_command_argument(c,switch)
    select case(adjustl(trim(switch)))
    case('-g') ! geo file
      call get_command_argument(c+1,File_grd) ; c = c + 1
    case('-i') ! icc file
      call get_command_argument(c+1,File_icc) ; c = c + 1
    case('-s') ! solution file
      call get_command_argument(c+1,File_sol) ; c = c + 1
      sol = .true.
    case('-o') ! output file
      call get_command_argument(c+1,File_out) ; c = c + 1
    case('-p') ! type of patch to be post-processed
      call get_command_argument(c+1,string) ; c = c + 1
      read(string,*)patch
    case('-cell') ! dependent variables are saved cell centered instead of interpolated at nodes
      cell = .true.
    case('-ls') ! solution with level set
      level_set = .true.
    case('-Fr') ! Froude number
      call get_command_argument(c+1,string) ; c = c + 1
      read(string,*)Fr
      rFr2 = 1._R_P/(Fr*Fr)
    case('-zfs') ! z of free surface
      call get_command_argument(c+1,string) ; c = c + 1
      read(string,*)zfs
    case('-nt') ! no turbulent model
      laminar = .true.
       zeroeq = .false.
        oneeq = .false.
        twoeq = .false.
    case('-eq') ! turbulent equations model
      call get_command_argument(c+1,string)
      c = c + 1
      read(string,*)e
      select case(e)
      case(0)
      laminar = .false.
        zeroeq = .true.
        oneeq = .false.
        twoeq = .false.
      case(1)
      laminar = .false.
        zeroeq = .false.
        oneeq = .true.
        twoeq = .false.
      case(2)
      laminar = .false.
        zeroeq = .false.
        oneeq = .false.
        twoeq = .true.
      endselect
    case('-forces') ! forces computing
      forces = .true.
    case('-yplus') ! yplus computing
      yp = .true.
    case('-Re') ! Reynolds number
      call get_command_argument(c+1,string) ; c = c + 1
      read(string,*)Re
    case('-metrics') ! metrics info
      metrics = .true.
    case('-ascii') ! ascii output
        binary = .false.
    case('-tec') ! Tecplot file format
      call get_command_argument(c+1,string) ; c = c + 1
      string = Upper_Case(adjustl(trim(string)))
      tec = (adjustl(trim(string))=='YES')
    case('-vtk') ! VKT file format
      call get_command_argument(c+1,string) ; c = c + 1
      string = Upper_Case(adjustl(trim(string)))
      vtk = (adjustl(trim(string))=='YES')
    case('-proc') ! global blocks numeration actived
      call get_command_argument(c+1,string) ; c = c + 1
      read(string,*)myrank
      global_blk_num = .true.
    case('-os') ! OS type
      call get_command_argument(c+1,os_type) ; c = c + 1
      os_type = Upper_Case(os_type)
    case default
      write(stderr,'(A)') ' Unknown switch '
      write(stderr,'(A)') adjustl(trim(switch))
      write(stderr,*)
      call print_usage
      stop
    endselect
  enddo
  if ( laminar .and. (zeroeq .or. oneeq .or. twoeq) ) then
     write(stderr,'(A)') ' Incompatible switches, laminar disables turbulent model'
     stop
  end if
  ! set OS type
  call OS%init(c_id=os_type)
  if (adjustl(trim(File_grd))=='unset'.OR.adjustl(trim(File_icc))=='unset') then
    write(stderr,'(A)') ' File name of grd and icc files must be provided!'
    write(stderr,*)
    call print_usage
    stop
  endif
  if (adjustl(trim(File_out))=='unset') then
    ! the name of grd file is used as output file base name
    File_out = adjustl(trim(File_grd))
  endif
  ! converting the directory separators of files names according to the OS using
  File_grd = adjustl(trim(File_grd)) ; File_grd = string_OS_sep(File_grd) ! GRD file name.
  File_icc = adjustl(trim(File_icc)) ; File_icc = string_OS_sep(File_icc) ! ICC file name.
  File_sol = adjustl(trim(File_sol)) ; File_sol = string_OS_sep(File_sol) ! Solution file name.
  File_out = adjustl(trim(File_out)) ; File_out = string_OS_sep(File_out) ! Output file name.
  if (vtk) then
    vtkbdir  = adjustl(trim(basedir(File_out)))//adjustl(trim(basename(File_out)))//'-vts_files'//OS%sep
    vtkbfile = adjustl(trim(vtkbdir))//adjustl(trim(basename(File_out)))
  endif
           write(stdout,'(A)')' GRD File '//trim(File_grd)
           write(stdout,'(A)')' ICC File '//trim(File_icc)
  if (sol) write(stdout,'(A)')' SOL File '//trim(File_sol)
           write(stdout,'(A)')' OUT File '//trim(File_out)
  if (sol) then
    Nvar=8
    if (binary) then
      tecvarname = 'x y z icc u v w p'
    else
      tecvarname = ' VARIABLES ="x" "y" "z" "icc" "u" "v" "w" "p"'
    endif
    if (zeroeq) then
      Nvar = Nvar + 1 ! visc must be saved
      if (binary) then
        tecvarname = trim(tecvarname)//' visc'
      else
        tecvarname = trim(tecvarname)//' "visc"'
      endif
    elseif (oneeq) then
      Nvar = Nvar + 2 ! visc and vitl must be saved
      if (binary) then
        tecvarname = trim(tecvarname)//' visc vitl'
      else
        tecvarname = trim(tecvarname)//' "visc" "vitl"'
      endif
    elseif (twoeq) then
      Nvar = Nvar + 3 ! visc, ken and eps must be saved
      if (binary) then
        tecvarname = trim(tecvarname)//' visc ken eps'
      else
        tecvarname = trim(tecvarname)//' "visc" "ken" "eps"'
      endif
    endif
    if (level_set) then
      Nvar = Nvar + 2 ! f and f0 must be saved
      if (binary) then
        tecvarname = trim(tecvarname)//' f f0'
      else
        tecvarname = trim(tecvarname)//' "f" "f0"'
      endif
    endif
    if (forces) then
      Nvar = Nvar + 18 ! fx, fy, fz, mx, my and mz (total, pressure and viscous parts)
      if (binary) then
        tecvarname = trim(tecvarname)//' fx fy fz fx_p fy_p fz_p fx_v fy_v fz_v mx my mz mx_p my_p mz_p mx_v my_v mz_v'
      else
        tecvarname = trim(tecvarname)//' "fx" "fy" "fz" "fx_p" "fy_p" "fz_p" "fx_v" "fy_v" "fz_v"'//&
                                       ' "mx" "my" "mz" "mx_p" "my_p" "mz_p" "mx_v" "my_v" "mz_v"'
      endif
    endif
    if (yp) then
      Nvar = Nvar + 1 ! yplus
      if (binary) then
        tecvarname = trim(tecvarname)//' yplus'
      else
        tecvarname = trim(tecvarname)//' "yplus"'
      endif
    endif
  else
    Nvar=4
    if (binary) then
      tecvarname = 'x y z icc'
    else
      tecvarname = ' VARIABLES ="x" "y" "z" "icc"'
    endif
  endif
  if (metrics) then
    Nvar = Nvar + 4 ! Nx,Ny,Nz and S
    if (binary) then
      tecvarname = trim(tecvarname)//' Nx Ny Nz S'
    else
      tecvarname = trim(tecvarname)//' "Nx" "Ny" "Nz" "S"'
    endif
  endif
  ! inquiring the presence of "mb.par"
  inquire(file="mb.par",exist=is_file,iostat=err)
  if (is_file) then
    unitfree = Get_Unit() ; open(unit=unitfree,file="mb.par",action='READ')
    do e=1,8 ! skip the first 8 records
      read(unitfree,*)
    enddo
    read(unitfree,*) e      ! number grid levels
    do e=1,e+9 ! skip the e+9 records
      read(unitfree,*)
    enddo
    read(unitfree,*) Re     ! Reynolds number
    read(unitfree,*) Fr     ! Foude number
    read(unitfree,*)
    read(unitfree,*) string ! turbulence model
    read(unitfree,*)
    read(unitfree,*) zfs    ! quote of free surface
    close(unitfree)
    if (Fr>0._R_P) then
      level_set = .true.
      rFr2 = 1._R_P/(Fr*Fr)
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
    end if
    if (string(1:3)=='dde'.or.string(1:3)=='DDE') then
       spall = .true.
       des   = .true.
       ddes  = .true.
    end if
    zeroeq = (sgs.OR.balom)
    oneeq  = spall
    twoeq  = (lambr.OR.chang)
    write(stdout,'(A)') ' Found mb.par'
    write(stdout,'(A)') ' Reynolds number '//str(n=Re)
    write(stdout,'(A)') ' Froude   number '//str(n=Fr)
    if (level_set) then
      write(stdout,'(A)') ' Free sruface present'
    else
      write(stdout,'(A)') ' Free sruface not present'
    endif
    if (zeroeq) then
      write(stdout,'(A)') ' Zero equation turbulence model'
    elseif (oneeq) then
      write(stdout,'(A)') ' One equation turbulence model'
    elseif(twoeq) then
      write(stdout,'(A)') ' Two equations turbulence model'
    endif
  endif
  ! checking that all the info for forces computing have been provided
  if (forces) then
    if (.not.sol) then
      write(stderr,'(A)')' Trying to compute forces without providing solution file!'
      stop
    endif
    if (Re<0._R_P) then
      write(stderr,'(A)')' Trying to compute forces without providing Reynolds number!'
      write(stderr,'(A)')' The Reynolds number must be provide either by command line argument or by mb.par'
      stop
    endif
    if (level_set) then
      if (Fr<0._R_P) then
        write(stderr,'(A)')' Trying to compute forces without providing Froude number for a simulation with free surface!'
        write(stderr,'(A)')' The Froude number must be provide either by command line argument or by mb.par'
        stop
      elseif (zfs>MaxR_P*0.9_R_P) then
        write(stderr,'(A)')' Trying to compute forces without providing quote of free surface for a simulation with free surface!'
        write(stderr,'(A)')' The quote of free surface must be provide either by command line argument or by mb.par'
        stop
      endif
    endif
  endif
  ! checking that all the info for yplus computing have been provided
  if (yp) then
    if (.not.sol) then
      write(stderr,'(A)')' Trying to compute y+ without providing solution file!'
      stop
    endif
    if (Re<0._R_P) then
      write(stderr,'(A)')' Trying to compute y+ without providing Reynolds number!'
      write(stderr,'(A)')' The Reynolds number must be provide either by command line argument or by mb.par'
      stop
    endif
  endif
  ! inquiring the presence of "proc.input" for global blocks numeration
  if (global_blk_num) then
    inquire(file='proc.input',exist=is_file)
    if (.not.is_file) then
      call File_Not_Found(myrank,'proc.input','parse_command_arguments')
    else
      unitfree = Get_Unit() ; open(unit=unitfree,file="proc.input",action='READ')
      read(unitfree,*) ! record skipped
      read(unitfree,*) ! record skipped
      read(unitfree,*) ! record skipped
      read(unitfree,*) Nbtot ! number of total blocks
      if (allocated(procmap)) deallocate(procmap) ; allocate(procmap(1:2,1:Nbtot)) ; procmap = 0_I_P
      read(unitfree,*) ! record skipped
      read(unitfree,*) ! record skipped
      read(unitfree,*) ! record skipped
      do b=1,Nbtot
        read(unitfree,*) c,procmap(1,b),c,procmap(2,b)
      end do
      close(unitfree)
      nprocs = maxval(procmap(2,:))
      Ng     = maxval(procmap(1,:))
      ! computing the local (of myrank) number of blocks
      Nb = 0
      do b=1,Nbtot
        if (procmap(2,b)==myrank) Nb = Nb + 1
      enddo
      if (allocated(blockmap)) deallocate(blockmap) ; allocate(blockmap(1:2,1:Nb)) ; blockmap = 0_I_P
      c = 0
      do b=1,Nbtot
        if (procmap(2,b)==myrank) then
          c = c + 1
          blockmap(1,c) = procmap(1,b)
          blockmap(2,c) = b
        endif
      enddo
    endif
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine parse_command_arguments

  subroutine compute_metrics(b)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Subroutine for computing block metrics.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P), intent(IN):: b                 ! Actual block number.
  type(Type_Vector)::        NFS,s1,s2,nd,db   ! Dummy vector variables.
  real(R_P)::                signi,signj,signk ! Dummy variables for checking the directions of normals.
  real(R_P)::                Vx,Vy,Vz          ! Dummy variables for computing volume.
  real(R_P)::                xp,yp,zp          ! Dummy variables for computing face coordinates.
  real(R_P)::                xm,ym,zm          ! Dummy variables for computing face coordinates.
  integer(I_P)::             i,j,k             ! Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! computing faces normals
  ! positioning at the middle of the block
  i = max(1,Ni(b)/2)
  j = max(1,Nj(b)/2)
  k = max(1,Nk(b)/2)
  ! checking the direction of i normals
  s1 = node(i,j  ,k) - node(i,j-1,k-1)
  s2 = node(i,j-1,k) - node(i,j,  k-1)
  nd = s1.cross.s2
  s1 = 0.25_R_P*(node(i,  j,k)+node(i,  j-1,k)+node(i,  j,k-1)+node(i,  j-1,k-1))
  s2 = 0.25_R_P*(node(i-1,j,k)+node(i-1,j-1,k)+node(i-1,j,k-1)+node(i-1,j-1,k-1))
  db = s1 - s2
  signi = sign(1._R_P,(nd.dot.db))
  ! checking the direction of j normals
  s1 = node(i,j,k  ) - node(i-1,j,k-1)
  s2 = node(i,j,k-1) - node(i-1,j,k  )
  nd = s1.cross.s2
  s1 = 0.25_R_P*(node(i,j,  k)+node(i-1,j,  k)+node(i,j,  k-1)+node(i-1,j,  k-1))
  s2 = 0.25_R_P*(node(i,j-1,k)+node(i-1,j-1,k)+node(i,j-1,k-1)+node(i-1,j-1,k-1))
  db = s1 - s2
  signj = sign(1._R_P,(nd.dot.db))
  ! checking the direction of k normals
  s1 = node(i,  j,k) - node(i-1,j-1,k)
  s2 = node(i-1,j,k) - node(i,  j-1,k)
  nd = s1.cross.s2
  s1 = 0.25_R_P*(node(i,j,k  )+node(i-1,j,k  )+node(i,j-1,k  )+node(i-1,j-1,k  ))
  s2 = 0.25_R_P*(node(i,j,k-1)+node(i-1,j,k-1)+node(i,j-1,k-1)+node(i-1,j-1,k-1))
  db = s1 - s2
  signk = sign(1._R_P,(nd.dot.db))
  !$OMP PARALLEL DEFAULT(NONE)                        &
  !$OMP PRIVATE(i,j,k,NFS,Vx,Vy,Vz,xp,yp,zp,xm,ym,zm) &
  !$OMP SHARED(b,Ni,Nj,Nk,gc,signi,signj,signk,node,Si,Sj,Sk,NFi,NFj,NFk,NFiS,NFjS,NFkS,volume)
  !$OMP DO
  do k=1-gc(5,b),Nk(b)+gc(6,b)
    do j=1-gc(3,b),Nj(b)+gc(4,b)
      do i=0-gc(1,b),Ni(b)+gc(2,b)
        NFS = face_normal4(pt1 = node(i,j-1,k-1), &
                           pt2 = node(i,j  ,k-1), &
                           pt3 = node(i,j  ,k  ), &
                           pt4 = node(i,j-1,k  ))
        NFS = NFS*signi
        NFiS(i,j,k) = NFS
        NFi (i,j,k) = normalize(NFS)
        Si  (i,j,k) =    normL2(NFS)
      enddo
    enddo
  enddo
  !$OMP DO
  do k=1-gc(5,b),Nk(b)+gc(6,b)
    do j=0-gc(3,b),Nj(b)+gc(4,b)
      do i=1-gc(1,b),Ni(b)+gc(2,b)
        NFS = face_normal4(pt1 = node(i-1,j,k-1), &
                           pt2 = node(i-1,j,k  ), &
                           pt3 = node(i  ,j,k  ), &
                           pt4 = node(i  ,j,k-1))
        NFS = NFS*signj
        NFjS(i,j,k) = NFS
        NFj (i,j,k) = normalize(NFS)
        Sj  (i,j,k) =    normL2(NFS)
      enddo
    enddo
  enddo
  !$OMP DO
  do k=0-gc(5,b),Nk(b)+gc(6,b)
    do j=1-gc(3,b),Nj(b)+gc(4,b)
      do i=1-gc(1,b),Ni(b)+gc(2,b)
        NFS = face_normal4(pt1 = node(i-1,j-1,k), &
                           pt2 = node(i  ,j-1,k), &
                           pt3 = node(i  ,j  ,k), &
                           pt4 = node(i-1,j  ,k))
        NFS = NFS*signk
        NFkS(i,j,k) = NFS
        NFk (i,j,k) = normalize(NFS)
        Sk  (i,j,k) =    normL2(NFS)
      enddo
    enddo
  enddo
  ! computing cells volumes
  !$OMP DO
  do k=1-gc(5,b),Nk(b)+gc(6,b)
    do j=1-gc(3,b),Nj(b)+gc(4,b)
      do i=1-gc(1,b),Ni(b)+gc(2,b)
        Vx = 0._R_P
        Vy = 0._R_P
        Vz = 0._R_P

        xp = 0.25_R_P*(node(i  ,j  ,k  )%x + node(i  ,j  ,k-1)%x + &
                       node(i  ,j-1,k  )%x + node(i  ,j-1,k-1)%x)
        yp = 0.25_R_P*(node(i  ,j  ,k  )%y + node(i  ,j  ,k-1)%y + &
                       node(i  ,j-1,k  )%y + node(i  ,j-1,k-1)%y)
        zp = 0.25_R_P*(node(i  ,j  ,k  )%z + node(i  ,j  ,k-1)%z + &
                       node(i  ,j-1,k  )%z + node(i  ,j-1,k-1)%z)
        xm = 0.25_R_P*(node(i-1,j  ,k  )%x + node(i-1,j  ,k-1)%x + &
                       node(i-1,j-1,k  )%x + node(i-1,j-1,k-1)%x)
        ym = 0.25_R_P*(node(i-1,j  ,k  )%y + node(i-1,j  ,k-1)%y + &
                       node(i-1,j-1,k  )%y + node(i-1,j-1,k-1)%y)
        zm = 0.25_R_P*(node(i-1,j  ,k  )%z + node(i-1,j  ,k-1)%z + &
                       node(i-1,j-1,k  )%z + node(i-1,j-1,k-1)%z)

        Vx = Vx + xp*NFiS(i,j,k)%x - xm*NFiS(i-1,j,k)%x
        Vy = Vy + yp*NFiS(i,j,k)%y - ym*NFiS(i-1,j,k)%y
        Vz = Vz + zp*NFiS(i,j,k)%z - zm*NFiS(i-1,j,k)%z

        xp = 0.25_R_P*(node(i  ,j  ,k  )%x + node(i  ,j  ,k-1)%x + &
                       node(i-1,j  ,k  )%x + node(i-1,j  ,k-1)%x)
        yp = 0.25_R_P*(node(i  ,j  ,k  )%y + node(i  ,j  ,k-1)%y + &
                       node(i-1,j  ,k  )%y + node(i-1,j  ,k-1)%y)
        zp = 0.25_R_P*(node(i  ,j  ,k  )%z + node(i  ,j  ,k-1)%z + &
                       node(i-1,j  ,k  )%z + node(i-1,j  ,k-1)%z)
        xm = 0.25_R_P*(node(i  ,j-1,k  )%x + node(i  ,j-1,k-1)%x + &
                       node(i-1,j-1,k  )%x + node(i-1,j-1,k-1)%x)
        ym = 0.25_R_P*(node(i  ,j-1,k  )%y + node(i  ,j-1,k-1)%y + &
                       node(i-1,j-1,k  )%y + node(i-1,j-1,k-1)%y)
        zm = 0.25_R_P*(node(i  ,j-1,k  )%z + node(i  ,j-1,k-1)%z + &
                       node(i-1,j-1,k  )%z + node(i-1,j-1,k-1)%z)

        Vx = Vx + xp*NFjS(i,j,k)%x - xm*NFjS(i,j-1,k)%x
        Vy = Vy + yp*NFjS(i,j,k)%y - ym*NFjS(i,j-1,k)%y
        Vz = Vz + zp*NFjS(i,j,k)%z - zm*NFjS(i,j-1,k)%z

        xp = 0.25_R_P*(node(i  ,j  ,k  )%x + node(i  ,j-1,k  )%x + &
                       node(i-1,j  ,k  )%x + node(i-1,j-1,k  )%x)
        yp = 0.25_R_P*(node(i  ,j  ,k  )%y + node(i  ,j-1,k  )%y + &
                       node(i-1,j  ,k  )%y + node(i-1,j-1,k  )%y)
        zp = 0.25_R_P*(node(i  ,j  ,k  )%z + node(i  ,j-1,k  )%z + &
                       node(i-1,j  ,k  )%z + node(i-1,j-1,k  )%z)
        xm = 0.25_R_P*(node(i  ,j  ,k-1)%x + node(i  ,j-1,k-1)%x + &
                       node(i-1,j  ,k-1)%x + node(i-1,j-1,k-1)%x)
        ym = 0.25_R_P*(node(i  ,j  ,k-1)%y + node(i  ,j-1,k-1)%y + &
                       node(i-1,j  ,k-1)%y + node(i-1,j-1,k-1)%y)
        zm = 0.25_R_P*(node(i  ,j  ,k-1)%z + node(i  ,j-1,k-1)%z + &
                       node(i-1,j  ,k-1)%z + node(i-1,j-1,k-1)%z)

        Vx = Vx + xp*NFkS(i,j,k)%x - xm*NFkS(i,j,k-1)%x
        Vy = Vy + yp*NFkS(i,j,k)%y - ym*NFkS(i,j,k-1)%y
        Vz = Vz + zp*NFkS(i,j,k)%z - zm*NFkS(i,j,k-1)%z

        volume(i,j,k) = max(Vx,Vy,Vz)
      enddo
    enddo
  enddo
  !$OMP END PARALLEL
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine compute_metrics

  subroutine bc_metrics_correction(b)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Subroutine for correcting the metrics of natural (and negative volume) boundary conditions cells.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P), intent(IN):: b          ! Actual block number.
  logical::                  bc_correct ! Flag for inquiring if the bc metrics must be corrected.
  logical::                  bc_wall    ! Flag for inquiring if the bc is "wall-type": different corrections must be used.
  real(R_P)::                tm         ! Tangential metrics parameter (-1 for wall-type bc).
  real(R_P)::                sn         ! Normal     metrics coefficient correction.
  integer(I_P)::             i,j,k      ! counters.
  integer(I_P), parameter::  wall         = -1
  integer(I_P), parameter::  simmetry     = -2
  integer(I_P), parameter::  movingwall   = -10
  integer(I_P), parameter::  passivewall  = -11
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  !$OMP PARALLEL DEFAULT(NONE)                  &
  !$OMP PRIVATE(i,j,k,bc_correct,bc_wall,tm,sn) &
  !$OMP SHARED(b,Ni,Nj,Nk,rcc,icc,NFi,NFj,NFk,NFiS,NFjS,NFkS,volume)
  ! left i
  !$OMP DO
  do k=1,Nk(b)
    do j=1,Nj(b)
      i = nint(rcc(icc(0,j,k)))
      bc_correct = ((i<0).OR.(volume(0,j,k)<(0.2_R_P*volume(1,j,k))))
      bc_wall    = ((i==wall).OR.(i==simmetry).OR.(i==movingwall).OR.(i==passivewall))
      tm = 1._R_P
      if (bc_wall) tm = -1._R_P
      if (bc_correct) then
         ! normal metrics
         sn = 2._R_P*(NFiS(1,j,k).dot.NFi(0,j,k))
         NFiS( -1,j,  k  ) = -NFiS(1,j,k)+sn*NFi(0,j,k)
         ! tangential metrics
         NFjS(  0,j  ,k  ) = tm*NFjS(1,j  ,k  )
         NFjS(  0,j-1,k  ) = tm*NFjS(1,j-1,k  )
         NFjS(  0,j  ,k-1) = tm*NFjS(1,j  ,k-1)
         NFjS(  0,j-1,k-1) = tm*NFjS(1,j-1,k-1)

         NFkS(  0,j  ,k  ) = tm*NFkS(1,j  ,k  )
         NFkS(  0,j-1,k  ) = tm*NFkS(1,j-1,k  )
         NFkS(  0,j  ,k-1) = tm*NFkS(1,j  ,k-1)
         NFkS(  0,j-1,k-1) = tm*NFkS(1,j-1,k-1)
         ! volume
         volume(0,j,  k  ) = volume(1,j,k)
      end if
    enddo
  enddo
  ! right i
  !$OMP DO
  do k=1,Nk(b)
    do j=1,Nj(b)
      i = nint(rcc(icc(Ni(b)+1,j,k)))
      bc_correct = ((i<0).OR.(volume(Ni(b)+1,j,k)<(0.2_R_P*volume(Ni(b),j,k))))
      bc_wall    = ((i==wall).OR.(i==simmetry).OR.(i==movingwall).OR.(i==passivewall))
      tm = 1._R_P
      if (bc_wall) tm = -1._R_P
      if (bc_correct) then
         ! normal metrics
         sn = 2._R_P*(NFiS(Ni(b)-1,j,k).dot.NFi(Ni(b),j,k))
         NFiS(  Ni(b)+1,j,  k  ) = -NFiS(Ni(b)-1,j,k)+sn*NFi(Ni(b),j,k)
         ! tangential metrics
         NFjS(  Ni(b)+1,j  ,k  ) = tm*NFjS(Ni(b),j  ,k  )
         NFjS(  Ni(b)+1,j-1,k  ) = tm*NFjS(Ni(b),j-1,k  )
         NFjS(  Ni(b)+1,j  ,k-1) = tm*NFjS(Ni(b),j  ,k-1)
         NFjS(  Ni(b)+1,j-1,k-1) = tm*NFjS(Ni(b),j-1,k-1)

         NFkS(  Ni(b)+1,j  ,k  ) = tm*NFkS(Ni(b),j  ,k  )
         NFkS(  Ni(b)+1,j-1,k  ) = tm*NFkS(Ni(b),j-1,k  )
         NFkS(  Ni(b)+1,j  ,k-1) = tm*NFkS(Ni(b),j  ,k-1)
         NFkS(  Ni(b)+1,j-1,k-1) = tm*NFkS(Ni(b),j-1,k-1)
         ! volume
         volume(Ni(b)+1,j,  k  ) = volume(Ni(b),j,k)
      end if
    enddo
  enddo
  ! left j
  !$OMP DO
  do k=1,Nk(b)
    do i=1,Ni(b)
      j = nint(rcc(icc(i,0,k)))
      bc_correct = ((j<0).OR.(volume(i,0,k)<(0.2_R_P*volume(i,1,k))))
      bc_wall    = ((j==wall).OR.(j==simmetry).OR.(j==movingwall).OR.(j==passivewall))
      tm = 1._R_P
      if (bc_wall) tm = -1._R_P
      if (bc_correct) then
         ! normal metrics
         sn = 2._R_P*(NFjS(i,1,k).dot.NFj(i,0,k))
         NFjS(  i, -1,k  ) = -NFjS(i,1,k)+sn*NFj(i,0,k)
         ! tangential metrics
         NFiS(  i  ,0,k  ) = tm*NFiS(i  ,1,k  )
         NFiS(  i-1,0,k  ) = tm*NFiS(i-1,1,k  )
         NFiS(  i  ,0,k-1) = tm*NFiS(i  ,1,k-1)
         NFiS(  i-1,0,k-1) = tm*NFiS(i-1,1,k-1)

         NFkS(  i  ,0,k  ) = tm*NFkS(i  ,1,k  )
         NFkS(  i-1,0,k  ) = tm*NFkS(i-1,1,k  )
         NFkS(  i  ,0,k-1) = tm*NFkS(i  ,1,k-1)
         NFkS(  i-1,0,k-1) = tm*NFkS(i-1,1,k-1)
         ! volume
         volume(i,  0,k  ) = volume(i,1,k)
      end if
    enddo
  enddo
  ! right j
  !$OMP DO
  do k=1,Nk(b)
    do i=1,Ni(b)
      j = nint(rcc(icc(i,Nj(b)+1,k)))
      bc_correct = ((j<0).OR.(volume(i,Nj(b)+1,k)<(0.2_R_P*volume(i,Nj(b),k))))
      bc_wall    = ((j==wall).OR.(j==simmetry).OR.(j==movingwall).OR.(j==passivewall))
      tm = 1._R_P
      if (bc_wall) tm = -1._R_P
      if (bc_correct) then
         ! normal metrics
         sn = 2._R_P*(NFjS(i,Nj(b)-1,k).dot.NFj(i,Nj(b),k))
         NFjS(  i,Nj(b)+1,  k  ) = -NFjS(i,Nj(b)-1,k)+sn*NFj(i,Nj(b),k)
         ! tangential metrics
         NFiS(  i  ,Nj(b)+1,k  ) = tm*NFiS(i  ,Nj(b),k  )
         NFiS(  i-1,Nj(b)+1,k  ) = tm*NFiS(i-1,Nj(b),k  )
         NFiS(  i  ,Nj(b)+1,k-1) = tm*NFiS(i  ,Nj(b),k-1)
         NFiS(  i-1,Nj(b)+1,k-1) = tm*NFiS(i-1,Nj(b),k-1)

         NFkS(  i  ,Nj(b)+1,k  ) = tm*NFkS(i  ,Nj(b),k  )
         NFkS(  i-1,Nj(b)+1,k  ) = tm*NFkS(i-1,Nj(b),k  )
         NFkS(  i  ,Nj(b)+1,k-1) = tm*NFkS(i  ,Nj(b),k-1)
         NFkS(  i-1,Nj(b)+1,k-1) = tm*NFkS(i-1,Nj(b),k-1)
         ! volume
         volume(i,Nj(b)+1,  k  ) = volume(i,Nj(b),k)
      end if
    enddo
  enddo
  ! left k
  !$OMP DO
  do j=1,Nj(b)
    do i=1,Ni(b)
      k = nint(rcc(icc(i,j,0)))
      bc_correct = ((k<0).OR.(volume(i,j,0)<(0.2_R_P*volume(i,j,1))))
      bc_wall    = ((k==wall).OR.(k==simmetry).OR.(k==movingwall).OR.(k==passivewall))
      tm = 1._R_P
      if (bc_wall) tm = -1._R_P
      if (bc_correct) then
         ! normal metrics
         sn = 2._R_P*(NFkS(i,j,1).dot.NFk(i,j,0))
         NFkS(  i,  j, -1) = -NFkS(i,j,1)+sn*NFk(i,j,0)
         ! tangential metrics
         NFiS(  i  ,j  ,0) = tm*NFiS(i  ,j  ,1)
         NFiS(  i-1,j  ,0) = tm*NFiS(i-1,j  ,1)
         NFiS(  i  ,j-1,0) = tm*NFiS(i  ,j-1,1)
         NFiS(  i-1,j-1,0) = tm*NFiS(i-1,j-1,1)

         NFjS(  i  ,j  ,0) = tm*NFjS(i  ,j  ,1)
         NFjS(  i-1,j  ,0) = tm*NFjS(i-1,j  ,1)
         NFjS(  i  ,j-1,0) = tm*NFjS(i  ,j-1,1)
         NFjS(  i-1,j-1,0) = tm*NFjS(i-1,j-1,1)
         ! volume
         volume(i,  j,  0) = volume(i,j,1)
      end if
    enddo
  enddo
  ! right k
  !$OMP DO
  do j=1,Nj(b)
    do i=1,Ni(b)
      k = nint(rcc(icc(i,j,Nk(b)+1)))
      bc_correct = ((k<0).OR.(volume(i,j,Nk(b)+1)<(0.2_R_P*volume(i,j,Nk(b)))))
      bc_wall    = ((k==wall).OR.(k==simmetry).OR.(k==movingwall).OR.(k==passivewall))
      tm = 1._R_P
      if (bc_wall) tm = -1._R_P
      if (bc_correct) then
         ! normal metrics
         sn = 2._R_P*(NFkS(i,j,Nk(b)-1).dot.NFk(i,j,Nk(b)))
         NFkS(  i,  j,  Nk(b)+1) = -NFkS(i,j,Nk(b)-1)+sn*NFk(i,j,Nk(b))
         ! tangential metrics
         NFiS(  i  ,j  ,Nk(b)+1) = tm*NFiS(i  ,j  ,Nk(b))
         NFiS(  i-1,j  ,Nk(b)+1) = tm*NFiS(i-1,j  ,Nk(b))
         NFiS(  i  ,j-1,Nk(b)+1) = tm*NFiS(i  ,j-1,Nk(b))
         NFiS(  i-1,j-1,Nk(b)+1) = tm*NFiS(i-1,j-1,Nk(b))

         NFjS(  i  ,j  ,Nk(b)+1) = tm*NFjS(i  ,j  ,Nk(b))
         NFjS(  i-1,j  ,Nk(b)+1) = tm*NFjS(i-1,j  ,Nk(b))
         NFjS(  i  ,j-1,Nk(b)+1) = tm*NFjS(i  ,j-1,Nk(b))
         NFjS(  i-1,j-1,Nk(b)+1) = tm*NFjS(i-1,j-1,Nk(b))
         ! volume
         volume(i,  j,  Nk(b)+1) = volume(i,j,Nk(b))
      end if
    enddo
  enddo
  !$OMP END PARALLEL
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine bc_metrics_correction

  function patch_save(ni1,ni2,nj1,nj2,nk1,nk2,ci1,ci2,cj1,cj2,ck1,ck2,face,b) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Function for saving the patch.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P),      intent(IN):: ni1,ni2                 ! First and last node i indexes.
  integer(I_P),      intent(IN):: nj1,nj2                 ! First and last node j indexes.
  integer(I_P),      intent(IN):: nk1,nk2                 ! First and last node k indexes.
  integer(I_P),      intent(IN):: ci1,ci2                 ! First and last cell i indexes.
  integer(I_P),      intent(IN):: cj1,cj2                 ! First and last cell j indexes.
  integer(I_P),      intent(IN):: ck1,ck2                 ! First and last cell k indexes.
  integer(I_P),      intent(IN):: face                    ! Face where patch is defined: 1,2,3,4,5,6.
  integer(I_P),      intent(IN):: b                       ! Actual block number.
  integer(I_P)::                  err                     ! Error traping flag: 0 no errors, >0 error occours.
  real(R_P), allocatable::        vari(:,:,:)             ! Interpolated generic variables.
  integer(I_P)::                  tecnull(1:Nvar)         ! Tecplot null array.
  integer(I_P)::                  tecvarloc(1:Nvar)       ! Tecplot array of variables location.
  integer(I_P)::                  nnode,ncell             ! Number of nodes and cells.
  real(R_P)::                     fsumpx,fsumpy,fsumpz, & ! |
                                  fsumvx,fsumvy,fsumvz, & ! | Dummies vars for computing sum of forces in OpenMP parallel blocks.
                                  msumpx,msumpy,msumpz, & ! | OpenMP doesn't support reduction on derived type vars.
                                  msumvx,msumvy,msumvz    ! |
  type(Type_Vector)::             NdS                     ! Normal "tilde" for viscous part of forces (or distance for y+).
  type(Type_Vector)::             fco                     ! Pacth center coordinates.
  real(R_P)::                     hyd                     ! Hydrostatic correction of pressure term.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.cell) then ! variables must be interpolated at nodes, allocating dummy variable
    if (allocated(vari)) deallocate(vari) ; allocate(vari(ni1:ni2,nj1:nj2,nk1:nk2))
  endif
  if (forces) then ! computing the forces acting on the patch
    fsumpx = 0._R_P ; fsumpy = 0._R_P ; fsumpz = 0._R_P
    fsumvx = 0._R_P ; fsumvy = 0._R_P ; fsumvz = 0._R_P
    msumpx = 0._R_P ; msumpy = 0._R_P ; msumpz = 0._R_P
    msumvx = 0._R_P ; msumvy = 0._R_P ; msumvz = 0._R_P
    select case(face)
    case(1,2)
      !$OMP PARALLEL DEFAULT(NONE)                                                                        &
      !$OMP PRIVATE(i,j,k,NdS,hyd)                                                                        &
      !$OMP SHARED(ni1,ni2,nj1,nj2,nk1,nk2,ci1,ci2,cj1,cj2,ck1,ck2,face,patch,rFr2,zfs,                   &
      !$OMP        node,NFiS,NFjS,NFkS,volume,icc,ticc,f0,momentum,pressure,f_p,f_v,m_p,m_v,level_set,Re) &
      !$OMP REDUCTION(+: fsumpx,fsumpy,fsumpz,fsumvx,fsumvy,fsumvz,Ssum)
      !$OMP DO
      do k=ck1,ck2
        do j=cj1,cj2
          ! checking if this is an active cell
          if (ticc(ni1-1+face,j,k)/=patch) cycle
          if ( icc(ci1,j,k)/=0) cycle
          if(level_set) then
            if (f0(ni1-1+face,j,k)>0._R_P) cycle
          endif

          ! patch center coordinates
          fco = 0.25_R_P*(node(ni1,j,k)+node(ni1,j-1,k)+node(ni1,j,k-1)+node(ni1,j-1,k-1))

          ! pressure part of forces
          hyd = 0._R_P
          if(level_set) hyd = (fco%z - zfs)*rFr2
          f_p(ci1,j,k) = -(1.5_R_P*pressure(ci1,j,k) - 0.5_R_P*pressure(ci1+3-2*face,j,k) - hyd)*NFiS(ni1,j,k)
          ! viscous part of forces
          NdS = (2._R_P*NFiS(ni1,j,k) + NFiS(ni1+1,j,k) + NFiS(ni1-1,j,k))/(2._R_P*(volume(ni1,j,k)+volume(ni1+1,j,k)))
          f_v(ci1,j,k) = (momentum(ci1-1+face,j,k) - momentum(ci1-2+face,j,k))*(NFiS(ni1,j,k).dot.NdS)/Re

          if (face==2) then
            f_p(ci1,j,k) = -f_p(ci1,j,k)
            f_v(ci1,j,k) = -f_v(ci1,j,k)
          endif
          m_p(ci1,j,k) = fco.cross.f_p(ci1,j,k)
          m_v(ci1,j,k) = fco.cross.f_v(ci1,j,k)
          ! updating global sum (integral)
          fsumpx = fsumpx + f_p(ci1,j,k)%x ; fsumpy = fsumpy + f_p(ci1,j,k)%y ; fsumpz = fsumpz + f_p(ci1,j,k)%z
          fsumvx = fsumvx + f_v(ci1,j,k)%x ; fsumvy = fsumvy + f_v(ci1,j,k)%y ; fsumvz = fsumvz + f_v(ci1,j,k)%z
          msumpx = msumpx + m_p(ci1,j,k)%x ; msumpy = msumpy + m_p(ci1,j,k)%y ; msumpz = msumpz + m_p(ci1,j,k)%z
          msumvx = msumvx + m_v(ci1,j,k)%x ; msumvy = msumvy + m_v(ci1,j,k)%y ; msumvz = msumvz + m_v(ci1,j,k)%z
          Ssum = Ssum + Si(ni1,j,k)
        enddo
      enddo
      !$OMP END PARALLEL
      if (.not.cell) then ! extrapolating ghost cell values
        f_p(ci1,nj1  ,:    ) = f_p(ci1,cj1,:  ) ; f_v(ci1,nj1  ,:    ) = f_v(ci1,cj1,:  )
        f_p(ci1,nj2+1,:    ) = f_p(ci1,nj2,:  ) ; f_v(ci1,nj2+1,:    ) = f_v(ci1,nj2,:  )
        f_p(ci1,:    ,nk1  ) = f_p(ci1,:  ,ck1) ; f_v(ci1,:    ,nk1  ) = f_v(ci1,:  ,ck1)
        f_p(ci1,:    ,nk2+1) = f_p(ci1,:  ,nk2) ; f_v(ci1,:    ,nk2+1) = f_v(ci1,:  ,nk2)
        ! corners
        f_p(ci1,nj1,nk1) = f_p(ci1,cj1,ck1) ; f_v(ci1,nj1,nk1) = f_v(ci1,cj1,ck1)
        f_p(ci1,nj1,nk2) = f_p(ci1,cj1,ck2) ; f_v(ci1,nj1,nk2) = f_v(ci1,cj1,ck2)
        f_p(ci1,nj2,nk1) = f_p(ci1,cj2,ck1) ; f_v(ci1,nj2,nk1) = f_v(ci1,cj2,ck1)
        f_p(ci1,nj2,nk2) = f_p(ci1,cj2,ck2) ; f_v(ci1,nj2,nk2) = f_v(ci1,cj2,ck2)
      endif
    case(3,4)
      !$OMP PARALLEL DEFAULT(NONE)                                                                &
      !$OMP PRIVATE(i,j,k,NdS,hyd)                                                                &
      !$OMP SHARED(ni1,ni2,nj1,nj2,nk1,nk2,ci1,ci2,cj1,cj2,ck1,ck2,face,patch,rFr2,zfs,           &
      !$OMP        node,NFiS,NFjS,NFkS,volume,icc,ticc,f0,momentum,pressure,f_p,f_v,level_set,Re) &
      !$OMP REDUCTION(+: fsumpx,fsumpy,fsumpz,fsumvx,fsumvy,fsumvz,Ssum)
      !$OMP DO
      do k=ck1,ck2
        do i=ci1,ci2
          ! checking if this is an active cell
          if (ticc(i,nj1-1+(face-2),k)/=patch) cycle
          if ( icc(i,cj1,k)/=0) cycle
          if(level_set) then
            if (f0(i,nj1-1+(face-2),k)>0._R_P) cycle
          endif

          ! patch center coordinates
          fco = 0.25_R_P*(node(i,nj1,k)+node(i-1,nj1,k)+node(i,nj1,k-1)+node(i-1,nj1,k-1))

          ! pressure part of forces
          hyd = 0._R_P
          if(level_set) hyd = (fco%z - zfs)*rFr2
          f_p(i,cj1,k) = -(1.5_R_P*pressure(i,cj1,k) - 0.5_R_P*pressure(i,cj1+3-2*(face-2),k) - hyd)*NFjS(i,nj1,k)
          ! viscous part of forces
          NdS = (2._R_P*NFjS(i,nj1,k) + NFjS(i,nj1+1,k) + NFjS(i,nj1-1,k))/(2._R_P*(volume(i,nj1,k)+volume(i,nj1+1,k)))
          f_v(i,cj1,k) = (momentum(i,cj1-1+(face-2),k) - momentum(i,cj1-2+(face-2),k))*(NFjS(i,nj1,k).dot.NdS)/Re

          if (face==4) then
            f_p(i,cj1,k) = -f_p(i,cj1,k)
            f_v(i,cj1,k) = -f_v(i,cj1,k)
          endif
          m_p(i,cj1,k) = fco.cross.f_p(i,cj1,k)
          m_v(i,cj1,k) = fco.cross.f_v(i,cj1,k)
          ! updating global sum (integral)
          fsumpx = fsumpx + f_p(i,cj1,k)%x ; fsumpy = fsumpy + f_p(i,cj1,k)%y ; fsumpz = fsumpz + f_p(i,cj1,k)%z
          fsumvx = fsumvx + f_v(i,cj1,k)%x ; fsumvy = fsumvy + f_v(i,cj1,k)%y ; fsumvz = fsumvz + f_v(i,cj1,k)%z
          msumpx = msumpx + m_p(i,cj1,k)%x ; msumpy = msumpy + m_p(i,cj1,k)%y ; msumpz = msumpz + m_p(i,cj1,k)%z
          msumvx = msumvx + m_v(i,cj1,k)%x ; msumvy = msumvy + m_v(i,cj1,k)%y ; msumvz = msumvz + m_v(i,cj1,k)%z
          Ssum = Ssum + Sj(i,nj1,k)
        enddo
      enddo
      !$OMP END PARALLEL
      if (.not.cell) then ! extrapolating ghost cell values
        f_p(ni1  ,cj1,:    ) = f_p(ci1,cj1,:  ) ; f_v(ni1  ,cj1,:    ) = f_v(ci1,cj1,:  )
        f_p(ni2+1,cj1,:    ) = f_p(ni2,cj1,:  ) ; f_v(ni2+1,cj1,:    ) = f_v(ni2,cj1,:  )
        f_p(:    ,cj1,nk1  ) = f_p(:  ,cj1,ck1) ; f_v(:    ,cj1,nk1  ) = f_v(:  ,cj1,ck1)
        f_p(:    ,cj1,nk2+1) = f_p(:  ,cj1,nk2) ; f_v(:    ,cj1,nk2+1) = f_v(:  ,cj1,nk2)
        ! corners
        f_p(ni1,cj1,nk1) = f_p(ci1,cj1,ck1) ; f_v(ni1,cj1,nk1) = f_v(ci1,cj1,ck1)
        f_p(ni1,cj1,nk2) = f_p(ci1,cj1,ck2) ; f_v(ni1,cj1,nk2) = f_v(ci1,cj1,ck2)
        f_p(ni2,cj1,nk1) = f_p(ci2,cj1,ck1) ; f_v(ni2,cj1,nk1) = f_v(ci2,cj1,ck1)
        f_p(ni2,cj1,nk2) = f_p(ci2,cj1,ck2) ; f_v(ni2,cj1,nk2) = f_v(ci2,cj1,ck2)
      endif
    case(5,6)
      !$OMP PARALLEL DEFAULT(NONE)                                                                &
      !$OMP PRIVATE(i,j,k,NdS,hyd)                                                                &
      !$OMP SHARED(ni1,ni2,nj1,nj2,nk1,nk2,ci1,ci2,cj1,cj2,ck1,ck2,face,patch,rFr2,zfs,           &
      !$OMP        node,NFiS,NFjS,NFkS,volume,icc,ticc,f0,momentum,pressure,f_p,f_v,level_set,Re) &
      !$OMP REDUCTION(+: fsumpx,fsumpy,fsumpz,fsumvx,fsumvy,fsumvz,Ssum)
      !$OMP DO
      do j=cj1,cj2
        do i=ci1,ci2
          ! checking if this is an active cell
          if (ticc(i,j,nk1-1+(face-4))/=patch) cycle
          if ( icc(i,j,ck1)/=0) cycle
          if(level_set) then
            if (f0(i,j,nk1-1+(face-4))>0._R_P) cycle
          endif

          ! patch center coordinates
          fco = 0.25_R_P*(node(i,j,nk1)+node(i-1,j,nk1)+node(i,j-1,nk1)+node(i-1,j-1,nk1))

          ! pressure part of forces
          hyd = 0._R_P
          if(level_set) hyd = (fco%z - zfs)*rFr2
          f_p(i,j,ck1) = -(1.5_R_P*pressure(i,j,ck1) - 0.5_R_P*pressure(i,j,ck1+3-2*(face-4)) - hyd)*NFkS(i,j,nk1)
          ! viscous part of forces
          NdS = (2._R_P*NFkS(i,j,nk1) + NFkS(i,j,nk1+1) + NFkS(i,j,nk1-1))/(2._R_P*(volume(i,j,nk1)+volume(i,j,nk1+1)))
          f_v(i,j,ck1) = (momentum(i,j,ck1-1+(face-4)) - momentum(i,j,ck1-2+(face-4)))*(NFkS(i,j,nk1).dot.NdS)/Re

          if (face==6) then
            f_p(i,j,ck1) = -f_p(i,j,ck1)
            f_v(i,j,ck1) = -f_v(i,j,ck1)
          endif
          m_p(i,j,ck1) = fco.cross.f_p(i,j,ck1)
          m_v(i,j,ck1) = fco.cross.f_v(i,j,ck1)
          ! updating global sum (integral)
          fsumpx = fsumpx + f_p(i,j,ck1)%x ; fsumpy = fsumpy + f_p(i,j,ck1)%y ; fsumpz = fsumpz + f_p(i,j,ck1)%z
          fsumvx = fsumvx + f_v(i,j,ck1)%x ; fsumvy = fsumvy + f_v(i,j,ck1)%y ; fsumvz = fsumvz + f_v(i,j,ck1)%z
          msumpx = msumpx + m_p(i,j,ck1)%x ; msumpy = msumpy + m_p(i,j,ck1)%y ; msumpz = msumpz + m_p(i,j,ck1)%z
          msumvx = msumvx + m_v(i,j,ck1)%x ; msumvy = msumvy + m_v(i,j,ck1)%y ; msumvz = msumvz + m_v(i,j,ck1)%z
          Ssum = Ssum + Sk(i,j,nk1)
        enddo
      enddo
      !$OMP END PARALLEL
      if (.not.cell) then ! extrapolating ghost cell values
        f_p(ni1  ,:    ,ck1) = f_p(ci1,:  ,ck1) ; f_v(ni1  ,:    ,ck1) = f_v(ci1,:  ,ck1)
        f_p(ni2+1,:    ,ck1) = f_p(ni2,:  ,ck1) ; f_v(ni2+1,:    ,ck1) = f_v(ni2,:  ,ck1)
        f_p(:    ,nj1  ,ck1) = f_p(:  ,cj1,ck1) ; f_v(:    ,nj1  ,ck1) = f_v(:  ,cj1,ck1)
        f_p(:    ,nj2+1,ck1) = f_p(:  ,nj2,ck1) ; f_v(:    ,nj2+1,ck1) = f_v(:  ,nj2,ck1)
        ! corners
        f_p(ni1,nj1,ck1) = f_p(ci1,cj1,ck1) ; f_v(ni1,nj1,ck1) = f_v(ci1,cj1,ck1)
        f_p(ni1,nj2,ck1) = f_p(ci1,cj2,ck1) ; f_v(ni1,nj2,ck1) = f_v(ci1,cj2,ck1)
        f_p(ni2,nj1,ck1) = f_p(ci2,cj1,ck1) ; f_v(ni2,nj1,ck1) = f_v(ci2,cj1,ck1)
        f_p(ni2,nj2,ck1) = f_p(ci2,cj2,ck1) ; f_v(ni2,nj2,ck1) = f_v(ci2,cj2,ck1)
      endif
    endselect
    fsum_p = fsum_p + (ex*fsumpx+ey*fsumpy+ez*fsumpz)
    fsum_v = fsum_v + (ex*fsumvx+ey*fsumvy+ez*fsumvz)
    msum_p = msum_p + (ex*msumpx+ey*msumpy+ez*msumpz)
    msum_v = msum_v + (ex*msumvx+ey*msumvy+ez*msumvz)
    ! saving on the Riccardo Broglia forces.RB scratch file
    write(unit_for_RB_scr)b,face,(ci2-ci1+1),(cj2-cj1+1),(ck2-ck1+1)
    write(unit_for_RB_scr)f_p(ci1:ci2,cj1:cj2,ck1:ck2)%x
    write(unit_for_RB_scr)f_p(ci1:ci2,cj1:cj2,ck1:ck2)%y
    write(unit_for_RB_scr)f_p(ci1:ci2,cj1:cj2,ck1:ck2)%z
    write(unit_for_RB_scr)f_v(ci1:ci2,cj1:cj2,ck1:ck2)%x
    write(unit_for_RB_scr)f_v(ci1:ci2,cj1:cj2,ck1:ck2)%y
    write(unit_for_RB_scr)f_v(ci1:ci2,cj1:cj2,ck1:ck2)%z
    write(unit_for_RB_scr)m_p(ci1:ci2,cj1:cj2,ck1:ck2)%x
    write(unit_for_RB_scr)m_p(ci1:ci2,cj1:cj2,ck1:ck2)%y
    write(unit_for_RB_scr)m_p(ci1:ci2,cj1:cj2,ck1:ck2)%z
    write(unit_for_RB_scr)m_v(ci1:ci2,cj1:cj2,ck1:ck2)%x
    write(unit_for_RB_scr)m_v(ci1:ci2,cj1:cj2,ck1:ck2)%y
    write(unit_for_RB_scr)m_v(ci1:ci2,cj1:cj2,ck1:ck2)%z
    ! saving on the forces-grp(#) file
    if (Ng>0) then
       write(unit_g_for(blockmap(1,b)),'(A)',iostat=err)trim(str(n=(fsumvx+fsumpx)))//' '// &
                                                        trim(str(n=(fsumvy+fsumpy)))//' '// &
                                                        trim(str(n=(fsumvz+fsumpz)))//' '// &
                                                        trim(str(n=(       fsumpx)))//' '// &
                                                        trim(str(n=(       fsumpy)))//' '// &
                                                        trim(str(n=(       fsumpz)))//' '// &
                                                        trim(str(n=(fsumvx       )))//' '// &
                                                        trim(str(n=(fsumvy       )))//' '// &
                                                        trim(str(n=(fsumvz       )))//' '// &
                                                        trim(str(n=(msumvx+msumpx)))//' '// &
                                                        trim(str(n=(msumvy+msumpy)))//' '// &
                                                        trim(str(n=(msumvz+msumpz)))//' '// &
                                                        trim(str(n=(       msumpx)))//' '// &
                                                        trim(str(n=(       msumpy)))//' '// &
                                                        trim(str(n=(       msumpz)))//' '// &
                                                        trim(str(n=(msumvx       )))//' '// &
                                                        trim(str(n=(msumvy       )))//' '// &
                                                        trim(str(n=(msumvz       )))//' '// &
      'Fx,Fy,Fz,Fx_p,Fy_p,Fz_p,Fx_v,Fy_v,Fz_v,Mx,My,Mz,Mx_p,My_p,Mz_p,Mx_v,My_v,Mz_v"'
    endif
  endif
  if (yp) then ! computing y+ on the patch
    select case(face)
    case(1,2)
      !$OMP PARALLEL DEFAULT(NONE) &
      !$OMP PRIVATE(i,j,k,NdS)     &
      !$OMP SHARED(ni1,ni2,nj1,nj2,nk1,nk2,ci1,ci2,cj1,cj2,ck1,ck2,face,patch,node,NFiS,icc,ticc,momentum,f0,Re,yplus)
      !$OMP DO
      do k=ck1,ck2
        do j=cj1,cj2
          ! checking if this is an active cell
          if (ticc(ni1-1+face,j,k)/=patch) cycle
          if ( icc(ci1,j,k)/=0) cycle
          if(level_set) then
            if (f0(ni1-1+face,j,k)>0._R_P) cycle
          endif

          ! distance from the patch
          NdS = 0.25_R_P*(node(ni1,j,k)+node(ni1,j-1,k)+node(ni1,j,k-1)+node(ni1,j-1,k-1)) - &
                0.25_R_P*(node(ci1,j,k)+node(ci1,j-1,k)+node(ci1,j,k-1)+node(ci1,j-1,k-1))
          yplus(ci1,j,k) = sqrt(Re*normL2(momentum(ci1,j,k).ortho.NFiS(ni1,j,k))*(abs(NdS.dot.NFiS(ni1,j,k))))
        enddo
      enddo
      !$OMP END PARALLEL
      if (.not.cell) then ! extrapolating ghost cell values
        yplus(ci1,nj1  ,:    ) = yplus(ci1,cj1,:  )
        yplus(ci1,nj2+1,:    ) = yplus(ci1,nj2,:  )
        yplus(ci1,:    ,nk1  ) = yplus(ci1,:  ,ck1)
        yplus(ci1,:    ,nk2+1) = yplus(ci1,:  ,nk2)
      endif
    case(3,4)
      !$OMP PARALLEL DEFAULT(NONE) &
      !$OMP PRIVATE(i,j,k,NdS)     &
      !$OMP SHARED(ni1,ni2,nj1,nj2,nk1,nk2,ci1,ci2,cj1,cj2,ck1,ck2,face,patch,node,NFjS,icc,ticc,momentum,f0,Re,yplus)
      !$OMP DO
      do k=ck1,ck2
        do i=ci1,ci2
          ! checking if this is an active cell
          if (ticc(i,nj1-1+(face-2),k)/=patch) cycle
          if ( icc(i,cj1,k)/=0) cycle
          if(level_set) then
            if (f0(i,nj1-1+(face-2),k)>0._R_P) cycle
          endif

          ! distance from the patch
          NdS = 0.25_R_P*(node(i,nj1,k)+node(i-1,nj1,k)+node(i,nj1,k-1)+node(i-1,nj1,k-1)) - &
                0.25_R_P*(node(i,cj1,k)+node(i-1,cj1,k)+node(i,cj1,k-1)+node(i-1,cj1,k-1))
          yplus(i,cj1,k) = sqrt(Re*normL2(momentum(i,cj1,k).ortho.NFjS(i,nj1,k))*(abs(NdS.dot.NFjS(i,nj1,k))))
        enddo
      enddo
      !$OMP END PARALLEL
      if (.not.cell) then ! extrapolating ghost cell values
        yplus(ni1  ,cj1,:    ) = yplus(ci1,cj1,:  )
        yplus(ni2+1,cj1,:    ) = yplus(ni2,cj1,:  )
        yplus(:    ,cj1,nk1  ) = yplus(:  ,cj1,ck1)
        yplus(:    ,cj1,nk2+1) = yplus(:  ,cj1,nk2)
      endif
    case(5,6)
      !$OMP PARALLEL DEFAULT(NONE) &
      !$OMP PRIVATE(i,j,k,NdS)     &
      !$OMP SHARED(ni1,ni2,nj1,nj2,nk1,nk2,ci1,ci2,cj1,cj2,ck1,ck2,face,patch,node,NFkS,icc,ticc,momentum,f0,Re,yplus)
      !$OMP DO
      do j=cj1,cj2
        do i=ci1,ci2
          ! checking if this is an active cell
          if (ticc(i,j,nk1-1+(face-4))/=patch) cycle
          if ( icc(i,j,ck1)/=0) cycle
          if(level_set) then
            if (f0(i,j,nk1-1+(face-4))>0._R_P) cycle
          endif

          ! distance from the patch
          NdS = 0.25_R_P*(node(i,j,nk1)+node(i-1,j,nk1)+node(i,j-1,nk1)+node(i-1,j-1,nk1)) - &
                0.25_R_P*(node(i,j,ck1)+node(i-1,j,ck1)+node(i,j-1,ck1)+node(i-1,j-1,ck1))
          yplus(i,j,ck1) = sqrt(Re*normL2(momentum(i,j,ck1).ortho.NFkS(i,j,nk1))*(abs(NdS.dot.NFkS(i,j,nk1))))
        enddo
      enddo
      !$OMP END PARALLEL
      if (.not.cell) then ! extrapolating ghost cell values
        yplus(ni1  ,:    ,ck1) = yplus(ci1,:  ,ck1)
        yplus(ni2+1,:    ,ck1) = yplus(ni2,:  ,ck1)
        yplus(:    ,nj1  ,ck1) = yplus(:  ,cj1,ck1)
        yplus(:    ,nj2+1,ck1) = yplus(:  ,nj2,ck1)
      endif
    endselect
  endif
  nnode = (ni2-ni1+1)*(nj2-nj1+1)*(nk2-nk1+1) ! computing number of nodes
  ncell = (ci2-ci1+1)*(cj2-cj1+1)*(ck2-ck1+1) ! computing number of cells
  ! patch geometry
  if (tec) then
    if (binary) then
      tecnull(1:Nvar)= 0
      tecvarloc = 1
      if (cell) tecvarloc(4:Nvar)= 0
      err = teczne112('ptc_'//trim(strz(3,patch))//'_n_'//trim(strz(5,Np))// &
                      '_blk_'//trim(strz(4,blockmap(2,b)))//                 &
                      '_fac_'//trim(strz(1,face))//                          &
                      '_grp_'//trim(strz(3,blockmap(1,b)))//tecendrec,       &
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
                      1,                                                     & !1=>block,0=>point
                      0,                                                     &
                      0,                                                     &
                      0,                                                     &
                      0,                                                     &
                      0,                                                     &
                      tecnull,                                               &
                      tecvarloc,                                             &
                      tecnull,                                               &
                      0)
    else
      if (cell) then
        write(unit_out,'(A)',iostat=err)' ZONE  T="ptc_'//trim(strz(3,patch))//'_n_'//trim(strz(5,Np))// &
                                        '_blk_'//trim(strz(4,blockmap(2,b)))//                           &
                                        '_fac_'//trim(strz(1,face))//                                    &
                                        '_grp_'//trim(strz(3,blockmap(1,b)))//'"'//                      &
                                        ', I='//trim(str(no_sign=.true.,n=ni2-ni1+1))//                  &
                                        ', J='//trim(str(no_sign=.true.,n=nj2-nj1+1))//                  &
                                        ', K='//trim(str(no_sign=.true.,n=nk2-nk1+1))//                  &
                                        ', DATAPACKING=BLOCK'//                                          &
                                        ', VARLOCATION=([1-3]=NODAL,[4-'//trim(str(.true.,Nvar))//']=CELLCENTERED)'
      else
        write(unit_out,'(A)',iostat=err)' ZONE  T="ptc_'//trim(strz(3,patch))//'_n_'//trim(strz(5,Np))// &
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
    err = tec_var(n=nnode,var=node(ni1:ni2,nj1:nj2,nk1:nk2)%x,d=1)
    err = tec_var(n=nnode,var=node(ni1:ni2,nj1:nj2,nk1:nk2)%y,d=1)
    err = tec_var(n=nnode,var=node(ni1:ni2,nj1:nj2,nk1:nk2)%z,d=1)
  endif
  if (vtk) then
    if (binary) then
      err = VTK_INI_XML(output_format = 'binary',                                                                      &
                        filename      = adjustl(trim(vtkbfile))//".ptc"//trim(strz(2,patch))//'.n'//trim(strz(5,Np))// &
                                        '_blk_'//trim(strz(4,blockmap(2,b)))//                                         &
                                        '_fac_'//trim(strz(1,face))//                                                  &
                                        '_grp_'//trim(strz(3,blockmap(1,b)))//'.vts',                                  &
                        mesh_topology = 'StructuredGrid',                                                              &
                        nx1 = ni1, nx2 = ni2,                                                                          &
                        ny1 = nj1, ny2 = nj2,                                                                          &
                        nz1 = nk1, nz2 = nk2)
    else
      err = VTK_INI_XML(output_format = 'ascii',                                                                       &
                        filename      = adjustl(trim(vtkbfile))//".ptc"//trim(strz(2,patch))//'.n'//trim(strz(5,Np))// &
                                        '_blk_'//trim(strz(4,blockmap(2,b)))//                                         &
                                        '_fac_'//trim(strz(1,face))//                                                  &
                                        '_grp_'//trim(strz(3,blockmap(1,b)))//'.vts',                                  &
                        mesh_topology = 'StructuredGrid',                                                              &
                        nx1 = ni1, nx2 = ni2,                                                                          &
                        ny1 = nj1, ny2 = nj2,                                                                          &
                        nz1 = nk1, nz2 = nk2)
    endif
    err = VTK_GEO_XML(nx1 = ni1, nx2 = ni2,                                 &
                      ny1 = nj1, ny2 = nj2,                                 &
                      nz1 = nk1, nz2 = nk2,                                 &
                      NN = nnode,                                           &
                      X=reshape(node(ni1:ni2,nj1:nj2,nk1:nk2)%x,(/nnode/)), &
                      Y=reshape(node(ni1:ni2,nj1:nj2,nk1:nk2)%y,(/nnode/)), &
                      Z=reshape(node(ni1:ni2,nj1:nj2,nk1:nk2)%z,(/nnode/)))
  endif
  ! patch variables
  ! icc
  if (cell) then
    if (tec) err = tec_var(n=ncell,var=real(ricc(ci1:ci2,cj1:cj2,ck1:ck2),R_P), d=1)
    if (vtk) err = VTK_DAT_XML(var_location = 'cell', var_block_action = 'open')
    if (vtk) err = VTK_VAR_XML(NC_NN=ncell,varname='icc',var=reshape(ricc(ci1:ci2,cj1:cj2,ck1:ck2),(/ncell/)))
  else
    if (tec) err = tec_var(n=nnode,var=real(vicc(ni1:ni2,nj1:nj2,nk1:nk2),R_P), d=1)
    if (vtk) err = VTK_DAT_XML(var_location = 'node', var_block_action = 'open')
    if (vtk) err = VTK_VAR_XML(NC_NN=nnode,varname='icc',var=reshape(vicc(ni1:ni2,nj1:nj2,nk1:nk2),(/nnode/)))
  endif
  ! solution variables
  if (sol) then
    if (cell) then
      if (tec) err = tec_var(n=ncell,var=momentum(ci1:ci2,cj1:cj2,ck1:ck2)%x,d=1)
      if (tec) err = tec_var(n=ncell,var=momentum(ci1:ci2,cj1:cj2,ck1:ck2)%y,d=1)
      if (tec) err = tec_var(n=ncell,var=momentum(ci1:ci2,cj1:cj2,ck1:ck2)%z,d=1)
      if (tec) err = tec_var(n=ncell,var=pressure(ci1:ci2,cj1:cj2,ck1:ck2)  ,d=1)
      if (vtk) err = VTK_VAR_XML(NC_NN=ncell,varname='u',var=reshape(momentum(ci1:ci2,cj1:cj2,ck1:ck2)%X,(/ncell/)))
      if (vtk) err = VTK_VAR_XML(NC_NN=ncell,varname='v',var=reshape(momentum(ci1:ci2,cj1:cj2,ck1:ck2)%y,(/ncell/)))
      if (vtk) err = VTK_VAR_XML(NC_NN=ncell,varname='w',var=reshape(momentum(ci1:ci2,cj1:cj2,ck1:ck2)%z,(/ncell/)))
      if (vtk) err = VTK_VAR_XML(NC_NN=ncell,varname='p',var=reshape(pressure(ci1:ci2,cj1:cj2,ck1:ck2)  ,(/ncell/)))
      if (zeroeq) then
        if (tec) err = tec_var(n=ncell,var=visc(ci1:ci2,cj1:cj2,ck1:ck2),d=1)
         if (vtk)  err = VTK_VAR_XML(NC_NN=ncell,varname='visc',var=reshape(visc(ci1:ci2,cj1:cj2,ck1:ck2),(/ncell/)))
      elseif (oneeq) then
        if (tec) err = tec_var(n=ncell,var=visc(ci1:ci2,cj1:cj2,ck1:ck2),d=1)
        if (tec) err = tec_var(n=ncell,var=vitl(ci1:ci2,cj1:cj2,ck1:ck2),d=1)
        if (vtk) err = VTK_VAR_XML(NC_NN=ncell,varname='visc',var=reshape(visc(ci1:ci2,cj1:cj2,ck1:ck2),(/ncell/)))
        if (vtk) err = VTK_VAR_XML(NC_NN=ncell,varname='vitl',var=reshape(visc(ci1:ci2,cj1:cj2,ck1:ck2),(/ncell/)))
      elseif (twoeq) then
        if (tec) err = tec_var(n=ncell,var=visc(ci1:ci2,cj1:cj2,ck1:ck2),d=1)
        if (tec) err = tec_var(n=ncell,var=ken (ci1:ci2,cj1:cj2,ck1:ck2),d=1)
        if (tec) err = tec_var(n=ncell,var=eps (ci1:ci2,cj1:cj2,ck1:ck2),d=1)
        if (vtk) err = VTK_VAR_XML(NC_NN=ncell,varname='visc',var=reshape(visc(ci1:ci2,cj1:cj2,ck1:ck2),(/ncell/)))
        if (vtk) err = VTK_VAR_XML(NC_NN=ncell,varname='ken', var=reshape(ken (ci1:ci2,cj1:cj2,ck1:ck2),(/ncell/)))
        if (vtk) err = VTK_VAR_XML(NC_NN=ncell,varname='eps', var=reshape(eps (ci1:ci2,cj1:cj2,ck1:ck2),(/ncell/)))
      endif
      if (level_set) then
        if (tec) err = tec_var(n=ncell,var=f (ci1:ci2,cj1:cj2,ck1:ck2),d=1)
        if (tec) err = tec_var(n=ncell,var=f0(ci1:ci2,cj1:cj2,ck1:ck2),d=1)
        if (vtk) err = VTK_VAR_XML(NC_NN=ncell,varname='f', var=reshape(f (ci1:ci2,cj1:cj2,ck1:ck2),(/ncell/)))
        if (vtk) err = VTK_VAR_XML(NC_NN=ncell,varname='f0',var=reshape(f0(ci1:ci2,cj1:cj2,ck1:ck2),(/ncell/)))
      endif
      if (forces) then
        if (tec) err = tec_var(n=ncell,var=(f_p(ci1:ci2,cj1:cj2,ck1:ck2)%x + &
                                            f_v(ci1:ci2,cj1:cj2,ck1:ck2)%x),d=1)
        if (tec) err = tec_var(n=ncell,var=(f_p(ci1:ci2,cj1:cj2,ck1:ck2)%y + &
                                            f_v(ci1:ci2,cj1:cj2,ck1:ck2)%y),d=1)
        if (tec) err = tec_var(n=ncell,var=(f_p(ci1:ci2,cj1:cj2,ck1:ck2)%z + &
                                            f_v(ci1:ci2,cj1:cj2,ck1:ck2)%z),d=1)
        if (tec) err = tec_var(n=ncell,var= f_p(ci1:ci2,cj1:cj2,ck1:ck2)%x,d=1)
        if (tec) err = tec_var(n=ncell,var= f_p(ci1:ci2,cj1:cj2,ck1:ck2)%y,d=1)
        if (tec) err = tec_var(n=ncell,var= f_p(ci1:ci2,cj1:cj2,ck1:ck2)%z,d=1)
        if (tec) err = tec_var(n=ncell,var= f_v(ci1:ci2,cj1:cj2,ck1:ck2)%x,d=1)
        if (tec) err = tec_var(n=ncell,var= f_v(ci1:ci2,cj1:cj2,ck1:ck2)%y,d=1)
        if (tec) err = tec_var(n=ncell,var= f_v(ci1:ci2,cj1:cj2,ck1:ck2)%z,d=1)
        if (tec) err = tec_var(n=ncell,var=(m_p(ci1:ci2,cj1:cj2,ck1:ck2)%x + &
                                            m_v(ci1:ci2,cj1:cj2,ck1:ck2)%x),d=1)
        if (tec) err = tec_var(n=ncell,var=(m_p(ci1:ci2,cj1:cj2,ck1:ck2)%y + &
                                            m_v(ci1:ci2,cj1:cj2,ck1:ck2)%y),d=1)
        if (tec) err = tec_var(n=ncell,var=(m_p(ci1:ci2,cj1:cj2,ck1:ck2)%z + &
                                            m_v(ci1:ci2,cj1:cj2,ck1:ck2)%z),d=1)
        if (tec) err = tec_var(n=ncell,var= m_p(ci1:ci2,cj1:cj2,ck1:ck2)%x,d=1)
        if (tec) err = tec_var(n=ncell,var= m_p(ci1:ci2,cj1:cj2,ck1:ck2)%y,d=1)
        if (tec) err = tec_var(n=ncell,var= m_p(ci1:ci2,cj1:cj2,ck1:ck2)%z,d=1)
        if (tec) err = tec_var(n=ncell,var= m_v(ci1:ci2,cj1:cj2,ck1:ck2)%x,d=1)
        if (tec) err = tec_var(n=ncell,var= m_v(ci1:ci2,cj1:cj2,ck1:ck2)%y,d=1)
        if (tec) err = tec_var(n=ncell,var= m_v(ci1:ci2,cj1:cj2,ck1:ck2)%z,d=1)
        if (vtk) err = VTK_VAR_XML(NC_NN=ncell,varname='fx',  var=reshape((f_p(ci1:ci2,cj1:cj2,ck1:ck2)%x + &
                                                                           f_v(ci1:ci2,cj1:cj2,ck1:ck2)%x),(/ncell/)))
        if (vtk) err = VTK_VAR_XML(NC_NN=ncell,varname='fy',  var=reshape((f_p(ci1:ci2,cj1:cj2,ck1:ck2)%y + &
                                                                           f_v(ci1:ci2,cj1:cj2,ck1:ck2)%y),(/ncell/)))
        if (vtk) err = VTK_VAR_XML(NC_NN=ncell,varname='fz',  var=reshape((f_p(ci1:ci2,cj1:cj2,ck1:ck2)%z + &
                                                                           f_v(ci1:ci2,cj1:cj2,ck1:ck2)%z),(/ncell/)))
        if (vtk) err = VTK_VAR_XML(NC_NN=ncell,varname='fx_p',var=reshape( f_p(ci1:ci2,cj1:cj2,ck1:ck2)%x,(/ncell/)))
        if (vtk) err = VTK_VAR_XML(NC_NN=ncell,varname='fy_p',var=reshape( f_p(ci1:ci2,cj1:cj2,ck1:ck2)%y,(/ncell/)))
        if (vtk) err = VTK_VAR_XML(NC_NN=ncell,varname='fz_p',var=reshape( f_p(ci1:ci2,cj1:cj2,ck1:ck2)%z,(/ncell/)))
        if (vtk) err = VTK_VAR_XML(NC_NN=ncell,varname='fx_v',var=reshape( f_v(ci1:ci2,cj1:cj2,ck1:ck2)%x,(/ncell/)))
        if (vtk) err = VTK_VAR_XML(NC_NN=ncell,varname='fy_v',var=reshape( f_v(ci1:ci2,cj1:cj2,ck1:ck2)%y,(/ncell/)))
        if (vtk) err = VTK_VAR_XML(NC_NN=ncell,varname='fz_v',var=reshape( f_v(ci1:ci2,cj1:cj2,ck1:ck2)%z,(/ncell/)))
        if (vtk) err = VTK_VAR_XML(NC_NN=ncell,varname='fx',  var=reshape((m_p(ci1:ci2,cj1:cj2,ck1:ck2)%x + &
                                                                           m_v(ci1:ci2,cj1:cj2,ck1:ck2)%x),(/ncell/)))
        if (vtk) err = VTK_VAR_XML(NC_NN=ncell,varname='fy',  var=reshape((m_p(ci1:ci2,cj1:cj2,ck1:ck2)%y + &
                                                                           m_v(ci1:ci2,cj1:cj2,ck1:ck2)%y),(/ncell/)))
        if (vtk) err = VTK_VAR_XML(NC_NN=ncell,varname='fz',  var=reshape((m_p(ci1:ci2,cj1:cj2,ck1:ck2)%z + &
                                                                           m_v(ci1:ci2,cj1:cj2,ck1:ck2)%z),(/ncell/)))
        if (vtk) err = VTK_VAR_XML(NC_NN=ncell,varname='fx_p',var=reshape( m_p(ci1:ci2,cj1:cj2,ck1:ck2)%x,(/ncell/)))
        if (vtk) err = VTK_VAR_XML(NC_NN=ncell,varname='fy_p',var=reshape( m_p(ci1:ci2,cj1:cj2,ck1:ck2)%y,(/ncell/)))
        if (vtk) err = VTK_VAR_XML(NC_NN=ncell,varname='fz_p',var=reshape( m_p(ci1:ci2,cj1:cj2,ck1:ck2)%z,(/ncell/)))
        if (vtk) err = VTK_VAR_XML(NC_NN=ncell,varname='fx_v',var=reshape( m_v(ci1:ci2,cj1:cj2,ck1:ck2)%x,(/ncell/)))
        if (vtk) err = VTK_VAR_XML(NC_NN=ncell,varname='fy_v',var=reshape( m_v(ci1:ci2,cj1:cj2,ck1:ck2)%y,(/ncell/)))
        if (vtk) err = VTK_VAR_XML(NC_NN=ncell,varname='fz_v',var=reshape( m_v(ci1:ci2,cj1:cj2,ck1:ck2)%z,(/ncell/)))
      endif
      if (yp) then
        if (tec) err = tec_var(n=ncell,var= yplus(ni1:ni2,cj1:cj2,ck1:ck2),d=1)
        if (vtk) err = VTK_VAR_XML(NC_NN=ncell,varname='yplus',var=reshape(yplus(ci1:ci2,cj1:cj2,ck1:ck2),(/ncell/)))
      endif
    else
      call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                            var  = momentum(ni1:ni2+1,nj1:nj2+1,nk1:nk2+1)%x,         &
                            vari = vari(    ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
      if (tec) err = tec_var(n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
      if (vtk) err = VTK_VAR_XML(NC_NN=nnode,varname='u',var=reshape(vari(ni1:ni2,nj1:nj2,nk1:nk2),(/nnode/)))
      call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                            var  = momentum(ni1:ni2+1,nj1:nj2+1,nk1:nk2+1)%y,         &
                            vari = vari(    ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
      if (tec) err = tec_var(n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
      if (vtk) err = VTK_VAR_XML(NC_NN=nnode,varname='v',var=reshape(vari(ni1:ni2,nj1:nj2,nk1:nk2),(/nnode/)))
      call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                            var  = momentum(ni1:ni2+1,nj1:nj2+1,nk1:nk2+1)%z,         &
                            vari = vari(    ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
      if (tec) err = tec_var(n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
      if (vtk) err = VTK_VAR_XML(NC_NN=nnode,varname='w',var=reshape(vari(ni1:ni2,nj1:nj2,nk1:nk2),(/nnode/)))
      call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                            var  = pressure(ni1:ni2+1,nj1:nj2+1,nk1:nk2+1),           &
                            vari = vari(    ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
      if (tec) err = tec_var(n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
      if (vtk) err = VTK_VAR_XML(NC_NN=nnode,varname='p',var=reshape(vari(ni1:ni2,nj1:nj2,nk1:nk2),(/nnode/)))
      if (zeroeq) then
        call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                              var  = visc(ni1:ni2+1,nj1:nj2+1,nk1:nk2+1),               &
                              vari = vari(ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
        if (tec) err = tec_var(n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
        if (vtk) err = VTK_VAR_XML(NC_NN=nnode,varname='visc',var=reshape(vari(ni1:ni2,nj1:nj2,nk1:nk2),(/nnode/)))
      elseif (oneeq) then
        call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                              var  = visc(ni1:ni2+1,nj1:nj2+1,nk1:nk2+1),               &
                              vari = vari(ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
        if (tec) err = tec_var(n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
        if (vtk) err = VTK_VAR_XML(NC_NN=nnode,varname='visc',var=reshape(vari(ni1:ni2,nj1:nj2,nk1:nk2),(/nnode/)))
        call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                              var  = vitl(ni1:ni2+1,nj1:nj2+1,nk1:nk2+1),               &
                              vari = vari(ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
        if (tec) err = tec_var(n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
        if (vtk) err = VTK_VAR_XML(NC_NN=nnode,varname='vitl',var=reshape(vari(ni1:ni2,nj1:nj2,nk1:nk2),(/nnode/)))
      elseif (twoeq) then
        call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                              var  = visc(ni1:ni2+1,nj1:nj2+1,nk1:nk2+1),               &
                              vari = vari(ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
        if (tec) err = tec_var(n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
        if (vtk) err = VTK_VAR_XML(NC_NN=nnode,varname='visc',var=reshape(vari(ni1:ni2,nj1:nj2,nk1:nk2),(/nnode/)))
        call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                              var  = ken (ni1:ni2+1,nj1:nj2+1,nk1:nk2+1),               &
                              vari = vari(ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
        if (tec) err = tec_var(n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
        if (vtk) err = VTK_VAR_XML(NC_NN=nnode,varname='ken',var=reshape(vari(ni1:ni2,nj1:nj2,nk1:nk2),(/nnode/)))
        call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                              var  = eps (ni1:ni2+1,nj1:nj2+1,nk1:nk2+1),               &
                              vari = vari(ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
        if (tec) err = tec_var(n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
        if (vtk) err = VTK_VAR_XML(NC_NN=nnode,varname='eps',var=reshape(vari(ni1:ni2,nj1:nj2,nk1:nk2),(/nnode/)))
      endif
      if (level_set) then
        call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                              var  = f   (ni1:ni2+1,nj1:nj2+1,nk1:nk2+1),               &
                              vari = vari(ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
        if (tec) err = tec_var(n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
        if (vtk) err = VTK_VAR_XML(NC_NN=nnode,varname='f',var=reshape(vari(ni1:ni2,nj1:nj2,nk1:nk2),(/nnode/)))
        call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                              var  = f0  (ni1:ni2+1,nj1:nj2+1,nk1:nk2+1),               &
                              vari = vari(ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
        if (tec) err = tec_var(n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
        if (vtk) err = VTK_VAR_XML(NC_NN=nnode,varname='f0',var=reshape(vari(ni1:ni2,nj1:nj2,nk1:nk2),(/nnode/)))
      endif
      if (forces) then
        call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                              var  = f_p (ni1:ni2+1,nj1:nj2+1,nk1:nk2+1)%x +            &
                                     f_v (ni1:ni2+1,nj1:nj2+1,nk1:nk2+1)%x,             &
                              vari = vari(ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
        if (tec) err = tec_var(n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
        if (vtk) err = VTK_VAR_XML(NC_NN=nnode,varname='fx',var=reshape(vari(ni1:ni2,nj1:nj2,nk1:nk2),(/nnode/)))
        call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                              var  = f_p (ni1:ni2+1,nj1:nj2+1,nk1:nk2+1)%y +            &
                                     f_v (ni1:ni2+1,nj1:nj2+1,nk1:nk2+1)%y,             &
                              vari = vari(ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
        if (tec) err = tec_var(n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
        if (vtk) err = VTK_VAR_XML(NC_NN=nnode,varname='fy',var=reshape(vari(ni1:ni2,nj1:nj2,nk1:nk2),(/nnode/)))
        call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                              var  = f_p (ni1:ni2+1,nj1:nj2+1,nk1:nk2+1)%z +            &
                                     f_v (ni1:ni2+1,nj1:nj2+1,nk1:nk2+1)%z,             &
                              vari = vari(ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
        if (tec) err = tec_var(n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
        if (vtk) err = VTK_VAR_XML(NC_NN=nnode,varname='fz',var=reshape(vari(ni1:ni2,nj1:nj2,nk1:nk2),(/nnode/)))
        call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                              var  = f_p (ni1:ni2+1,nj1:nj2+1,nk1:nk2+1)%x,             &
                              vari = vari(ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
        if (tec) err = tec_var(n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
        if (vtk) err = VTK_VAR_XML(NC_NN=nnode,varname='fx_p',var=reshape(vari(ni1:ni2,nj1:nj2,nk1:nk2),(/nnode/)))
        call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                              var  = f_p (ni1:ni2+1,nj1:nj2+1,nk1:nk2+1)%y,             &
                              vari = vari(ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
        if (tec) err = tec_var(n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
        if (vtk) err = VTK_VAR_XML(NC_NN=nnode,varname='fy_p',var=reshape(vari(ni1:ni2,nj1:nj2,nk1:nk2),(/nnode/)))
        call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                              var  = f_p (ni1:ni2+1,nj1:nj2+1,nk1:nk2+1)%z,             &
                              vari = vari(ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
        if (tec) err = tec_var(n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
        if (vtk) err = VTK_VAR_XML(NC_NN=nnode,varname='fz_p',var=reshape(vari(ni1:ni2,nj1:nj2,nk1:nk2),(/nnode/)))
        call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                              var  = f_v (ni1:ni2+1,nj1:nj2+1,nk1:nk2+1)%x,             &
                              vari = vari(ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
        if (tec) err = tec_var(n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
        if (vtk) err = VTK_VAR_XML(NC_NN=nnode,varname='fx_v',var=reshape(vari(ni1:ni2,nj1:nj2,nk1:nk2),(/nnode/)))
        call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                              var  = f_v (ni1:ni2+1,nj1:nj2+1,nk1:nk2+1)%y,             &
                              vari = vari(ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
        if (tec) err = tec_var(n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
        if (vtk) err = VTK_VAR_XML(NC_NN=nnode,varname='fy_v',var=reshape(vari(ni1:ni2,nj1:nj2,nk1:nk2),(/nnode/)))
        call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                              var  = f_v (ni1:ni2+1,nj1:nj2+1,nk1:nk2+1)%z,             &
                              vari = vari(ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
        if (tec) err = tec_var(n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
        if (vtk) err = VTK_VAR_XML(NC_NN=nnode,varname='fz_v',var=reshape(vari(ni1:ni2,nj1:nj2,nk1:nk2),(/nnode/)))
        call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                              var  = m_p (ni1:ni2+1,nj1:nj2+1,nk1:nk2+1)%x +            &
                                     m_v (ni1:ni2+1,nj1:nj2+1,nk1:nk2+1)%x,             &
                              vari = vari(ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
        if (tec) err = tec_var(n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
        if (vtk) err = VTK_VAR_XML(NC_NN=nnode,varname='mx',var=reshape(vari(ni1:ni2,nj1:nj2,nk1:nk2),(/nnode/)))
        call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                              var  = m_p (ni1:ni2+1,nj1:nj2+1,nk1:nk2+1)%y +            &
                                     m_v (ni1:ni2+1,nj1:nj2+1,nk1:nk2+1)%y,             &
                              vari = vari(ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
        if (tec) err = tec_var(n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
        if (vtk) err = VTK_VAR_XML(NC_NN=nnode,varname='my',var=reshape(vari(ni1:ni2,nj1:nj2,nk1:nk2),(/nnode/)))
        call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                              var  = m_p (ni1:ni2+1,nj1:nj2+1,nk1:nk2+1)%z +            &
                                     m_v (ni1:ni2+1,nj1:nj2+1,nk1:nk2+1)%z,             &
                              vari = vari(ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
        if (tec) err = tec_var(n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
        if (vtk) err = VTK_VAR_XML(NC_NN=nnode,varname='mz',var=reshape(vari(ni1:ni2,nj1:nj2,nk1:nk2),(/nnode/)))
        call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                              var  = m_p (ni1:ni2+1,nj1:nj2+1,nk1:nk2+1)%x,             &
                              vari = vari(ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
        if (tec) err = tec_var(n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
        if (vtk) err = VTK_VAR_XML(NC_NN=nnode,varname='mx_p',var=reshape(vari(ni1:ni2,nj1:nj2,nk1:nk2),(/nnode/)))
        call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                              var  = m_p (ni1:ni2+1,nj1:nj2+1,nk1:nk2+1)%y,             &
                              vari = vari(ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
        if (tec) err = tec_var(n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
        if (vtk) err = VTK_VAR_XML(NC_NN=nnode,varname='my_p',var=reshape(vari(ni1:ni2,nj1:nj2,nk1:nk2),(/nnode/)))
        call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                              var  = m_p (ni1:ni2+1,nj1:nj2+1,nk1:nk2+1)%z,             &
                              vari = vari(ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
        if (tec) err = tec_var(n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
        if (vtk) err = VTK_VAR_XML(NC_NN=nnode,varname='mz_p',var=reshape(vari(ni1:ni2,nj1:nj2,nk1:nk2),(/nnode/)))
        call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                              var  = m_v (ni1:ni2+1,nj1:nj2+1,nk1:nk2+1)%x,             &
                              vari = vari(ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
        if (tec) err = tec_var(n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
        if (vtk) err = VTK_VAR_XML(NC_NN=nnode,varname='mx_v',var=reshape(vari(ni1:ni2,nj1:nj2,nk1:nk2),(/nnode/)))
        call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                              var  = m_v (ni1:ni2+1,nj1:nj2+1,nk1:nk2+1)%y,             &
                              vari = vari(ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
        if (tec) err = tec_var(n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
        if (vtk) err = VTK_VAR_XML(NC_NN=nnode,varname='my_v',var=reshape(vari(ni1:ni2,nj1:nj2,nk1:nk2),(/nnode/)))
        call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                              var  = m_v (ni1:ni2+1,nj1:nj2+1,nk1:nk2+1)%z,             &
                              vari = vari(ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
        if (tec) err = tec_var(n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
        if (vtk) err = VTK_VAR_XML(NC_NN=nnode,varname='mz_v',var=reshape(vari(ni1:ni2,nj1:nj2,nk1:nk2),(/nnode/)))
      endif
      if (yp) then
        call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                              var  = yplus(ni1:ni2+1,nj1:nj2+1,nk1:nk2+1),              &
                              vari = vari( ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
        if (tec) err = tec_var(n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
        if (vtk) err = VTK_VAR_XML(NC_NN=nnode,varname='yplus',var=reshape(vari(ni1:ni2,nj1:nj2,nk1:nk2),(/nnode/)))
      endif
    endif
  endif
  if (metrics) then
    if (cell) then
      select case(face)
      case(1,2)
        if (tec) err = tec_var(n=ncell,var=NFi(ni1,cj1:cj2,ck1:ck2)%x,d=1)
        if (tec) err = tec_var(n=ncell,var=NFi(ni1,cj1:cj2,ck1:ck2)%y,d=1)
        if (tec) err = tec_var(n=ncell,var=NFi(ni1,cj1:cj2,ck1:ck2)%z,d=1)
        if (tec) err = tec_var(n=ncell,var= Si(ni1,cj1:cj2,ck1:ck2)  ,d=1)
        if (vtk) err = VTK_VAR_XML(NC_NN=ncell,varname='NFx',var=reshape(NFi(ni1,cj1:cj2,ck1:ck2)%x,(/ncell/)))
        if (vtk) err = VTK_VAR_XML(NC_NN=ncell,varname='NFy',var=reshape(NFi(ni1,cj1:cj2,ck1:ck2)%y,(/ncell/)))
        if (vtk) err = VTK_VAR_XML(NC_NN=ncell,varname='NFz',var=reshape(NFi(ni1,cj1:cj2,ck1:ck2)%z,(/ncell/)))
        if (vtk) err = VTK_VAR_XML(NC_NN=ncell,varname='S',  var=reshape( Si(ni1,cj1:cj2,ck1:ck2)  ,(/ncell/)))
      case(3,4)
        if (tec) err = tec_var(n=ncell,var=NFj(ci1:ci2,nj1,ck1:ck2)%x,d=1)
        if (tec) err = tec_var(n=ncell,var=NFj(ci1:ci2,nj1,ck1:ck2)%y,d=1)
        if (tec) err = tec_var(n=ncell,var=NFj(ci1:ci2,nj1,ck1:ck2)%z,d=1)
        if (tec) err = tec_var(n=ncell,var= Sj(ci1:ci2,nj1,ck1:ck2)  ,d=1)
        if (vtk) err = VTK_VAR_XML(NC_NN=ncell,varname='NFx',var=reshape(NFj(ci1:ci2,nj1,ck1:ck2)%x,(/ncell/)))
        if (vtk) err = VTK_VAR_XML(NC_NN=ncell,varname='NFy',var=reshape(NFj(ci1:ci2,nj1,ck1:ck2)%y,(/ncell/)))
        if (vtk) err = VTK_VAR_XML(NC_NN=ncell,varname='NFz',var=reshape(NFj(ci1:ci2,nj1,ck1:ck2)%z,(/ncell/)))
        if (vtk) err = VTK_VAR_XML(NC_NN=ncell,varname='S',  var=reshape( Sj(ci1:ci2,nj1,ck1:ck2)  ,(/ncell/)))
      case(5,6)
        if (tec) err = tec_var(n=ncell,var=NFk(ci1:ci2,cj1:cj2,nk1)%x,d=1)
        if (tec) err = tec_var(n=ncell,var=NFk(ci1:ci2,cj1:cj2,nk1)%y,d=1)
        if (tec) err = tec_var(n=ncell,var=NFk(ci1:ci2,cj1:cj2,nk1)%z,d=1)
        if (tec) err = tec_var(n=ncell,var= Sk(ci1:ci2,cj1:cj2,nk1)  ,d=1)
        if (vtk) err = VTK_VAR_XML(NC_NN=ncell,varname='NFx',var=reshape(NFk(ci1:ci2,cj1:cj2,nk1)%x,(/ncell/)))
        if (vtk) err = VTK_VAR_XML(NC_NN=ncell,varname='NFy',var=reshape(NFk(ci1:ci2,cj1:cj2,nk1)%y,(/ncell/)))
        if (vtk) err = VTK_VAR_XML(NC_NN=ncell,varname='NFz',var=reshape(NFk(ci1:ci2,cj1:cj2,nk1)%z,(/ncell/)))
        if (vtk) err = VTK_VAR_XML(NC_NN=ncell,varname='S',  var=reshape( Sk(ci1:ci2,cj1:cj2,nk1)  ,(/ncell/)))
      endselect
    else
      select case(face)
      case(1,2)
        call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                              var  =  NFi(ni1:ni2+1,nj1:nj2+1,nk1:nk2+1)%x,             &
                              vari = vari(ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
        if (tec) err = tec_var(n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
        if (vtk) err = VTK_VAR_XML(NC_NN=nnode,varname='NFx',var=reshape(vari(ni1:ni2,nj1:nj2,nk1:nk2),(/nnode/)))
        call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                              var  =  NFi(ni1:ni2+1,nj1:nj2+1,nk1:nk2+1)%y,             &
                              vari = vari(ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
        if (tec) err = tec_var(n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
        if (vtk) err = VTK_VAR_XML(NC_NN=nnode,varname='NFy',var=reshape(vari(ni1:ni2,nj1:nj2,nk1:nk2),(/nnode/)))
        call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                              var  =  NFi(ni1:ni2+1,nj1:nj2+1,nk1:nk2+1)%z,             &
                              vari = vari(ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
        if (tec) err = tec_var(n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
        if (vtk) err = VTK_VAR_XML(NC_NN=nnode,varname='NFz',var=reshape(vari(ni1:ni2,nj1:nj2,nk1:nk2),(/nnode/)))
        call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                              var  =   Si(ni1:ni2+1,nj1:nj2+1,nk1:nk2+1),               &
                              vari = vari(ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
        if (tec) err = tec_var(n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
        if (vtk) err = VTK_VAR_XML(NC_NN=nnode,varname='S',var=reshape(vari(ni1:ni2,nj1:nj2,nk1:nk2),(/nnode/)))
      case(3,4)
        call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                              var  =  NFj(ni1:ni2+1,nj1:nj2+1,nk1:nk2+1)%x,             &
                              vari = vari(ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
        if (tec) err = tec_var(n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
        if (vtk) err = VTK_VAR_XML(NC_NN=nnode,varname='NFx',var=reshape(vari(ni1:ni2,nj1:nj2,nk1:nk2),(/nnode/)))
        call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                              var  =  NFj(ni1:ni2+1,nj1:nj2+1,nk1:nk2+1)%y,             &
                              vari = vari(ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
        if (tec) err = tec_var(n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
        if (vtk) err = VTK_VAR_XML(NC_NN=nnode,varname='NFy',var=reshape(vari(ni1:ni2,nj1:nj2,nk1:nk2),(/nnode/)))
        call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                              var  =  NFj(ni1:ni2+1,nj1:nj2+1,nk1:nk2+1)%z,             &
                              vari = vari(ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
        if (tec) err = tec_var(n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
        if (vtk) err = VTK_VAR_XML(NC_NN=nnode,varname='NFz',var=reshape(vari(ni1:ni2,nj1:nj2,nk1:nk2),(/nnode/)))
        call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                              var  =   Sj(ni1:ni2+1,nj1:nj2+1,nk1:nk2+1),               &
                              vari = vari(ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
        if (tec) err = tec_var(n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
        if (vtk) err = VTK_VAR_XML(NC_NN=nnode,varname='S',var=reshape(vari(ni1:ni2,nj1:nj2,nk1:nk2),(/nnode/)))
      case(5,6)
        call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                              var  =  NFk(ni1:ni2+1,nj1:nj2+1,nk1:nk2+1)%x,             &
                              vari = vari(ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
        if (tec) err = tec_var(n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
        if (vtk) err = VTK_VAR_XML(NC_NN=nnode,varname='NFx',var=reshape(vari(ni1:ni2,nj1:nj2,nk1:nk2),(/nnode/)))
        call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                              var  =  NFk(ni1:ni2+1,nj1:nj2+1,nk1:nk2+1)%y,             &
                              vari = vari(ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
        if (tec) err = tec_var(n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
        if (vtk) err = VTK_VAR_XML(NC_NN=nnode,varname='NFy',var=reshape(vari(ni1:ni2,nj1:nj2,nk1:nk2),(/nnode/)))
        call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                              var  =  NFk(ni1:ni2+1,nj1:nj2+1,nk1:nk2+1)%z,             &
                              vari = vari(ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
        if (tec) err = tec_var(n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
        if (vtk) err = VTK_VAR_XML(NC_NN=nnode,varname='NFz',var=reshape(vari(ni1:ni2,nj1:nj2,nk1:nk2),(/nnode/)))
        call varinterpolation(ni1=ni1,ni2=ni2,nj1=nj1,nj2=nj2,nk1=nk1,nk2=nk2,face=face,&
                              var  =   Sk(ni1:ni2+1,nj1:nj2+1,nk1:nk2+1),               &
                              vari = vari(ni1:ni2  ,nj1:nj2  ,nk1:nk2  ))
        if (tec) err = tec_var(n=nnode,var=vari(ni1:ni2,nj1:nj2,nk1:nk2),d=1)
        if (vtk) err = VTK_VAR_XML(NC_NN=nnode,varname='S',var=reshape(vari(ni1:ni2,nj1:nj2,nk1:nk2),(/nnode/)))
      endselect
    endif
  endif
  if (vtk) then
    if (cell) then
      err = VTK_DAT_XML(var_location = 'cell', var_block_action = 'close')
    else
      err = VTK_DAT_XML(var_location = 'node', var_block_action = 'close')
    endif
    err = VTK_GEO_XML()
    err = VTK_END_XML()
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction patch_save

  subroutine varinterpolation(ni1,ni2,nj1,nj2,nk1,nk2,face,var,vari)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Subroutine for interpolating celle centered variable into node centered one.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P), intent(IN)::  ni1,ni2,nj1,nj2,nk1,nk2             ! Nodes bounds.
  integer(I_P), intent(IN)::  face                                ! Actual face.
  real(R_P),    intent(IN)::  var (ni1:ni2+1,nj1:nj2+1,nk1:nk2+1) ! Cell centered variable.
  real(R_P),    intent(OUT):: vari(ni1:ni2  ,nj1:nj2  ,nk1:nk2  ) ! Node centered interpolated variable.
  integer(I_P)::              i,j,k                               ! Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  do k=nk1,nk2
    do j=nj1,nj2
      do i=ni1,ni2
        select case(face)
        case(1,2) ! ni1=ni2 => interpolating patch ni1+1=ci1
          vari(i,j,k) = var(i+1,j+1,k+1) &
                      + var(i+1,j,  k+1) &
                      + var(i+1,j+1,k  ) &
                      + var(i+1,j  ,k  )
        case(3,4) ! nj1=nj2 => interpolating patch nj1+1=cj1
          vari(i,j,k) = var(i+1,j+1,k+1) &
                      + var(i  ,j+1,k+1) &
                      + var(i+1,j+1,k  ) &
                      + var(i  ,j+1,k  )
        case(5,6) ! nk1=nk2 => interpolating patch nk1+1=ck1
          vari(i,j,k) = var(i+1,j+1,k+1) &
                      + var(i  ,j+1,k+1) &
                      + var(i+1,j  ,k+1) &
                      + var(i  ,j  ,k+1)
        endselect
        vari(i,j,k) = 0.25_R_P*vari(i,j,k)
      enddo
    enddo
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine varinterpolation

  function tec_var(n,var,d) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Interface function for saving variables into Tecplot file.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P), intent(IN):: n        ! Number of var elements.
  real(R_P),    intent(IN):: var(1:n) ! Variable to be saved.
  integer(I_P), intent(IN):: d        ! Double precision output (1 yes, 0 no).
  integer(I_P)::             err      ! Error traping flag: 0 no errors, >0 error occours.
  integer(I_P)::             e        ! Element counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (binary) then
    err = tecdat112(n,var,d)
  else
    write(unit_out,FR_P,iostat=err)(var(e),e=1,n)
  endif
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction tec_var
endprogram XnPatches
