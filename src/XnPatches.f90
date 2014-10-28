program XnPatches
!-----------------------------------------------------------------------------------------------------------------------------------
! Code Name XnPatches
! Author    Stefano Zaghi
! Version   0.0.4
! Date      28-10-2014
!
! XnPatches loads Xnavis mesh and solution files and produces post-processed Tecplot file of selected patches.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision                                                                    ! Integers and reals precision definition.
USE Block_Variables                                                                 ! Block variables definition.
USE Data_Type_Command_Line_Interface                                                ! Definition of Type_Command_Line_Interface.
USE Data_Type_PostProcess                                                           ! Definition of Type_PostProcess.
USE Lib_TEC                                                                         ! Library for I/O Tecplot files.
!USE Lib_VTK_IO                                                                      ! Library for I/O VTK files.
USE, intrinsic:: ISO_FORTRAN_ENV, only: stdout => OUTPUT_UNIT, stderr => ERROR_UNIT ! Standard output/error logical units.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
type(Type_PostProcess)::     pp                !< Post-processor data.
integer(I4P)::               Nb = 0_I4P        !< Number of blocks.
integer(I4P),  allocatable:: gc(:,:)           !< Number of ghost cells.
integer(I4P),  allocatable:: Ni(:),Nj(:),Nk(:) !< Number of cells of each block.
real(R4P),     allocatable:: rcc(:)            !< rcc array.
integer(I4P)::               Nrcc = 0_I4P      !< rcc array dimension.
integer(I4P)::               Np = 0_I_P        !< Number of patches.
integer(I4P)::               error             !< Error traping flag: 0 no errors, >0 error occours.
integer(I4P)::               ni1,ni2           !< First and last node i indexes.
integer(I4P)::               nj1,nj2           !< First and last node j indexes.
integer(I4P)::               nk1,nk2           !< First and last node k indexes.
integer(I4P)::               ci1,ci2           !< First and last cell i indexes.
integer(I4P)::               cj1,cj2           !< First and last cell j indexes.
integer(I4P)::               ck1,ck2           !< First and last cell k indexes.
logical::                    patch_found       !< Inquiring flag for patches searching.
integer(I4P)::               i,j,k,b,v         !< Counters.
real(R8P)::                  t                 !< Solution time.
! misc variables
!character(100)::          vtkbdir             ! VTK base name of output directory.
!character(100)::          vtkbfile            ! VTK base name of output files.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
! parsing command line arguments
call parse_command_arguments(pp)

!if (vtk) then
!  ! making VTS_files directory
!  err = make_dir(directory=vtkbdir)
!endif

! loading mesh dimensions
read(pp%unit_grd)Nb
if (.not.pp%global_blk_num) then
  call pp%set_blockmap(Nb=Nb)
endif
read(pp%unit_icc)b
if (b/=Nb) then
  write(stderr,'(A)')'+--> Inconsistent ICC file with GRD one:'
  write(stderr,'(A)')'|--> Nb GRD file: '//trim(str(.true.,Nb))
  write(stderr,'(A)')'|--> Nb ICC file: '//trim(str(.true.,b))
  stop
endif
if (pp%sol) then
  read(pp%unit_sol)t
endif
write(stdout,'(A)')'+--> Number of block of input files: '//trim(str(.true.,Nb))
if (allocated(gc)) deallocate(gc) ; allocate(gc(1:6,1:Nb)) ; gc = 2
if (allocated(Ni)) deallocate(Ni) ; allocate(Ni(    1:Nb))
if (allocated(Nj)) deallocate(Nj) ; allocate(Nj(    1:Nb))
if (allocated(Nk)) deallocate(Nk) ; allocate(Nk(    1:Nb))
if (pp%verbose) write(stdout,'(A)')'+--> Blocks dimensions'
do b=1,Nb
  read(pp%unit_grd,err=1)Ni(b),Nj(b),Nk(b),gc(1,b),gc(2,b),gc(3,b),gc(4,b),gc(5,b),gc(6,b)
  1 continue
  if (pp%verbose) write(stdout,'(A)')'|-->   b,Ni,Nk,Nk,gc(6): '//trim(str(.true.,b))//', '//       &
                                                                  trim(str(.true.,Ni(b)))//', '//   &
                                                                  trim(str(.true.,Nj(b)))//', '//   &
                                                                  trim(str(.true.,Nk(b)))//', '//   &
                                                                  trim(str(.true.,gc(1,b)))//', '// &
                                                                  trim(str(.true.,gc(2,b)))//', '// &
                                                                  trim(str(.true.,gc(3,b)))//', '// &
                                                                  trim(str(.true.,gc(4,b)))//', '// &
                                                                  trim(str(.true.,gc(5,b)))//', '// &
                                                                  trim(str(.true.,gc(6,b)))
  read(pp%unit_icc,err=2)i,j,k,v,v,v,v,v,v
  2 continue
  if (i/=Ni(b).OR.j/=Nj(b).OR.k/=Nk(b)) then
    write(stderr,'(A)')'+--> Inconsistent ICC file with GRD one:'
    write(stderr,'(A)')'|--> Block: '//trim(str(.true.,b))
    write(stderr,'(A)')'|--> Ni,Nj,Nk GRD file: '//trim(str(.true.,Ni(b)))//','//&
                                                   trim(str(.true.,Nj(b)))//','//&
                                                   trim(str(.true.,Nk(b)))
    write(stderr,'(A)')'|--> Ni,Nj,Nk ICC file: '//trim(str(.true.,i))//','//trim(str(.true.,j))//','//trim(str(.true.,k))
    stop
  endif
  if (pp%sol) then
    read(pp%unit_sol)i,j,k
    if (i/=Ni(b).OR.j/=Nj(b).OR.k/=Nk(b)) then
      write(stderr,'(A)')'+--> Inconsistent SOL file with GRD one:'
      write(stderr,'(A)')'|--> Block: '//trim(str(.true.,b))
      write(stderr,'(A)')'|--> Ni,Nj,Nk GRD file: '//trim(str(.true.,Ni(b)))//','//&
                                                     trim(str(.true.,Nj(b)))//','//&
                                                     trim(str(.true.,Nk(b)))
      write(stderr,'(A)')'|--> Ni,Nj,Nk SOL file: '//trim(str(.true.,i))//','//trim(str(.true.,j))//','//trim(str(.true.,k))
      stop
    endif
  endif
enddo
do b=1,Nb
  read(pp%unit_icc)(((v,i=1-gc(1,b),Ni(b)+gc(2,b)),j=1-gc(3,b),Nj(b)+gc(4,b)),k=1-gc(5,b),Nk(b)+gc(6,b))
enddo
read(pp%unit_icc)Nrcc
if (pp%verbose) write(stdout,'(A)')'+--> RCC dimension '//trim(str(.true.,Nrcc))
if (allocated(rcc)) deallocate(rcc) ; allocate(rcc(1:Nrcc))
read(pp%unit_icc)(rcc(v),v=1,Nrcc)
! rewind icc file at icc pointer record
rewind(pp%unit_icc)
read(pp%unit_icc)b
do b=1,Nb
  read(pp%unit_icc,err=3)i,j,k,v,v,v,v,v,v
  3 continue
enddo

Np = 0 ! initialize patch number
do b=1,Nb ! post processing patches
  ! allocating mesh variables
  call block_allocate(pp=pp,gc=gc(:,b),Ni=Ni(b),Nj=Nj(b),Nk=Nk(b)) ! allocating block variables
  ! loading block nodes coordinates
  read(pp%unit_grd)(((node(i,j,k)%x,i=0-gc(1,b),Ni(b)+gc(2,b)),j=0-gc(3,b),Nj(b)+gc(4,b)),k=0-gc(5,b),Nk(b)+gc(6,b))
  read(pp%unit_grd)(((node(i,j,k)%y,i=0-gc(1,b),Ni(b)+gc(2,b)),j=0-gc(3,b),Nj(b)+gc(4,b)),k=0-gc(5,b),Nk(b)+gc(6,b))
  read(pp%unit_grd)(((node(i,j,k)%z,i=0-gc(1,b),Ni(b)+gc(2,b)),j=0-gc(3,b),Nj(b)+gc(4,b)),k=0-gc(5,b),Nk(b)+gc(6,b))
  read(pp%unit_icc)(((icc( i,j,k),  i=1-gc(1,b),Ni(b)+gc(2,b)),j=1-gc(3,b),Nj(b)+gc(4,b)),k=1-gc(5,b),Nk(b)+gc(6,b))
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
  if (.not.pp%cell) then ! ricc must be interpolated at nodes
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
  if (pp%metrics.or.pp%forces.or.pp%yp) then
    ! computing the metrics of cells
    call compute_metrics(gc=gc(:,b),Ni=Ni(b),Nj=Nj(b),Nk=Nk(b))
    ! correcting the metrics of boundary conditions cells
    call bc_metrics_correction(Ni=Ni(b),Nj=Nj(b),Nk=Nk(b),rcc=rcc)
  endif
  if (pp%sol) then
    read(pp%unit_sol)(((momentum(i,j,k)%x,i=0,Ni(b)+1),j=0,Nj(b)+1),k=0,Nk(b)+1)
    read(pp%unit_sol)(((momentum(i,j,k)%y,i=0,Ni(b)+1),j=0,Nj(b)+1),k=0,Nk(b)+1)
    read(pp%unit_sol)(((momentum(i,j,k)%z,i=0,Ni(b)+1),j=0,Nj(b)+1),k=0,Nk(b)+1)
    read(pp%unit_sol)(((pressure(i,j,k)  ,i=0,Ni(b)+1),j=0,Nj(b)+1),k=0,Nk(b)+1)
    if (pp%level_set) then
      read(pp%unit_sol)(((f (i,j,k),i=0,Ni(b)+1),j=0,Nj(b)+1),k=0,Nk(b)+1)
      read(pp%unit_sol)(((f0(i,j,k),i=0,Ni(b)+1),j=0,Nj(b)+1),k=0,Nk(b)+1)
    endif
    if (pp%zeroeq.or.pp%oneeq.or.pp%twoeq) then
      read(pp%unit_sol)(((visc(i,j,k),i=0,Ni(b)+1),j=0,Nj(b)+1),k=0,Nk(b)+1)
    endif
    if (pp%oneeq) then
      read(pp%unit_sol)(((vitl(i,j,k),i=0,Ni(b)+1),j=0,Nj(b)+1),k=0,Nk(b)+1)
    endif
    if (pp%twoeq) then
      read(pp%unit_sol)(((ken(i,j,k),i=0,Ni(b)+1),j=0,Nj(b)+1),k=0,Nk(b)+1)
      read(pp%unit_sol)(((eps(i,j,k),i=0,Ni(b)+1),j=0,Nj(b)+1),k=0,Nk(b)+1)
    endif
  endif
  ! checking the presence of the patch
  patch_found = (any(ticc(0,:,:)==pp%patch)).OR.(any(ticc(Ni(b)+1,:      ,:      )==pp%patch)).OR. &
                (any(ticc(:,0,:)==pp%patch)).OR.(any(ticc(:      ,Nj(b)+1,:      )==pp%patch)).OR. &
                (any(ticc(:,:,0)==pp%patch)).OR.(any(ticc(:      ,:      ,Nk(b)+1)==pp%patch))
  patch_found = (patch_found.OR.((any(ticc(0,:,:)==pp%patch+10)).OR.(any(ticc(Ni(b)+1,:      ,:      )==pp%patch+10)).OR. &
                                 (any(ticc(:,0,:)==pp%patch+10)).OR.(any(ticc(:      ,Nj(b)+1,:      )==pp%patch+10)).OR. &
                                 (any(ticc(:,:,0)==pp%patch+10)).OR.(any(ticc(:      ,:      ,Nk(b)+1)==pp%patch+10))))
  if (patch_found) then
    if (any(ticc(0,:,:)==pp%patch).OR.any(ticc(0,:,:)==pp%patch+10)) then
      Np = Np + 1

      ci1 = 1+pp%s_offset ; ci2 = 1+pp%s_offset
      cj1 = 1             ; cj2 = Nj(b)
      ck1 = 1             ; ck2 = Nk(b)

      ni1 = 0+pp%s_offset ; ni2 = 0+pp%s_offset
      nj1 = cj1-1         ; nj2 = cj2
      nk1 = ck1-1         ; nk2 = ck2

      if (pp%forces) call compute_forces(pp=pp,b=b,face=1,     &
                                         ni1 = ni1, ni2 = ni2, &
                                         nj1 = nj1, nj2 = nj2, &
                                         nk1 = nk1, nk2 = nk2, &
                                         ci1 = ci1, ci2 = ci2, &
                                         cj1 = cj1, cj2 = cj2, &
                                         ck1 = ck1, ck2 = ck2)
      if (pp%yp) call compute_yp(pp=pp,face=1,         &
                                 ni1 = ni1, ni2 = ni2, &
                                 nj1 = nj1, nj2 = nj2, &
                                 nk1 = nk1, nk2 = nk2, &
                                 ci1 = ci1, ci2 = ci2, &
                                 cj1 = cj1, cj2 = cj2, &
                                 ck1 = ck1, ck2 = ck2)
      if (pp%tec) error = file_tec%save_patch(Np=Np,pp=pp,b=b,face=1,&
                                              ni1 = ni1, ni2 = ni2,  &
                                              nj1 = nj1, nj2 = nj2,  &
                                              nk1 = nk1, nk2 = nk2,  &
                                              ci1 = ci1, ci2 = ci2,  &
                                              cj1 = cj1, cj2 = cj2,  &
                                              ck1 = ck1, ck2 = ck2)

    endif
    if (any(ticc(Ni(b)+1,:,:)==pp%patch).OR.any(ticc(Ni(b)+1,:,:)==pp%patch+10)) then
      Np = Np + 1

      ci1 = Ni(b)-pp%s_offset ; ci2 = Ni(b)-pp%s_offset
      cj1 = 1                 ; cj2 = Nj(b)
      ck1 = 1                 ; ck2 = Nk(b)

      ni1 = Ni(b)-pp%s_offset ; ni2 = Ni(b)-pp%s_offset
      nj1 = cj1-1             ; nj2 = cj2
      nk1 = ck1-1             ; nk2 = ck2

      if (pp%forces) call compute_forces(pp=pp,b=b,face=2,     &
                                         ni1 = ni1, ni2 = ni2, &
                                         nj1 = nj1, nj2 = nj2, &
                                         nk1 = nk1, nk2 = nk2, &
                                         ci1 = ci1, ci2 = ci2, &
                                         cj1 = cj1, cj2 = cj2, &
                                         ck1 = ck1, ck2 = ck2)
      if (pp%yp) call compute_yp(pp=pp,face=2,         &
                                 ni1 = ni1, ni2 = ni2, &
                                 nj1 = nj1, nj2 = nj2, &
                                 nk1 = nk1, nk2 = nk2, &
                                 ci1 = ci1, ci2 = ci2, &
                                 cj1 = cj1, cj2 = cj2, &
                                 ck1 = ck1, ck2 = ck2)
      if (pp%tec) error = file_tec%save_patch(Np=Np,pp=pp,b=b,face=2,&
                                              ni1 = ni1, ni2 = ni2,  &
                                              nj1 = nj1, nj2 = nj2,  &
                                              nk1 = nk1, nk2 = nk2,  &
                                              ci1 = ci1, ci2 = ci2,  &
                                              cj1 = cj1, cj2 = cj2,  &
                                              ck1 = ck1, ck2 = ck2)
    endif
    if (any(ticc(:,0,:)==pp%patch).OR.any(ticc(:,0,:)==pp%patch+10)) then
      Np = Np + 1

      cj1 = 1+pp%s_offset ; cj2 = 1+pp%s_offset
      ci1 = 1             ; ci2 = Ni(b)
      ck1 = 1             ; ck2 = Nk(b)

      nj1 = 0+pp%s_offset ; nj2 = 0+pp%s_offset
      ni1 = ci1-1         ; ni2 = ci2
      nk1 = ck1-1         ; nk2 = ck2

      if (pp%forces) call compute_forces(pp=pp,b=b,face=3,     &
                                         ni1 = ni1, ni2 = ni2, &
                                         nj1 = nj1, nj2 = nj2, &
                                         nk1 = nk1, nk2 = nk2, &
                                         ci1 = ci1, ci2 = ci2, &
                                         cj1 = cj1, cj2 = cj2, &
                                         ck1 = ck1, ck2 = ck2)
      if (pp%yp) call compute_yp(pp=pp,face=3,         &
                                 ni1 = ni1, ni2 = ni2, &
                                 nj1 = nj1, nj2 = nj2, &
                                 nk1 = nk1, nk2 = nk2, &
                                 ci1 = ci1, ci2 = ci2, &
                                 cj1 = cj1, cj2 = cj2, &
                                 ck1 = ck1, ck2 = ck2)
      if (pp%tec) error = file_tec%save_patch(Np=Np,pp=pp,b=b,face=3,&
                                              ni1 = ni1, ni2 = ni2,  &
                                              nj1 = nj1, nj2 = nj2,  &
                                              nk1 = nk1, nk2 = nk2,  &
                                              ci1 = ci1, ci2 = ci2,  &
                                              cj1 = cj1, cj2 = cj2,  &
                                              ck1 = ck1, ck2 = ck2)
    endif
    if (any(ticc(:,Nj(b)+1,:)==pp%patch).OR.any(ticc(:,Nj(b)+1,:)==pp%patch+10)) then
      Np = Np + 1

      cj1 = Nj(b)-pp%s_offset ; cj2 = Nj(b)-pp%s_offset
      ci1 = 1                 ; ci2 = Ni(b)
      ck1 = 1                 ; ck2 = Nk(b)

      nj1 = Nj(b)-pp%s_offset ; nj2 = Nj(b)-pp%s_offset
      ni1 = ci1-1             ; ni2 = ci2
      nk1 = ck1-1             ; nk2 = ck2

      if (pp%forces) call compute_forces(pp=pp,b=b,face=4,     &
                                         ni1 = ni1, ni2 = ni2, &
                                         nj1 = nj1, nj2 = nj2, &
                                         nk1 = nk1, nk2 = nk2, &
                                         ci1 = ci1, ci2 = ci2, &
                                         cj1 = cj1, cj2 = cj2, &
                                         ck1 = ck1, ck2 = ck2)
      if (pp%yp) call compute_yp(pp=pp,face=4,         &
                                 ni1 = ni1, ni2 = ni2, &
                                 nj1 = nj1, nj2 = nj2, &
                                 nk1 = nk1, nk2 = nk2, &
                                 ci1 = ci1, ci2 = ci2, &
                                 cj1 = cj1, cj2 = cj2, &
                                 ck1 = ck1, ck2 = ck2)
      if (pp%tec) error = file_tec%save_patch(Np=Np,pp=pp,b=b,face=4,&
                                              ni1 = ni1, ni2 = ni2,  &
                                              nj1 = nj1, nj2 = nj2,  &
                                              nk1 = nk1, nk2 = nk2,  &
                                              ci1 = ci1, ci2 = ci2,  &
                                              cj1 = cj1, cj2 = cj2,  &
                                              ck1 = ck1, ck2 = ck2)
    endif
    if (any(ticc(:,:,0)==pp%patch).OR.any(ticc(:,:,0)==pp%patch+10)) then
      Np = Np + 1

      ck1 = 1+pp%s_offset ; ck2 = 1+pp%s_offset
      ci1 = 1             ; ci2 = Ni(b)
      cj1 = 1             ; cj2 = Nj(b)

      nk1 = 0+pp%s_offset ; nk2 = 0+pp%s_offset
      ni1 = ci1-1         ; ni2 = ci2
      nj1 = cj1-1         ; nj2 = cj2

      if (pp%forces) call compute_forces(pp=pp,b=b,face=5,     &
                                         ni1 = ni1, ni2 = ni2, &
                                         nj1 = nj1, nj2 = nj2, &
                                         nk1 = nk1, nk2 = nk2, &
                                         ci1 = ci1, ci2 = ci2, &
                                         cj1 = cj1, cj2 = cj2, &
                                         ck1 = ck1, ck2 = ck2)
      if (pp%yp) call compute_yp(pp=pp,face=5,         &
                                 ni1 = ni1, ni2 = ni2, &
                                 nj1 = nj1, nj2 = nj2, &
                                 nk1 = nk1, nk2 = nk2, &
                                 ci1 = ci1, ci2 = ci2, &
                                 cj1 = cj1, cj2 = cj2, &
                                 ck1 = ck1, ck2 = ck2)
      if (pp%tec) error = file_tec%save_patch(Np=Np,pp=pp,b=b,face=5,&
                                              ni1 = ni1, ni2 = ni2,  &
                                              nj1 = nj1, nj2 = nj2,  &
                                              nk1 = nk1, nk2 = nk2,  &
                                              ci1 = ci1, ci2 = ci2,  &
                                              cj1 = cj1, cj2 = cj2,  &
                                              ck1 = ck1, ck2 = ck2)
    endif
    if (any(ticc(:,:,Nk(b)+1)==pp%patch).OR.any(ticc(:,:,Nk(b)+1)==pp%patch+10)) then
      Np = Np + 1

      ck1 = Nk(b)-pp%s_offset ; ck2 = Nk(b)-pp%s_offset
      ci1 = 1                 ; ci2 = Ni(b)
      cj1 = 1                 ; cj2 = Nj(b)

      nk1 = Nk(b)-pp%s_offset ; nk2 = Nk(b)-pp%s_offset
      ni1 = ci1-1             ; ni2 = ci2
      nj1 = cj1-1             ; nj2 = cj2

      if (pp%forces) call compute_forces(pp=pp,b=b,face=6,     &
                                         ni1 = ni1, ni2 = ni2, &
                                         nj1 = nj1, nj2 = nj2, &
                                         nk1 = nk1, nk2 = nk2, &
                                         ci1 = ci1, ci2 = ci2, &
                                         cj1 = cj1, cj2 = cj2, &
                                         ck1 = ck1, ck2 = ck2)
      if (pp%yp) call compute_yp(pp=pp,face=6,         &
                                 ni1 = ni1, ni2 = ni2, &
                                 nj1 = nj1, nj2 = nj2, &
                                 nk1 = nk1, nk2 = nk2, &
                                 ci1 = ci1, ci2 = ci2, &
                                 cj1 = cj1, cj2 = cj2, &
                                 ck1 = ck1, ck2 = ck2)
      if (pp%tec) error = file_tec%save_patch(Np=Np,pp=pp,b=b,face=6,&
                                              ni1 = ni1, ni2 = ni2,  &
                                              nj1 = nj1, nj2 = nj2,  &
                                              nk1 = nk1, nk2 = nk2,  &
                                              ci1 = ci1, ci2 = ci2,  &
                                              cj1 = cj1, cj2 = cj2,  &
                                              ck1 = ck1, ck2 = ck2)
    endif
  else
    cycle
  endif
enddo

if (pp%forces) error = pp%save_forces()
call pp%finalize()
if (pp%tec) call file_tec%finalize(pp=pp)
if (allocated(rcc)) deallocate(rcc)

!if (vtk) then
!  ! saving the multiblock VTM wrapper
!  err = VTM_INI_XML(filename=adjustl(trim(File_out))//"."//trim(strz(2,patch))//'.vtm')
!  err = VTM_BLK_XML(block_action='open')
!  err = VTM_WRF_XML(flist=(/(adjustl(trim(OS%basename(File_out)))//'-vts_files'//OS%sep//     &
!                             adjustl(trim(OS%basename(File_out)))//"."//trim(strz(2,patch))// &
!                            '_n_'//trim(strz(5,b))//'.vts',                                   &
!                             b=1,Np)/))
!  err = VTM_BLK_XML(block_action='close')
!  err = VTM_END_XML()
!endif

if (Np==0) then
  write(stdout,'(A)')'+--> No patches have been found!'
endif
stop
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @brief Procedure for parsing command line arguments.
  subroutine parse_command_arguments(pp)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_PostProcess), intent(INOUT):: pp    !< Post-processor data.
  type(Type_Command_Line_Interface)::     cli   !< Command Line Interface (CLI).
  integer(I4P)::                          error !< Error trapping flag.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(stdout,'(A)')'+--> XnPatches, post-processor for Xnavis '
  ! initializing CLI
  call cli%init(progname='XnPatches',                                                            &
                version ='v0.0.2',                                                               &
                help    =' The XnPatches Command Line Interface (CLI) has the following options',&
                examples=["XnPatches -g cc.01.grd -i cc.01 -o wall",                             &
                          "XnPatches -g cc.01.grd -i cc.01 -s sol.00.01 -forces -Re 1d6 -o sol", &
                          "XnPatches -g cc.01.grd -i cc.01 -s sol.00.01 -o sol"])
  ! setting CLAs
  call cli%add(pref='|-->',switch='-g',       help='Grid file (.grd)',required=.true.,act='store',error=error)
  call cli%add(pref='|-->',switch='-i',       help='ICC file',required=.true.,act='store',error=error)
  call cli%add(pref='|-->',switch='-o',       help='Output file name; default is basename of grd file with the proper extension',&
                                              required=.false.,act='store',def='unset',error=error)
  call cli%add(pref='|-->',switch='-s',       help='Solution file name; if passed the solution variables are saved',&
                                              required=.false.,act='store',def='unset',error=error)
  call cli%add(pref='|-->',switch='-p',       help='Boundary condition processed',required=.false.,act='store',&
                                              def='1',error=error)
  call cli%add(pref='|-->',switch='-cell',    help='Variables other than nodes coord. are cell centered',&
                                              required=.false.,act='store_true',def='.false.',error=error)
  call cli%add(pref='|-->',switch='-ls',      help='Solution with level set',&
                                              required=.false.,act='store_true',def='.false.',error=error)
  call cli%add(pref='|-->',switch='-nt',      help='No turbulent model used',&
                                              required=.false.,act='store_true',def='.false.',error=error)
  call cli%add(pref='|-->',switch='-eq',      help='Equations turbulent model',required=.false.,act='store',&
                                              def='1',choices='0,1,2',error=error)
  call cli%add(pref='|-->',switch='-stream',  help='Streamlines offset: consider the cell "patch+#value" (along i,j,k direction)'//&
                                                   ' instead of cell "patch" for visualizing limiting streamlines',                &
                                              required=.false.,act='store',def='0',error=error)
  call cli%add(pref='|-->',switch='-Re',      help='Reynolds number',required=.false.,act='store',def='-1.E0',error=error)
  call cli%add(pref='|-->',switch='-Fr',      help='Froude number',required=.false.,act='store',def='-1.E0',error=error)
  call cli%add(pref='|-->',switch='-zfs',     help='Quote of free surface',required=.false.,act='store',def='1.E100',error=error)
  call cli%add(pref='|-->',switch='-forces',  help='Compute forces, Reynolds number is necessary',&
                                              required=.false.,act='store_true',def='.false.',error=error)
  call cli%add(pref='|-->',switch='-yplus',   help='Compute y+, Reynolds number and solution are necessary',&
                                              required=.false.,act='store_true',def='.false.',error=error)
  call cli%add(pref='|-->',switch='-metrics', help='Save metrics variables',&
                                              required=.false.,act='store_true',def='.false.',error=error)
  call cli%add(pref='|-->',switch='-ascii',   help='write ascii output files',&
                                              required=.false.,act='store_true',def='.false.',error=error)
  call cli%add(pref='|-->',switch='-tec',     help='write output Tecplot files',&
                                              required=.false.,act='store',def='yes',choices='yes,no',error=error)
  call cli%add(pref='|-->',switch='-vtk',     help='write output VTK files',&
                                              required=.false.,act='store',def='no',choices='yes,no',error=error)
  call cli%add(pref='|-->',switch='-proc',    help='process number for global block numeration if proc.input is found',&
                                              required=.false.,act='store',def='-1',error=error)
  call cli%add(pref='|-->',switch='-os',      help='type of Operating System',&
                                              required=.false.,act='store',def='UIX',choices='UIX,WIN',error=error)
  call cli%add(pref='|-->',switch='-vb',      help='Verbose output',required=.false.,act='store_true',def='.false.',error=error)
  ! parsing CLI
  write(stdout,'(A)')'+--> Parsing Command Line Arguments'
  call cli%parse(error=error,pref='|-->'); if (error/=0) stop
  ! setting post-processing options accordingly to CLI
  call pp%set_from_cli(cli=cli)
  write(stdout,'(A)')'+--> Inquiring mb.par'
  call pp%set_from_mbpar()
  write(stdout,'(A)')'+--> Inquiring proc.input'
  call pp%set_from_procinput()
  write(stdout,'(A)')'+--> Initializing IO files'
  call pp%input_files_init()
  if (pp%tec) call file_tec%init(pp=pp)
  !if (vtk) then
  !  vtkbdir  = adjustl(trim(OS%basedir(File_out)))//adjustl(trim(OS%basename(File_out)))//'-vts_files'//OS%sep
  !  vtkbfile = adjustl(trim(vtkbdir))//adjustl(trim(OS%basename(File_out)))
  !endif
              write(stdout,'(A)')'+--> GRD File '//trim(pp%File_grd)
              write(stdout,'(A)')'|--> ICC File '//trim(pp%File_icc)
  if (pp%sol) write(stdout,'(A)')'|--> SOL File '//trim(pp%File_sol)
              write(stdout,'(A)')'|--> OUT File '//trim(pp%File_out)
  ! checking that all the info for forces computing have been provided
  if (pp%forces) then
    if (.not.pp%sol) then
      write(stderr,'(A)')'+--> Trying to compute forces without provide solution file!'
      stop
    endif
    if (pp%Re<0._R_P) then
      write(stderr,'(A)')'+--> Trying to compute forces without provide Reynolds number!'
      write(stderr,'(A)')'|--> The Reynolds number must be provide either by command line argument or by mb.par'
      stop
    endif
    if (pp%level_set) then
      if (pp%Fr<0._R_P) then
        write(stderr,'(A)')'+--> Trying to compute forces without provide Froude number for a simulation with free surface!'
        write(stderr,'(A)')'|--> The Froude number must be provide either by command line argument or by mb.par'
        stop
      elseif (pp%zfs>0.7D100) then
        write(stderr,'(A)')'+--> Trying to compute forces without provide quote of free surface for a simulation with free surface!'
        write(stderr,'(A)')'|--> The quote of free surface must be provide either by command line argument or by mb.par'
        stop
      endif
    endif
  endif
  ! checking that all the info for yplus computing have been provided
  if (pp%yp) then
    if (.not.pp%sol) then
      write(stderr,'(A)')'+-->  Trying to compute y+ without provide solution file!'
      stop
    endif
    if (pp%Re<0._R_P) then
      write(stderr,'(A)')'+--> Trying to compute y+ without provide Reynolds number!'
      write(stderr,'(A)')'|--> The Reynolds number must be provide either by command line argument or by mb.par'
      stop
    endif
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine parse_command_arguments
endprogram XnPatches
