!> @brief This module contains the block allocatable variables
module Block_Variables
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision          ! Integers and reals precision definition.
USE Data_Type_PostProcess ! Definition of Type_PostProcess.
USE Data_Type_Vector      ! Definition of Type_Vector.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
public
save
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
type(Type_Vector), allocatable:: node(:,:,:)      ! Nodes coordinates.
type(Type_Vector), allocatable:: NFi(:,:,:)       ! Face i normal versor.
type(Type_Vector), allocatable:: NFj(:,:,:)       ! Face j normal versor.
type(Type_Vector), allocatable:: NFk(:,:,:)       ! Face k normal versor.
type(Type_Vector), allocatable:: NFiS(:,:,:)      ! Face i normal versor with surface area module.
type(Type_Vector), allocatable:: NFjS(:,:,:)      ! Face j normal versor with surface area module.
type(Type_Vector), allocatable:: NFkS(:,:,:)      ! Face k normal versor with surface area module.
real(R8P),         allocatable:: Si(:,:,:)        ! Face i area.
real(R8P),         allocatable:: Sj(:,:,:)        ! Face j area.
real(R8P),         allocatable:: Sk(:,:,:)        ! Face k area.
real(R8P),         allocatable:: volume(:,:,:)    ! Volumes of cells.
integer(I4P),      allocatable:: icc(:,:,:)       ! Cell centered icc values.
integer(I4P),      allocatable:: ricc(:,:,:)      ! Cell centered rcc values.
integer(I4P),      allocatable:: vicc(:,:,:)      ! Node centered rcc values.
type(Type_Vector), allocatable:: momentum(:,:,:)  ! Momentum.
real(R8P),         allocatable:: pressure(:,:,:)  ! Pressure.
real(R8P),         allocatable:: f(:,:,:)         ! Level set function.
real(R8P),         allocatable:: f0(:,:,:)        ! Level 0 (level set).
real(R8P),         allocatable:: visc(:,:,:)      ! Viscosity.
real(R8P),         allocatable:: vitl(:,:,:)      ! Turbulent viscosity.
real(R8P),         allocatable:: ken(:,:,:)       ! Turbulent kinetic energy.
real(R8P),         allocatable:: eps(:,:,:)       ! Turbulent kinetic energy dissipation.
real(R8P),         allocatable:: vord(:,:,:)      ! Variable to identify vortices (lambda 2).
real(R8P),         allocatable:: qfac(:,:,:)      ! Variable to identify vortices (q factor).
real(R8P),         allocatable:: heli(:,:,:)      ! Helicity.
type(Type_Vector), allocatable:: vorticity(:,:,:) ! Vorticity.
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @brief Procedure for allocating block variables.
  subroutine block_allocate(pp,gc,Ni,Nj,Nk)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_PostProcess), intent(IN):: pp       !< Post-processor data.
  integer(I4P),           intent(IN):: gc(1:6)  !< Number of ghost cells.
  integer(I4P),           intent(IN):: Ni,Nj,Nk !< Number of cells.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(node)) deallocate(node) ; allocate(node(0-gc(1):Ni+gc(2),0-gc(3):Nj+gc(4),0-gc(5):Nk+gc(6)))
  if (pp%fcc) then
    if (allocated(icc))  deallocate(icc)  ;                   allocate(icc (1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)))
    if (allocated(ricc)) deallocate(ricc) ;                   allocate(ricc(0      :Ni+1    ,0      :Nj+1    ,0      :Nk+1    ))
    if (allocated(vicc)) deallocate(vicc) ; if (.not.pp%cell) allocate(vicc(0      :Ni+1    ,0      :Nj+1    ,0      :Nk+1    ))
  endif
  if (pp%sol) then
    if (allocated(momentum)) deallocate(momentum) ; allocate(momentum(1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)))
    if (allocated(pressure)) deallocate(pressure) ; allocate(pressure(1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)))
    if (pp%level_set) then
      if (allocated(f))  deallocate(f)  ; allocate(f (1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)))
      if (allocated(f0)) deallocate(f0) ; allocate(f0(1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)))
    endif
    if (pp%zeroeq.or.pp%oneeq.or.pp%twoeq) then
      if (allocated(visc)) deallocate(visc) ; allocate(visc(1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)))
    endif
    if (pp%oneeq) then
      if (allocated(vitl)) deallocate(vitl) ; allocate(vitl(1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)))
    endif
    if (pp%twoeq) then
      if (allocated(ken)) deallocate(ken) ; allocate(ken(1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)))
      if (allocated(eps)) deallocate(eps) ; allocate(eps(1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)))
    endif
    if (pp%vordet) then ! computing vord variable
      if (allocated(NFi))       deallocate(NFi)       ; allocate(NFi      (0-gc(1):Ni+gc(2),0-gc(3):Nj+gc(4),0-gc(5):Nk+gc(6)))
      if (allocated(NFj))       deallocate(NFj)       ; allocate(NFj      (0-gc(1):Ni+gc(2),0-gc(3):Nj+gc(4),0-gc(5):Nk+gc(6)))
      if (allocated(NFk))       deallocate(NFk)       ; allocate(NFk      (0-gc(1):Ni+gc(2),0-gc(3):Nj+gc(4),0-gc(5):Nk+gc(6)))
      if (allocated(NFiS))      deallocate(NFiS)      ; allocate(NFiS     (0-gc(1):Ni+gc(2),0-gc(3):Nj+gc(4),0-gc(5):Nk+gc(6)))
      if (allocated(NFjS))      deallocate(NFjS)      ; allocate(NFjS     (0-gc(1):Ni+gc(2),0-gc(3):Nj+gc(4),0-gc(5):Nk+gc(6)))
      if (allocated(NFkS))      deallocate(NFkS)      ; allocate(NFkS     (0-gc(1):Ni+gc(2),0-gc(3):Nj+gc(4),0-gc(5):Nk+gc(6)))
      if (allocated(Si))        deallocate(Si)        ; allocate(Si       (0-gc(1):Ni+gc(2),0-gc(3):Nj+gc(4),0-gc(5):Nk+gc(6)))
      if (allocated(Sj))        deallocate(Sj)        ; allocate(Sj       (0-gc(1):Ni+gc(2),0-gc(3):Nj+gc(4),0-gc(5):Nk+gc(6)))
      if (allocated(Sk))        deallocate(Sk)        ; allocate(Sk       (0-gc(1):Ni+gc(2),0-gc(3):Nj+gc(4),0-gc(5):Nk+gc(6)))
      if (allocated(volume))    deallocate(volume)    ; allocate(volume   (1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)))
      if (allocated(vord))      deallocate(vord)      ; allocate(vord     (1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)))
      if (allocated(qfac))      deallocate(qfac)      ; allocate(qfac     (1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)))
      if (allocated(heli))      deallocate(heli)      ; allocate(heli     (1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)))
      if (allocated(vorticity)) deallocate(vorticity) ; allocate(vorticity(1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)))
    endif
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine block_allocate

  !> @brief Procedure for interpolating cell centered variable into node centered one.
  subroutine varinterpolation(Ni,Nj,Nk,var,vari)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P), intent(IN)::  Ni,Nj,Nk                   !< Block dimensions.
  real(R_P),    intent(IN)::  var (0:Ni+1,0:Nj+1,0:Nk+1) !< Cell centered variable.
  real(R_P),    intent(OUT):: vari(0:Ni  ,0:Nj  ,0:Nk  ) !< Node centered interpolated variable.
  integer(I_P)::              i,j,k                      !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(i,j,k)         &
  !$OMP SHARED(Ni,Nj,Nk,var,vari)
  !$OMP DO
  do k=0,Nk
    do j=0,Nj
      do i=0,Ni
        vari(i,j,k) = var(i+1,j+1,k+1)  + var(i,j+1,k+1) &
                    + var(i+1,j  ,k+1)  + var(i,j,  k+1) &
                    + var(i+1,j+1,k  )  + var(i,j+1,k  ) &
                    + var(i+1,j  ,k  )  + var(i,j  ,k  )
        vari(i,j,k) = 0.125_R_P*vari(i,j,k)
      enddo
    enddo
  enddo
  !$OMP END PARALLEL
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine varinterpolation
endmodule Block_Variables
