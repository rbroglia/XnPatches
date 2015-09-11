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
integer(I4P),      allocatable:: ticc(:,:,:)      ! Cell centered rcc values (only type of cc).
integer(I4P),      allocatable:: vicc(:,:,:)      ! Node centered rcc values.
type(Type_Vector), allocatable:: momentum(:,:,:)  ! Momentum.
real(R8P),         allocatable:: pressure(:,:,:)  ! Pressure.
real(R8P),         allocatable:: f(:,:,:)         ! Level set function.
real(R8P),         allocatable:: f0(:,:,:)        ! Level 0 (level set).
real(R8P),         allocatable:: visc(:,:,:)      ! Viscosity.
real(R8P),         allocatable:: vitl(:,:,:)      ! Turbulent viscosity.
real(R8P),         allocatable:: ken(:,:,:)       ! Turbulent kinetic energy.
real(R8P),         allocatable:: eps(:,:,:)       ! Turbulent kinetic energy dissipation.
real(R8P),         allocatable:: yplus(:,:,:)     ! Estimation of y+.
type(Type_Vector), allocatable:: f_p(:,:,:)       ! Forces vector, pressure part.
type(Type_Vector), allocatable:: f_v(:,:,:)       ! Forces vector, viscous  part.
type(Type_Vector), allocatable:: m_p(:,:,:)       ! Moments vector, pressure part.
type(Type_Vector), allocatable:: m_v(:,:,:)       ! Moments vector, viscous  part.
type(Type_Vector), allocatable:: tau(:,:,:)       ! Ambiguos vector.
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
  if (allocated(icc))  deallocate(icc)  ;                   allocate(icc (1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)))
  if (allocated(ricc)) deallocate(ricc) ;                   allocate(ricc(0      :Ni+1    ,0      :Nj+1    ,0      :Nk+1    ))
  if (allocated(ticc)) deallocate(ticc) ;                   allocate(ticc(0      :Ni+1    ,0      :Nj+1    ,0      :Nk+1    ))
  if (allocated(vicc)) deallocate(vicc) ; if (.not.pp%cell) allocate(vicc(0      :Ni+1    ,0      :Nj+1    ,0      :Nk+1    ))
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
    if (pp%metrics.or.pp%forces.or.pp%yp.or.pp%tau) then
      if (allocated(NFi))    deallocate(NFi)    ; allocate(NFi   (0-gc(1):Ni+gc(2),0-gc(3):Nj+gc(4),0-gc(5):Nk+gc(6)))
      if (allocated(NFj))    deallocate(NFj)    ; allocate(NFj   (0-gc(1):Ni+gc(2),0-gc(3):Nj+gc(4),0-gc(5):Nk+gc(6)))
      if (allocated(NFk))    deallocate(NFk)    ; allocate(NFk   (0-gc(1):Ni+gc(2),0-gc(3):Nj+gc(4),0-gc(5):Nk+gc(6)))
      if (allocated(NFiS))   deallocate(NFiS)   ; allocate(NFiS  (0-gc(1):Ni+gc(2),0-gc(3):Nj+gc(4),0-gc(5):Nk+gc(6)))
      if (allocated(NFjS))   deallocate(NFjS)   ; allocate(NFjS  (0-gc(1):Ni+gc(2),0-gc(3):Nj+gc(4),0-gc(5):Nk+gc(6)))
      if (allocated(NFkS))   deallocate(NFkS)   ; allocate(NFkS  (0-gc(1):Ni+gc(2),0-gc(3):Nj+gc(4),0-gc(5):Nk+gc(6)))
      if (allocated(Si))     deallocate(Si)     ; allocate(Si    (0-gc(1):Ni+gc(2),0-gc(3):Nj+gc(4),0-gc(5):Nk+gc(6)))
      if (allocated(Sj))     deallocate(Sj)     ; allocate(Sj    (0-gc(1):Ni+gc(2),0-gc(3):Nj+gc(4),0-gc(5):Nk+gc(6)))
      if (allocated(Sk))     deallocate(Sk)     ; allocate(Sk    (0-gc(1):Ni+gc(2),0-gc(3):Nj+gc(4),0-gc(5):Nk+gc(6)))
      if (allocated(volume)) deallocate(volume) ; allocate(volume(1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6)))
    endif
    if (pp%forces) then
      if (allocated(f_p)) deallocate(f_p) ; allocate(f_p(0-gc(1):Ni+gc(2),0-gc(3):Nj+gc(4),0-gc(5):Nk+gc(6))) ; f_p = 0._R_P
      if (allocated(f_v)) deallocate(f_v) ; allocate(f_v(0-gc(1):Ni+gc(2),0-gc(3):Nj+gc(4),0-gc(5):Nk+gc(6))) ; f_v = 0._R_P
      if (allocated(m_p)) deallocate(m_p) ; allocate(m_p(0-gc(1):Ni+gc(2),0-gc(3):Nj+gc(4),0-gc(5):Nk+gc(6))) ; m_p = 0._R_P
      if (allocated(m_v)) deallocate(m_v) ; allocate(m_v(0-gc(1):Ni+gc(2),0-gc(3):Nj+gc(4),0-gc(5):Nk+gc(6))) ; m_v = 0._R_P
    endif
    if (pp%yp) then
      if (allocated(yplus)) deallocate(yplus) ; allocate(yplus(0-gc(1):Ni+gc(2),0-gc(3):Nj+gc(4),0-gc(5):Nk+gc(6)))
    endif
    if (pp%tau) then
      if (allocated(tau)) deallocate(tau) ; allocate(tau(0-gc(1):Ni+gc(2),0-gc(3):Nj+gc(4),0-gc(5):Nk+gc(6)))
    endif
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine block_allocate

  !> @brief Procedure for computing block metrics.
  subroutine compute_metrics(gc,Ni,Nj,Nk)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(IN):: gc(1:6)           !< Number of ghost cells.
  integer(I4P), intent(IN):: Ni,Nj,Nk          !< Number of cells.
  type(Type_Vector)::        NFS,s1,s2,nd,db   !< Dummy vector variables.
  real(R8P)::                signi,signj,signk !< Dummy variables for checking the directions of normals.
  real(R8P)::                Vx,Vy,Vz          !< Dummy variables for computing volume.
  real(R8P)::                xp,yp,zp          !< Dummy variables for computing face coordinates.
  real(R8P)::                xm,ym,zm          !< Dummy variables for computing face coordinates.
  integer(I4P)::             i,j,k             !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! computing faces normals
  ! positioning at the middle of the block
  i = max(1,Ni/2)
  j = max(1,Nj/2)
  k = max(1,Nk/2)
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
  !$OMP SHARED(Ni,Nj,Nk,gc,signi,signj,signk,node,Si,Sj,Sk,NFi,NFj,NFk,NFiS,NFjS,NFkS,volume)
  !$OMP DO
  do k=1-gc(5),Nk+gc(6)
    do j=1-gc(3),Nj+gc(4)
      do i=0-gc(1),Ni+gc(2)
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
  do k=1-gc(5),Nk+gc(6)
    do j=0-gc(3),Nj+gc(4)
      do i=1-gc(1),Ni+gc(2)
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
  do k=0-gc(5),Nk+gc(6)
    do j=1-gc(3),Nj+gc(4)
      do i=1-gc(1),Ni+gc(2)
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
  do k=1-gc(5),Nk+gc(6)
    do j=1-gc(3),Nj+gc(4)
      do i=1-gc(1),Ni+gc(2)
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

  !> @brief Subroutine for correcting the metrics of natural (and negative volume) boundary conditions cells.
  subroutine bc_metrics_correction(Ni,Nj,Nk,rcc)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(IN):: Ni,Nj,Nk           !< Number of cells.
  real(R4P),    intent(IN):: rcc(1:)            !< rcc array.
  logical::                  bc_correct         !< Flag for inquiring if the bc metrics must be corrected.
  logical::                  bc_wall            !< Flag for inquiring if the bc is "wall-type": different corrections must be used.
  real(R8P)::                tm                 !< Tangential metrics parameter (-1 for wall-type bc).
  real(R8P)::                sn                 !< Normal     metrics coefficient correction.
  integer(I4P)::             i,j,k              !< counters.
  integer(I4P), parameter::  wall         = -1  !< Wall boundary condition.
  integer(I4P), parameter::  simmetry     = -2  !< Simmetry boundary condition.
  integer(I4P), parameter::  movingwall   = -10 !< Moving wall boundary condition.
  integer(I4P), parameter::  passivewall  = -11 !< Passive wall boundary condition.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  !$OMP PARALLEL DEFAULT(NONE)                  &
  !$OMP PRIVATE(i,j,k,bc_correct,bc_wall,tm,sn) &
  !$OMP SHARED(Ni,Nj,Nk,rcc,icc,NFi,NFj,NFk,NFiS,NFjS,NFkS,volume)
  ! left i
  !$OMP DO
  do k=1,Nk
    do j=1,Nj
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
  do k=1,Nk
    do j=1,Nj
      i = nint(rcc(icc(Ni+1,j,k)))
      bc_correct = ((i<0).OR.(volume(Ni+1,j,k)<(0.2_R_P*volume(Ni,j,k))))
      bc_wall    = ((i==wall).OR.(i==simmetry).OR.(i==movingwall).OR.(i==passivewall))
      tm = 1._R_P
      if (bc_wall) tm = -1._R_P
      if (bc_correct) then
         ! normal metrics
         sn = 2._R_P*(NFiS(Ni-1,j,k).dot.NFi(Ni,j,k))
         NFiS(  Ni+1,j,  k  ) = -NFiS(Ni-1,j,k)+sn*NFi(Ni,j,k)
         ! tangential metrics
         NFjS(  Ni+1,j  ,k  ) = tm*NFjS(Ni,j  ,k  )
         NFjS(  Ni+1,j-1,k  ) = tm*NFjS(Ni,j-1,k  )
         NFjS(  Ni+1,j  ,k-1) = tm*NFjS(Ni,j  ,k-1)
         NFjS(  Ni+1,j-1,k-1) = tm*NFjS(Ni,j-1,k-1)

         NFkS(  Ni+1,j  ,k  ) = tm*NFkS(Ni,j  ,k  )
         NFkS(  Ni+1,j-1,k  ) = tm*NFkS(Ni,j-1,k  )
         NFkS(  Ni+1,j  ,k-1) = tm*NFkS(Ni,j  ,k-1)
         NFkS(  Ni+1,j-1,k-1) = tm*NFkS(Ni,j-1,k-1)
         ! volume
         volume(Ni+1,j,  k  ) = volume(Ni,j,k)
      end if
    enddo
  enddo
  ! left j
  !$OMP DO
  do k=1,Nk
    do i=1,Ni
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
  do k=1,Nk
    do i=1,Ni
      j = nint(rcc(icc(i,Nj+1,k)))
      bc_correct = ((j<0).OR.(volume(i,Nj+1,k)<(0.2_R_P*volume(i,Nj,k))))
      bc_wall    = ((j==wall).OR.(j==simmetry).OR.(j==movingwall).OR.(j==passivewall))
      tm = 1._R_P
      if (bc_wall) tm = -1._R_P
      if (bc_correct) then
         ! normal metrics
         sn = 2._R_P*(NFjS(i,Nj-1,k).dot.NFj(i,Nj,k))
         NFjS(  i,Nj+1,  k  ) = -NFjS(i,Nj-1,k)+sn*NFj(i,Nj,k)
         ! tangential metrics
         NFiS(  i  ,Nj+1,k  ) = tm*NFiS(i  ,Nj,k  )
         NFiS(  i-1,Nj+1,k  ) = tm*NFiS(i-1,Nj,k  )
         NFiS(  i  ,Nj+1,k-1) = tm*NFiS(i  ,Nj,k-1)
         NFiS(  i-1,Nj+1,k-1) = tm*NFiS(i-1,Nj,k-1)

         NFkS(  i  ,Nj+1,k  ) = tm*NFkS(i  ,Nj,k  )
         NFkS(  i-1,Nj+1,k  ) = tm*NFkS(i-1,Nj,k  )
         NFkS(  i  ,Nj+1,k-1) = tm*NFkS(i  ,Nj,k-1)
         NFkS(  i-1,Nj+1,k-1) = tm*NFkS(i-1,Nj,k-1)
         ! volume
         volume(i,Nj+1,  k  ) = volume(i,Nj,k)
      end if
    enddo
  enddo
  ! left k
  !$OMP DO
  do j=1,Nj
    do i=1,Ni
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
  do j=1,Nj
    do i=1,Ni
      k = nint(rcc(icc(i,j,Nk+1)))
      bc_correct = ((k<0).OR.(volume(i,j,Nk+1)<(0.2_R_P*volume(i,j,Nk))))
      bc_wall    = ((k==wall).OR.(k==simmetry).OR.(k==movingwall).OR.(k==passivewall))
      tm = 1._R_P
      if (bc_wall) tm = -1._R_P
      if (bc_correct) then
         ! normal metrics
         sn = 2._R_P*(NFkS(i,j,Nk-1).dot.NFk(i,j,Nk))
         NFkS(  i,  j,  Nk+1) = -NFkS(i,j,Nk-1)+sn*NFk(i,j,Nk)
         ! tangential metrics
         NFiS(  i  ,j  ,Nk+1) = tm*NFiS(i  ,j  ,Nk)
         NFiS(  i-1,j  ,Nk+1) = tm*NFiS(i-1,j  ,Nk)
         NFiS(  i  ,j-1,Nk+1) = tm*NFiS(i  ,j-1,Nk)
         NFiS(  i-1,j-1,Nk+1) = tm*NFiS(i-1,j-1,Nk)

         NFjS(  i  ,j  ,Nk+1) = tm*NFjS(i  ,j  ,Nk)
         NFjS(  i-1,j  ,Nk+1) = tm*NFjS(i-1,j  ,Nk)
         NFjS(  i  ,j-1,Nk+1) = tm*NFjS(i  ,j-1,Nk)
         NFjS(  i-1,j-1,Nk+1) = tm*NFjS(i-1,j-1,Nk)
         ! volume
         volume(i,  j,  Nk+1) = volume(i,j,Nk)
      end if
    enddo
  enddo
  !$OMP END PARALLEL
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine bc_metrics_correction

  !> @brief Procedure for computing forces variables.
  subroutine compute_forces(pp,b,face,ni1,ni2,nj1,nj2,nk1,nk2,ci1,ci2,cj1,cj2,ck1,ck2)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_PostProcess), intent(INOUT):: pp                      !< Post-processor data.
  integer(I4P),           intent(IN)::    b                       !< Actual block.
  integer(I4P),           intent(IN)::    face                    !< Face where patch is defined: 1,2,3,4,5,6.
  integer(I4P),           intent(IN)::    ni1,ni2                 !< First and last node i indexes.
  integer(I4P),           intent(IN)::    nj1,nj2                 !< First and last node j indexes.
  integer(I4P),           intent(IN)::    nk1,nk2                 !< First and last node k indexes.
  integer(I4P),           intent(IN)::    ci1,ci2                 !< First and last cell i indexes.
  integer(I4P),           intent(IN)::    cj1,cj2                 !< First and last cell j indexes.
  integer(I4P),           intent(IN)::    ck1,ck2                 !< First and last cell k indexes.
  real(R8P)::                             fsumpx,fsumpy,fsumpz, & !< |
                                          fsumvx,fsumvy,fsumvz, & !< | Dummies vars for computing sum of forces in OpenMP parallel
                                          msumpx,msumpy,msumpz, & !< | blocks. OpenMP doesn't support reduction on derived type var.
                                          msumvx,msumvy,msumvz    !< |
  type(Type_Vector)::                     NdS                     !< Normal "tilde" for viscous part of forces (or distance for y+).
  type(Type_Vector)::                     fco                     !< Pacth center coordinates.
  real(R8P)::                             hyd                     !< Hydrostatic correction of pressure term.
  integer(I4P)::                          i,j,k                   !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
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
        if (ticc(ni1-1+face,j,k)/=pp%patch) cycle
        if ( icc(ci1,j,k)/=0) cycle
        if(pp%level_set) then
          if (f0(ni1-1+face,j,k)>0._R_P) cycle
        endif

        ! patch center coordinates
        fco = 0.25_R_P*(node(ni1,j,k)+node(ni1,j-1,k)+node(ni1,j,k-1)+node(ni1,j-1,k-1))

        ! pressure part of forces
        hyd = 0._R_P
        if(pp%level_set) hyd = (fco%z - pp%zfs)*pp%rFr2
        f_p(ci1,j,k) = -(1.5_R_P*pressure(ci1,j,k) - 0.5_R_P*pressure(ci1+3-2*face,j,k) - hyd)*NFiS(ni1,j,k)
        ! viscous part of forces
        NdS = (2._R_P*NFiS(ni1,j,k) + NFiS(ni1+1,j,k) + NFiS(ni1-1,j,k))/(2._R_P*(volume(ni1,j,k)+volume(ni1+1,j,k)))
        f_v(ci1,j,k) = (momentum(ci1-1+face,j,k) - momentum(ci1-2+face,j,k))*(NFiS(ni1,j,k).dot.NdS)/pp%Re

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
        pp%Ssum = pp%Ssum + Si(ni1,j,k)
      enddo
    enddo
    !$OMP END PARALLEL
    if (.not.pp%cell) then ! extrapolating ghost cell values
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
        if (ticc(i,nj1-1+(face-2),k)/=pp%patch) cycle
        if ( icc(i,cj1,k)/=0) cycle
        if(pp%level_set) then
          if (f0(i,nj1-1+(face-2),k)>0._R_P) cycle
        endif

        ! patch center coordinates
        fco = 0.25_R_P*(node(i,nj1,k)+node(i-1,nj1,k)+node(i,nj1,k-1)+node(i-1,nj1,k-1))

        ! pressure part of forces
        hyd = 0._R_P
        if(pp%level_set) hyd = (fco%z - pp%zfs)*pp%rFr2
        f_p(i,cj1,k) = -(1.5_R_P*pressure(i,cj1,k) - 0.5_R_P*pressure(i,cj1+3-2*(face-2),k) - hyd)*NFjS(i,nj1,k)
        ! viscous part of forces
        NdS = (2._R_P*NFjS(i,nj1,k) + NFjS(i,nj1+1,k) + NFjS(i,nj1-1,k))/(2._R_P*(volume(i,nj1,k)+volume(i,nj1+1,k)))
        f_v(i,cj1,k) = (momentum(i,cj1-1+(face-2),k) - momentum(i,cj1-2+(face-2),k))*(NFjS(i,nj1,k).dot.NdS)/pp%Re

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
        pp%Ssum = pp%Ssum + Sj(i,nj1,k)
      enddo
    enddo
    !$OMP END PARALLEL
    if (.not.pp%cell) then ! extrapolating ghost cell values
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
        if (ticc(i,j,nk1-1+(face-4))/=pp%patch) cycle
        if ( icc(i,j,ck1)/=0) cycle
        if(pp%level_set) then
          if (f0(i,j,nk1-1+(face-4))>0._R_P) cycle
        endif

        ! patch center coordinates
        fco = 0.25_R_P*(node(i,j,nk1)+node(i-1,j,nk1)+node(i,j-1,nk1)+node(i-1,j-1,nk1))

        ! pressure part of forces
        hyd = 0._R_P
        if(pp%level_set) hyd = (fco%z - pp%zfs)*pp%rFr2
        f_p(i,j,ck1) = -(1.5_R_P*pressure(i,j,ck1) - 0.5_R_P*pressure(i,j,ck1+3-2*(face-4)) - hyd)*NFkS(i,j,nk1)
        ! viscous part of forces
        NdS = (2._R_P*NFkS(i,j,nk1) + NFkS(i,j,nk1+1) + NFkS(i,j,nk1-1))/(2._R_P*(volume(i,j,nk1)+volume(i,j,nk1+1)))
        f_v(i,j,ck1) = (momentum(i,j,ck1-1+(face-4)) - momentum(i,j,ck1-2+(face-4)))*(NFkS(i,j,nk1).dot.NdS)/pp%Re

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
        pp%Ssum = pp%Ssum + Sk(i,j,nk1)
      enddo
    enddo
    !$OMP END PARALLEL
    if (.not.pp%cell) then ! extrapolating ghost cell values
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
  pp%fsum_p = pp%fsum_p + (ex*fsumpx+ey*fsumpy+ez*fsumpz)
  pp%fsum_v = pp%fsum_v + (ex*fsumvx+ey*fsumvy+ez*fsumvz)
  pp%msum_p = pp%msum_p + (ex*msumpx+ey*msumpy+ez*msumpz)
  pp%msum_v = pp%msum_v + (ex*msumvx+ey*msumvy+ez*msumvz)
  if (pp%Ng>0) then
     write(pp%unit_g_for(pp%blockmap(1,b)),'(A)')trim(str(n=(fsumvx+fsumpx)))//' '// &
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
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine compute_forces

  !> @brief Procedure for computing yplus variables.
  subroutine compute_yp(pp,face,ni1,ni2,nj1,nj2,nk1,nk2,ci1,ci2,cj1,cj2,ck1,ck2)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_PostProcess), intent(IN)::    pp                      !< Post-processor data.
  integer(I4P),           intent(IN)::    face                    !< Face where patch is defined: 1,2,3,4,5,6.
  integer(I4P),           intent(IN)::    ni1,ni2                 !< First and last node i indexes.
  integer(I4P),           intent(IN)::    nj1,nj2                 !< First and last node j indexes.
  integer(I4P),           intent(IN)::    nk1,nk2                 !< First and last node k indexes.
  integer(I4P),           intent(IN)::    ci1,ci2                 !< First and last cell i indexes.
  integer(I4P),           intent(IN)::    cj1,cj2                 !< First and last cell j indexes.
  integer(I4P),           intent(IN)::    ck1,ck2                 !< First and last cell k indexes.
  type(Type_Vector)::                     NdS                     !< Normal "tilde" for viscous part of forces (or distance for y+).
  integer(I4P)::                          i,j,k                   !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  select case(face)
  case(1,2)
    !$OMP PARALLEL DEFAULT(NONE) &
    !$OMP PRIVATE(i,j,k,NdS)     &
    !$OMP SHARED(ni1,ni2,nj1,nj2,nk1,nk2,ci1,ci2,cj1,cj2,ck1,ck2,face,patch,node,NFiS,icc,ticc,momentum,f0,Re,yplus)
    !$OMP DO
    do k=ck1,ck2
      do j=cj1,cj2
        ! checking if this is an active cell
        if (ticc(ni1-1+face,j,k)/=pp%patch) cycle
        if ( icc(ci1,j,k)/=0) cycle
        if(pp%level_set) then
          if (f0(ni1-1+face,j,k)>0._R_P) cycle
        endif

        ! distance from the patch
        NdS = 0.25_R_P*(node(ni1,j,k)+node(ni1,j-1,k)+node(ni1,j,k-1)+node(ni1,j-1,k-1)) - &
              0.25_R_P*(node(ci1,j,k)+node(ci1,j-1,k)+node(ci1,j,k-1)+node(ci1,j-1,k-1))
        yplus(ci1,j,k) = sqrt(pp%Re*normL2(momentum(ci1,j,k).ortho.NFiS(ni1,j,k))*(abs(NdS.dot.NFiS(ni1,j,k))))
      enddo
    enddo
    !$OMP END PARALLEL
    if (.not.pp%cell) then ! extrapolating ghost cell values
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
        if (ticc(i,nj1-1+(face-2),k)/=pp%patch) cycle
        if ( icc(i,cj1,k)/=0) cycle
        if(pp%level_set) then
          if (f0(i,nj1-1+(face-2),k)>0._R_P) cycle
        endif

        ! distance from the patch
        NdS = 0.25_R_P*(node(i,nj1,k)+node(i-1,nj1,k)+node(i,nj1,k-1)+node(i-1,nj1,k-1)) - &
              0.25_R_P*(node(i,cj1,k)+node(i-1,cj1,k)+node(i,cj1,k-1)+node(i-1,cj1,k-1))
        yplus(i,cj1,k) = sqrt(pp%Re*normL2(momentum(i,cj1,k).ortho.NFjS(i,nj1,k))*(abs(NdS.dot.NFjS(i,nj1,k))))
      enddo
    enddo
    !$OMP END PARALLEL
    if (.not.pp%cell) then ! extrapolating ghost cell values
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
        if (ticc(i,j,nk1-1+(face-4))/=pp%patch) cycle
        if ( icc(i,j,ck1)/=0) cycle
        if(pp%level_set) then
          if (f0(i,j,nk1-1+(face-4))>0._R_P) cycle
        endif

        ! distance from the patch
        NdS = 0.25_R_P*(node(i,j,nk1)+node(i-1,j,nk1)+node(i,j-1,nk1)+node(i-1,j-1,nk1)) - &
              0.25_R_P*(node(i,j,ck1)+node(i-1,j,ck1)+node(i,j-1,ck1)+node(i-1,j-1,ck1))
        yplus(i,j,ck1) = sqrt(pp%Re*normL2(momentum(i,j,ck1).ortho.NFkS(i,j,nk1))*(abs(NdS.dot.NFkS(i,j,nk1))))
      enddo
    enddo
    !$OMP END PARALLEL
    if (.not.pp%cell) then ! extrapolating ghost cell values
      yplus(ni1  ,:    ,ck1) = yplus(ci1,:  ,ck1)
      yplus(ni2+1,:    ,ck1) = yplus(ni2,:  ,ck1)
      yplus(:    ,nj1  ,ck1) = yplus(:  ,cj1,ck1)
      yplus(:    ,nj2+1,ck1) = yplus(:  ,nj2,ck1)
    endif
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine compute_yp

  !> @brief Procedure for computing yplus variables.
  subroutine compute_tau(pp,face,ni1,ni2,nj1,nj2,nk1,nk2,ci1,ci2,cj1,cj2,ck1,ck2)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_PostProcess), intent(INOUT):: pp                      !< Post-processor data.
  integer(I4P),           intent(IN)::    face                    !< Face where patch is defined: 1,2,3,4,5,6.
  integer(I4P),           intent(IN)::    ni1,ni2                 !< First and last node i indexes.
  integer(I4P),           intent(IN)::    nj1,nj2                 !< First and last node j indexes.
  integer(I4P),           intent(IN)::    nk1,nk2                 !< First and last node k indexes.
  integer(I4P),           intent(IN)::    ci1,ci2                 !< First and last cell i indexes.
  integer(I4P),           intent(IN)::    cj1,cj2                 !< First and last cell j indexes.
  integer(I4P),           intent(IN)::    ck1,ck2                 !< First and last cell k indexes.
  type(Type_Vector)::                     NdS                     !< Normal "tilde" for viscous part of forces (or distance for y+).
  integer(I4P)::                          i,j,k                   !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
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
        if (ticc(ni1-1+face,j,k)/=pp%patch) cycle
        if ( icc(ci1,j,k)/=0) cycle
        if(pp%level_set) then
          if (f0(ni1-1+face,j,k)>0._R_P) cycle
        endif

        ! viscous part of forces
        NdS = (2._R_P*NFiS(ni1,j,k) + NFiS(ni1+1,j,k) + NFiS(ni1-1,j,k))/(2._R_P*(volume(ni1,j,k)+volume(ni1+1,j,k)))
        tau(ci1,j,k) = (momentum(ci1-1+face,j,k) - momentum(ci1-2+face,j,k))*(NFi(ni1,j,k).dot.NdS)/pp%Re
        tau(ci1,j,k) = (tau(ci1,j,k).ortho.(NFi(ni1,j,k)))

        if (face==2) then
          tau(ci1,j,k) = -tau(ci1,j,k)
        endif
      enddo
    enddo
    !$OMP END PARALLEL
    if (.not.pp%cell) then ! extrapolating ghost cell values
      tau(ci1,nj1  ,:    ) = tau(ci1,cj1,:  )
      tau(ci1,nj2+1,:    ) = tau(ci1,nj2,:  )
      tau(ci1,:    ,nk1  ) = tau(ci1,:  ,ck1)
      tau(ci1,:    ,nk2+1) = tau(ci1,:  ,nk2)
      ! corners
      tau(ci1,nj1,nk1) = tau(ci1,cj1,ck1)
      tau(ci1,nj1,nk2) = tau(ci1,cj1,ck2)
      tau(ci1,nj2,nk1) = tau(ci1,cj2,ck1)
      tau(ci1,nj2,nk2) = tau(ci1,cj2,ck2)
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
        if (ticc(i,nj1-1+(face-2),k)/=pp%patch) cycle
        if ( icc(i,cj1,k)/=0) cycle
        if(pp%level_set) then
          if (f0(i,nj1-1+(face-2),k)>0._R_P) cycle
        endif

        ! viscous part of forces
        NdS = (2._R_P*NFjS(i,nj1,k) + NFjS(i,nj1+1,k) + NFjS(i,nj1-1,k))/(2._R_P*(volume(i,nj1,k)+volume(i,nj1+1,k)))
        tau(i,cj1,k) = (momentum(i,cj1-1+(face-2),k) - momentum(i,cj1-2+(face-2),k))*(NFj(i,nj1,k).dot.NdS)/pp%Re
        tau(i,cj1,k) = (tau(i,cj1,k).ortho.(NFj(i,nj1,k)))

        if (face==4) then
          tau(i,cj1,k) = -tau(i,cj1,k)
        endif
      enddo
    enddo
    !$OMP END PARALLEL
    if (.not.pp%cell) then ! extrapolating ghost cell values
      tau(ni1  ,cj1,:    ) = tau(ci1,cj1,:  )
      tau(ni2+1,cj1,:    ) = tau(ni2,cj1,:  )
      tau(:    ,cj1,nk1  ) = tau(:  ,cj1,ck1)
      tau(:    ,cj1,nk2+1) = tau(:  ,cj1,nk2)
      ! corners
      tau(ni1,cj1,nk1) = tau(ci1,cj1,ck1)
      tau(ni1,cj1,nk2) = tau(ci1,cj1,ck2)
      tau(ni2,cj1,nk1) = tau(ci2,cj1,ck1)
      tau(ni2,cj1,nk2) = tau(ci2,cj1,ck2)
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
        if (ticc(i,j,nk1-1+(face-4))/=pp%patch) cycle
        if ( icc(i,j,ck1)/=0) cycle
        if(pp%level_set) then
          if (f0(i,j,nk1-1+(face-4))>0._R_P) cycle
        endif

        ! viscous part of forces
        NdS = (2._R_P*NFkS(i,j,nk1) + NFkS(i,j,nk1+1) + NFkS(i,j,nk1-1))/(2._R_P*(volume(i,j,nk1)+volume(i,j,nk1+1)))
        tau(i,j,ck1) = (momentum(i,j,ck1-1+(face-4)) - momentum(i,j,ck1-2+(face-4)))*(NFk(i,j,nk1).dot.NdS)/pp%Re
        tau(i,j,ck1) = (tau(i,j,ck1).ortho.(NFk(i,j,nk1)))

        if (face==6) then
          tau(i,j,ck1) = -tau(i,j,ck1)
        endif
      enddo
    enddo
    !$OMP END PARALLEL
    if (.not.pp%cell) then ! extrapolating ghost cell values
      tau(ni1  ,:    ,ck1) = tau(ci1,:  ,ck1)
      tau(ni2+1,:    ,ck1) = tau(ni2,:  ,ck1)
      tau(:    ,nj1  ,ck1) = tau(:  ,cj1,ck1)
      tau(:    ,nj2+1,ck1) = tau(:  ,nj2,ck1)
      ! corners
      tau(ni1,nj1,ck1) = tau(ci1,cj1,ck1)
      tau(ni1,nj2,ck1) = tau(ci1,cj2,ck1)
      tau(ni2,nj1,ck1) = tau(ci2,cj1,ck1)
      tau(ni2,nj2,ck1) = tau(ci2,cj2,ck1)
    endif
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine compute_tau

  !> @brief Procedure for interpolating cell centered variable into node centered one.
  subroutine varinterpolation(ni1,ni2,nj1,nj2,nk1,nk2,face,var,vari)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(IN)::  ni1,ni2,nj1,nj2,nk1,nk2             ! Nodes bounds.
  integer(I4P), intent(IN)::  face                                ! Actual face.
  real(R8P),    intent(IN)::  var (ni1:ni2+1,nj1:nj2+1,nk1:nk2+1) ! Cell centered variable.
  real(R8P),    intent(OUT):: vari(ni1:ni2  ,nj1:nj2  ,nk1:nk2  ) ! Node centered interpolated variable.
  integer(I4P)::              i,j,k                               ! Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  do k=nk1,nk2
    do j=nj1,nj2
      do i=ni1,ni2
        vari(i,j,k) =  var(i+1,j+1,k+1) + var(i+1,j+1,k)  &
                     + var(i  ,j+1,k+1) + var(i  ,j+1,k)  &
                     + var(i+1,j  ,k+1) + var(i+1,j  ,k)  &
                     + var(i  ,j  ,k+1) + var(i  ,j  ,k)
        vari(i,j,k) = 0.125_R_P*vari(i,j,k)
      enddo
    enddo
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine varinterpolation
endmodule Block_Variables
