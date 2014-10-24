!> @brief This module contains procedures for computing Xnavis mesh metrics.
module Lib_Metrics
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision     ! Integers and reals precision definition.
USE Block_Variables  ! Block variables definition.
USE Data_Type_Vector ! Definition of Type_Vector.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
public:: compute_metrics
public:: bc_metrics_correction
!-----------------------------------------------------------------------------------------------------------------------------------
contains
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
endmodule Lib_Metrics
