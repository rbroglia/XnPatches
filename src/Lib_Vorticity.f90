!> @brief This module contains procedures for computing vorticity-related variables.
module Lib_Vorticity
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision     ! Integers and reals precision definition.
USE Block_Variables  ! Block variables definition.
USE Data_Type_Vector ! Definition of Type_Vector.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
public:: compute_vorticity
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @brief Procedure for computing vorticity variables.
  subroutine compute_vorticity(gc,Ni,Nj,Nk)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(IN)::  gc(1:6)                                                        !< Number of ghost cells.
  integer(I4P), intent(IN)::  Ni,Nj,Nk                                                       !< Number of cells.
  real(R8P)::                 Fi(1:3,1:3,0-gc(1):Ni+gc(2),0-gc(3):Nj+gc(4),0-gc(5):Nk+gc(6)) !< Fluxes i direction.
  real(R8P)::                 Fj(1:3,1:3,0-gc(1):Ni+gc(2),0-gc(3):Nj+gc(4),0-gc(5):Nk+gc(6)) !< Fluxes j direction.
  real(R8P)::                 Fk(1:3,1:3,0-gc(1):Ni+gc(2),0-gc(3):Nj+gc(4),0-gc(5):Nk+gc(6)) !< Fluxes k direction.
  type(Type_Vector)::         um                                                             !< Dummy vector variables.
  real(R8P)::                 c(0:2),emin,emax,eval,fval,mu                                  !< Dummy reals.
  real(R8P), dimension(3,3):: IDEN,G,S,O                                                     !< Matrices.
  integer(I4P)::              i,j,k,ii,jj,kk,iter                                            !< Counters.
  real(R8P), parameter ::     eps6=1d-6, eps9=1d-9                                           !< Tolerances.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  Fi   = 0._R_P
  Fj   = 0._R_P
  Fk   = 0._R_P
  vord = 0._R_P
  qfac = 0._R_P
  heli = 0._R_P
  IDEN = 0._R_P
  do i=1,3
     IDEN(i,i) = 1._R_P
  end do
  ! extrapolating momentum from inner cells to ghost cells
  do k=1,Nk
    do j=1,Nj
      momentum(1- gc(1),j,k) = 2._R_P*momentum(0,   j,k)-momentum(1, j,k)
      momentum(Ni+gc(2),j,k) = 2._R_P*momentum(Ni+1,j,k)-momentum(Ni,j,k)
    enddo
  enddo
  do k=1,Nk
    do i=0,Ni+1
      momentum(i,1- gc(3),k) = 2._R_P*momentum(i,0,   k)-momentum(i,1, k)
      momentum(i,Nj+gc(4),k) = 2._R_P*momentum(i,Nj+1,k)-momentum(i,Nj,k)
    enddo
  enddo
  do j=0,Nj+1
    do i=0,Ni+1
      momentum(i,j,1- gc(5)) = 2._R_P*momentum(i,j,0   )-momentum(i,j,1 )
      momentum(i,j,Nk+gc(6)) = 2._R_P*momentum(i,j,Nk+1)-momentum(i,j,Nk)
    enddo
  enddo
  ! computing fluxes
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(i,j,k,um)      &
  !$OMP SHARED(Ni,Nj,Nk,Fi,Fj,Fk,NFiS,NFjS,NFkS,momentum)
  !$OMP DO
  do k=0,Nk+1
    do j=0,Nj+1
      do i=-1,Ni+1
        um = 0.5_R_P*(momentum(i,j,k)+momentum(i+1,j,k))
        Fi(1,1,i,j,k) = um%x*NFiS(i,j,k)%x
        Fi(1,2,i,j,k) = um%x*NFiS(i,j,k)%y
        Fi(1,3,i,j,k) = um%x*NFiS(i,j,k)%z
        Fi(2,1,i,j,k) = um%y*NFiS(i,j,k)%x
        Fi(2,2,i,j,k) = um%y*NFiS(i,j,k)%y
        Fi(2,3,i,j,k) = um%y*NFiS(i,j,k)%z
        Fi(3,1,i,j,k) = um%z*NFiS(i,j,k)%x
        Fi(3,2,i,j,k) = um%z*NFiS(i,j,k)%y
        Fi(3,3,i,j,k) = um%z*NFiS(i,j,k)%z
      enddo
    enddo
  enddo
  !$OMP DO
  do k=0,Nk+1
    do j=-1,Nj+1
      do i=0,Ni+1
        um = 0.5_R_P*(momentum(i,j,k)+momentum(i,j+1,k))
        Fj(1,1,i,j,k) = um%x*NFjS(i,j,k)%x
        Fj(1,2,i,j,k) = um%x*NFjS(i,j,k)%y
        Fj(1,3,i,j,k) = um%x*NFjS(i,j,k)%z
        Fj(2,1,i,j,k) = um%y*NFjS(i,j,k)%x
        Fj(2,2,i,j,k) = um%y*NFjS(i,j,k)%y
        Fj(2,3,i,j,k) = um%y*NFjS(i,j,k)%z
        Fj(3,1,i,j,k) = um%z*NFjS(i,j,k)%x
        Fj(3,2,i,j,k) = um%z*NFjS(i,j,k)%y
        Fj(3,3,i,j,k) = um%z*NFjS(i,j,k)%z
      enddo
    enddo
  enddo
  !$OMP DO
  do k=-1,Nk+1
    do j=0,Nj+1
      do i=0,Ni+1
        um = 0.5_R_P*(momentum(i,j,k)+momentum(i,j,k+1))
        Fk(1,1,i,j,k) = um%x*NFkS(i,j,k)%x
        Fk(1,2,i,j,k) = um%x*NFkS(i,j,k)%y
        Fk(1,3,i,j,k) = um%x*NFkS(i,j,k)%z
        Fk(2,1,i,j,k) = um%y*NFkS(i,j,k)%x
        Fk(2,2,i,j,k) = um%y*NFkS(i,j,k)%y
        Fk(2,3,i,j,k) = um%y*NFkS(i,j,k)%z
        Fk(3,1,i,j,k) = um%z*NFkS(i,j,k)%x
        Fk(3,2,i,j,k) = um%z*NFkS(i,j,k)%y
        Fk(3,3,i,j,k) = um%z*NFkS(i,j,k)%z
      enddo
    enddo
  enddo
  !$OMP END PARALLEL
  do k=0,Nk+1
    do j=0,Nj+1
      do i=0,Ni+1
         ! computing velocity gradient
         do jj=1,3
            do ii=1,3
               G(ii,jj) = FI(ii,jj,i,j,k)-FI(ii,jj,i-1,j,k)&
                        + FJ(ii,jj,i,j,k)-FJ(ii,jj,i,j-1,k)&
                        + FK(ii,jj,i,j,k)-FK(ii,jj,i,j,k-1)
               G(ii,jj) = G(ii,jj)/volume(i,j,k)
            enddo
         enddo

         ! computing vorticity vector
         um%x = G(3,2) - G(2,3)
         um%y = G(1,3) - G(3,1) ! control signs!!!!
         um%z = G(2,1) - G(1,2)

         ! tensor S^2 + O^2 (saved in G)
         S = 0._R_P
         O = 0._R_P
         do kk=1,3
            do jj=1,3
               S(jj,kk) = 0.5_R_P*(G(jj,kk)+G(kk,jj))
               O(jj,kk) = 0.5_R_P*(G(jj,kk)-G(kk,jj))
            enddo
         enddo
         G = matmul(S,S) + matmul(O,O)

         ! coefficients of characterist polynomial: lamda^3 + c(2)*lamda^2 + c(1)*lamda + c(0) = 0
         c(2) = -(G(1,1) + G(2,2) + G(3,3))
         c(1) = G(1,1)*G(2,2) + G(1,1)*G(3,3) + G(2,2)*G(3,3) - G(2,3)**2 - G(1,3)**2 - G(1,2)**2
         c(0) = G(1,1)*G(2,3)**2 + G(2,2)*G(1,3)**2 + G(3,3)*G(1,2)**2 - 2._R_P*G(2,3)*G(1,3)*G(1,2) - G(1,1)*G(2,2)*G(3,3)

         ! computing second eigenvalue of characteristic polynomial
         mu = sqrt(c(2)**2 - 3._R_P*c(1))
         emin = (-c(2)-mu)/3._R_P
         emax = (-c(2)+mu)/3._R_P
         do iter=1,100
            eval = 0.5_R_P*(emin+emax)
            fval = eval**3 + c(2)*eval**2 + c(1)*eval + c(0)
            if (fval<0._R_P) then
               emax = eval
            else
               emin = eval
            end if
            if (abs(fval)<eps9 .and.((emax-emin)/eval)<eps6) exit
         end do

         ! saving vordet, qfactor and helicity
         vord(     i,j,k) = eval
         qfac(     i,j,k) = 0.5_R_P*(dot_product(O(1,:),O(1,:)) + dot_product(O(2,:),O(2,:)) + dot_product(O(3,:),O(3,:)) - &
                                    (dot_product(S(1,:),S(1,:)) + dot_product(S(2,:),S(2,:)) + dot_product(S(3,:),S(3,:))))
         heli(     i,j,k) = momentum(i,j,k).dot.um / (max(eps6*eps6,normL2(momentum(i,j,k))*normL2(um)))
         vorticity(i,j,k) = um
      enddo
    enddo
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine compute_vorticity
endmodule Lib_Vorticity
