!> @brief Assignment procedures definitions.
#if defined self_type_ && defined ass_scal_
  !> @brief Procedure for assignment between a scalar (real R16P) and self.
  elemental subroutine assign_ScalR16P(self,scal)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  self_type_
  real(R16P),         intent(IN)::    scal
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ass_scal_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_ScalR16P

  !> @brief Procedure for assignment between a scalar (real R8P) and self.
  elemental subroutine assign_ScalR8P(self,scal)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  self_type_
  real(R8P),          intent(IN)::    scal
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ass_scal_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_ScalR8P

  !> @brief Procedure for assignment between a scalar (real R4P) and self.
  elemental subroutine assign_ScalR4P(self,scal)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  self_type_
  real(R4P),          intent(IN)::    scal
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ass_scal_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_ScalR4P

  !> @brief Procedure for assignment between a scalar (integer I8P) and self.
  elemental subroutine assign_ScalI8P(self,scal)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  self_type_
  integer(I8P),       intent(IN)::    scal
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ass_scal_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_ScalI8P

  !> @brief Procedure for assignment between a scalar (integer I4P) and self.
  elemental subroutine assign_ScalI4P(self,scal)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  self_type_
  integer(I4P),       intent(IN)::    scal
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ass_scal_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_ScalI4P

  !> @brief Procedure for assignment between a scalar (integer I2P) and self.
  elemental subroutine assign_ScalI2P(self,scal)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  self_type_
  integer(I2P),       intent(IN)::    scal
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ass_scal_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_ScalI2P

  !> @brief Procedure for assignment between a scalar (integer I1P) and self.
  elemental subroutine assign_ScalI1P(self,scal)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  self_type_
  integer(I1P),       intent(IN)::    scal
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ass_scal_
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine assign_ScalI1P
#endif
