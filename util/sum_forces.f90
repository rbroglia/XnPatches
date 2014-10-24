program sum_forces
!-----------------------------------------------------------------------------------------------------------------------------------
USE, intrinsic:: ISO_FORTRAN_ENV, only: stdout => OUTPUT_UNIT, stderr => ERROR_UNIT ! Standard output/error logical units.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
integer, parameter::       R8P = selected_real_kind(15,307) !< 15  digits, range \f$[10^{-307} , 10^{+307}  - 1]\f$; 64 bits.
integer, parameter::       I4P = selected_int_kind(9)       !< Range \f$[-2^{31},+2^{31} - 1]\f$, 10 digits plus sign; 32 bits.
integer, parameter::       ui = 10                          !< Unit file.
character(10), parameter:: FR8P = '(E23.15E3)'              !< Output format for kind=R8P variable.
real(R8P)::                fdata(1:19),ddata(1:19)          !< Forces and dummy data.
character(100)::           froot                            !< Root of file names.
character(100)::           fname                            !< Dummy file name.
integer(I4P)::             Np                               !< Number of procs(files).
integer(I4P)::             Nca = 0                          !< Number of command line arguments.
character(8)::             switch                           !< Switch string.
logical::                  is_file                          !< Flag for inquiring files presence.
integer(I4P)::             p,c                              !< Counter.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
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
  case('-fr')
    call get_command_argument(c+1,froot) ; c = c + 1
  case('-np')
    call get_command_argument(c+1,switch) ; c = c + 1
    read(switch,*)Np
  case default
    write(stderr,'(A)') ' Unknown switch '
    write(stderr,'(A)') adjustl(trim(switch))
    write(stderr,*)
    call print_usage
    stop
  endselect
enddo
fdata=0._R8P
do p=0,Np-1
  ddata=0._R8P
  write(switch,'(I3.3)')p ; fname = trim(froot)//'.p'//trim(switch)//'-forces.dat' ; inquire(file=trim(fname),exist=is_file)
  if (.not.is_file) cycle ; open(unit=ui,file=trim(fname)) ; read(ui,*)ddata ; close(ui) ; fdata=fdata+ddata
enddo
write(stdout,'(19('//FR8P//',1X))')fdata
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  subroutine print_usage()
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(stdout,*)
  write(stdout,'(A)') ' sum_forces'
  write(stdout,'(A)') ' Sum forces from XnPatches outputs'
  write(stdout,'(A)') ' Usage:'
  write(stdout,'(A)') '   sum_forces -fr root_of_file_names -nf #number_of_files'
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine print_usage
endprogram sum_forces
