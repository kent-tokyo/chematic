! Independent check of explicit-decimal HISTORY fields under the manual's
! F12.6, E12.4 and G12.4 input edit descriptors. No DL_POLY code is included.
program real_field_probe
  use iso_fortran_env, only: real64, error_unit
  implicit none
  character(len=12) :: field
  real(real64) :: x, y, z
  integer :: status
  do
    read(*,'(a12)',iostat=status) field
    if (status < 0) exit
    if (status /= 0) stop 1
    read(field,'(f12.6)',iostat=status) x
    if (status /= 0) stop 2
    read(field,'(e12.4)',iostat=status) y
    if (status /= 0) stop 3
    read(field,'(g12.4)',iostat=status) z
    if (status /= 0) stop 4
    write(*,'(3(es26.17e3,1x))') x, y, z
  end do
end program real_field_probe
