program main
  use iso_fortran_env, only: output_unit
  use MPI
  use OMP_LIB

  use SsTC_driver_comms
  use SsTC_driver

  implicit none

  character(len=10) :: ver = "0.0.0"

  call MPI_INIT(ierror)

  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)

  if (rank == 0) write (unit=output_unit, fmt="(A)") ""
  if (rank == 0) write (unit=output_unit, fmt="(A)") "   /$$$$$$            /$$$$$$$$  /$$$$$$  "
  if (rank == 0) write (unit=output_unit, fmt="(A)") "  /$$__  $$          |__  $$__/ /$$__  $$ "
  if (rank == 0) write (unit=output_unit, fmt="(A)") " | $$  \__/  /$$$$$$$   | $$   | $$  \__/ "
  if (rank == 0) write (unit=output_unit, fmt="(A)") " |  $$$$$$  /$$_____/   | $$   | $$       "
  if (rank == 0) write (unit=output_unit, fmt="(A)") "  \____  $$|  $$$$$$    | $$   | $$       "
  if (rank == 0) write (unit=output_unit, fmt="(A)") "  /$$  \ $$ \____  $$   | $$   | $$    $$ "
  if (rank == 0) write (unit=output_unit, fmt="(A)") " |  $$$$$$/ /$$$$$$$/   | $$   |  $$$$$$/ "
  if (rank == 0) write (unit=output_unit, fmt="(A)") "  \______/ |_______/    |__/    \______/  "
  if (rank == 0) write (unit=output_unit, fmt="(A)") ""

  if (rank == 0) write (unit=output_unit, fmt="(A)") &
    " Solid-state Task Constructor Driver (SsTC_driver) v"//trim(adjustl(ver))//" built."

  call MPI_FINALIZE(ierror)

end program main
