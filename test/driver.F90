program driver

  use MPI
  use OMP_LIB

  use, intrinsic :: iso_fortran_env, only: error_unit
  use testdrive, only: error_type
  !Tests:
  use Sampling_Suite, only: sample_no_sum_par_1, sample_no_sum_par_2, &
    sample_no_sum_par_3, sample_no_sum_par_4, &
    sample_sum_par_1, sample_sum_par_2, &
    sample_sum_par_3, sample_sum_par_4, &
    sample_no_sum_s_par_1, sample_no_sum_s_par_2, &
    sample_no_sum_s_par_3, sample_no_sum_s_par_4, &
    sample_sum_s_par_1, sample_sum_s_par_2, &
    sample_sum_s_par_3, sample_sum_s_par_4
  use Randomized_Suite, only: sample_rnd_no_sum, sample_rnd_sum
  use Utilities_Suite, only: test_kpath, test_kslice, test_crys_to_cart_cart_to_crys
  use Material_Properties_Suite, only: sample_and_sum_GaAs_optcond_in_peak

  implicit none

  integer :: ierror, rank

  type(error_type), allocatable :: error

  call MPI_INIT(ierror)

  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)

  if (rank == 0) write (error_unit, "(A)") "Suite: Sampling Suite:"

  if (rank == 0) write (error_unit, "(A)") &
    "Test (1/16): Testing sampling for given int. and cont. ind. OMP+MPI paralleization."
  call sample_no_sum_par_1(error)
  if (rank == 0) then
    if (allocated(error)) error stop "Test (1/16): FAILED."
  endif
  if (rank == 0) write (error_unit, "(A)") "Test (1/16): PASSED."

  if (rank == 0) write (error_unit, "(A)") &
    "Test (2/16): Testing sampling for given int. and cont. ind. MPI paralleization."
  call sample_no_sum_par_2(error)
  if (rank == 0) then
    if (allocated(error)) error stop "Test (2/16): FAILED."
  endif
  if (rank == 0) write (error_unit, "(A)") "Test (2/16): PASSED."

  if (rank == 0) write (error_unit, "(A)") &
    "Test (3/16): Testing sampling for given int. and cont. ind. OMP paralleization."
  call sample_no_sum_par_3(error)
  if (rank == 0) then
    if (allocated(error)) error stop "Test (3/16): FAILED."
  endif
  if (rank == 0) write (error_unit, "(A)") "Test (3/16): PASSED."

  if (rank == 0) write (error_unit, "(A)") &
    "Test (4/16): Testing sampling for given int. and cont. ind. No paralleization."
  call sample_no_sum_par_4(error)
  if (rank == 0) then
    if (allocated(error)) error stop "Test (4/16): FAILED."
  endif
  if (rank == 0) write (error_unit, "(A)") "Test (4/16): PASSED."

  if (rank == 0) write (error_unit, "(A)") &
    "Test (5/16): Testing sampling for given int. and cont. ind. OMP+MPI paralleization. Implied sum."
  call sample_sum_par_1(error)
  if (rank == 0) then
    if (allocated(error)) error stop "Test (5/16): FAILED."
  endif
  if (rank == 0) write (error_unit, "(A)") "Test (5/16): PASSED."

  if (rank == 0) write (error_unit, "(A)") &
    "Test (6/16): Testing sampling for given int. and cont. ind. MPI paralleization. Implied sum."
  call sample_sum_par_2(error)
  if (rank == 0) then
    if (allocated(error)) error stop "Test (6/16): FAILED."
  endif
  if (rank == 0) write (error_unit, "(A)") "Test (6/16): PASSED."

  if (rank == 0) write (error_unit, "(A)") &
    "Test (7/16): Testing sampling for given int. and cont. ind. OMP paralleization. Implied sum."
  call sample_sum_par_3(error)
  if (rank == 0) then
    if (allocated(error)) error stop "Test (7/16): FAILED."
  endif
  if (rank == 0) write (error_unit, "(A)") "Test (7/16): PASSED."

  if (rank == 0) write (error_unit, "(A)") &
    "Test (8/16): Testing sampling for given int. and cont. ind. No paralleization. Implied sum."
  call sample_sum_par_4(error)
  if (rank == 0) then
    if (allocated(error)) error stop "Test (8/16): FAILED."
  endif
  if (rank == 0) write (error_unit, "(A)") "Test (8/16): PASSED."

  if (rank == 0) write (error_unit, "(A)") &
    "Test (9/16): Testing sampling when int. and cont. ind. are not given. OMP+MPI paralleization."
  call sample_no_sum_s_par_1(error)
  if (rank == 0) then
    if (allocated(error)) error stop "Test (9/16): FAILED."
  endif
  if (rank == 0) write (error_unit, "(A)") "Test (9/16): PASSED."

  if (rank == 0) write (error_unit, "(A)") &
    "Test (10/16): Testing sampling when int. and cont. ind. are not given. MPI paralleization."
  call sample_no_sum_s_par_2(error)
  if (rank == 0) then
    if (allocated(error)) error stop "Test (10/16): FAILED."
  endif
  if (rank == 0) write (error_unit, "(A)") "Test (10/16): PASSED."

  if (rank == 0) write (error_unit, "(A)") &
    "Test (11/16): Testing sampling when int. and cont. ind. are not given. OMP paralleization."
  call sample_no_sum_s_par_3(error)
  if (rank == 0) then
    if (allocated(error)) error stop "Test (11/16): FAILED."
  endif
  if (rank == 0) write (error_unit, "(A)") "Test (11/16): PASSED."

  if (rank == 0) write (error_unit, "(A)") &
    "Test (12/16): Testing sampling when int. and cont. ind. are not given. No paralleization."
  call sample_no_sum_s_par_4(error)
  if (rank == 0) then
    if (allocated(error)) error stop "Test (12/16): FAILED."
  endif
  if (rank == 0) write (error_unit, "(A)") "Test (12/16): PASSED."

  if (rank == 0) write (error_unit, "(A)") &
    "Test (13/16): Testing sampling when int. and cont. ind. are not given. OMP+MPI paralleization. Implied sum."
  call sample_sum_s_par_1(error)
  if (rank == 0) then
    if (allocated(error)) error stop "Test (13/16): FAILED."
  endif
  if (rank == 0) write (error_unit, "(A)") "Test (13/16): PASSED."

  if (rank == 0) write (error_unit, "(A)") &
    "Test (14/16): Testing sampling when int. and cont. ind. are not given. MPI paralleization. Implied sum."
  call sample_sum_s_par_2(error)
  if (rank == 0) then
    if (allocated(error)) error stop "Test (14/16): FAILED."
  endif
  if (rank == 0) write (error_unit, "(A)") "Test (14/16): PASSED."

  if (rank == 0) write (error_unit, "(A)") &
    "Test (15/16): Testing sampling when int. and cont. ind. are not given. OMP paralleization. Implied sum."
  call sample_sum_s_par_3(error)
  if (rank == 0) then
    if (allocated(error)) error stop "Test (15/16): FAILED."
  endif
  if (rank == 0) write (error_unit, "(A)") "Test (15/16): PASSED."

  if (rank == 0) write (error_unit, "(A)") &
    "Test (16/16): Testing sampling when int. and cont. ind. are not given. No paralleization. Implied sum."
  call sample_sum_s_par_4(error)
  if (rank == 0) then
    if (allocated(error)) error stop "Test (16/16): FAILED."
  endif
  if (rank == 0) write (error_unit, "(A)") "Test (16/16): PASSED."

  if (rank == 0) write (error_unit, "(A)") "Suite: Randomized Suite:"

  if (rank == 0) write (error_unit, "(A)") &
    "Test (1/2): Testing randomized sampling for random nk, int. and cont. indices. OMP+MPI paralleization."
  call sample_rnd_no_sum(error)
  if (rank == 0) then
    if (allocated(error)) error stop "Test (1/2): FAILED."
  endif
  if (rank == 0) write (error_unit, "(A)") "Test (1/2): PASSED."

  if (rank == 0) write (error_unit, "(A)") &
    "Test (2/2): Testing randomized sampling for random nk, int. and cont. indices. OMP+MPI paralleization. Implied sum."
  call sample_rnd_sum(error)
  if (rank == 0) then
    if (allocated(error)) error stop "Test (2/2): FAILED."
  endif
  if (rank == 0) write (error_unit, "(A)") "Test (2/2): PASSED."

  if (rank == 0) write (error_unit, "(A)") "Suite: Utilities Suite:"

  if (rank == 0) write (error_unit, "(A)") &
    "Test (1/3): Testing kpath utility."
  call test_kpath(error)
  if (rank == 0) then
    if (allocated(error)) error stop "Test (1/3): FAILED."
  endif
  if (rank == 0) write (error_unit, "(A)") "Test (1/3): PASSED."

  if (rank == 0) write (error_unit, "(A)") &
    "Test (2/3): Testing kslice utility."
  call test_kslice(error)
  if (rank == 0) then
    if (allocated(error)) error stop "Test (2/3): FAILED."
  endif
  if (rank == 0) write (error_unit, "(A)") "Test (2/3): PASSED."

  if (rank == 0) write (error_unit, "(A)") &
    "Test (3/3): Testing crys_to cart and cart_to_crys utilities."
  call test_crys_to_cart_cart_to_crys(error)
  if (rank == 0) then
    if (allocated(error)) error stop "Test (3/3): FAILED."
  endif
  if (rank == 0) write (error_unit, "(A)") "Test (3/3): PASSED."

  if (rank == 0) write (error_unit, "(A)") "Suite: Materials Properties Suite:"
  if (rank == 0) write (error_unit, "(A)") &
    "Test (1/1): Testing consistency on the calculation of the optical conductivity in GaAs."
  call sample_and_sum_GaAs_optcond_in_peak(error)
  if (rank == 0) then
    if (allocated(error)) error stop "Test (1/1): FAILED."
  endif
  if (rank == 0) write (error_unit, "(A)") "Test (1/1): PASSED."

  call MPI_FINALIZE(ierror)

end program driver
