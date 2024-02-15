module Randomized_Suite
  use, intrinsic :: iso_fortran_env, only: error_unit
  use MPI
  use OMP_LIB

  use SsTC_driver_kinds, only: wp => dp
  use SsTC_driver_definitions, only: cmplx_1
  use SsTC_driver, only: task_specifier
  use WannInt, only: crystal

  use testdrive, only: error_type

  implicit none
  private

  real(wp) :: tol = 1.0E4_wp

  public :: sample_rnd_no_sum, sample_rnd_sum

contains

  subroutine sample_rnd_no_sum(error)
    type(error_type), allocatable, intent(out) :: error

    type(crystal) :: dummy
    type(task_specifier) :: tsk
    complex(wp), allocatable :: result(:, :, :)

    integer :: nk, n_int_ind, n_cont_ind
    integer, allocatable :: int_ind(:), cont_ind(:)
    real(wp), allocatable :: cont_ind_start(:), cont_ind_end(:)

    real(wp) :: rnk, rn_int_ind, rn_cont_ind
    real(wp), allocatable :: rint_ind(:), rcont_ind(:)

    real(wp), allocatable :: kpath(:, :)
    integer :: i, kpart(3)

    real(wp) :: re, im

    integer :: ierror, rank

    !Params for sampling tests.
    integer, parameter :: minnk = 5000, maxnk = 10000, &
                          minii = 1, maxii = 2, &
                          minci = 1, maxci = 2, &
                          minii_range = 1, minci_range = 3, &
                          maxii_range = 2, maxci_range = 50

    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)

    call random_seed()

    !randomize nk, number of int_ind and cont_ind.
    call random_number(rnk)
    nk = nint(real(minnk, wp) + real(maxnk - minnk, wp)*rnk)
    call random_number(rn_int_ind)
    n_int_ind = nint(real(minii, wp) + real(maxii - minii, wp)*rn_int_ind)
    call random_number(rn_cont_ind)
    n_cont_ind = nint(real(minci, wp) + real(maxci - minci, wp)*rn_cont_ind)

    !DBG: Manually set variables for error reproducibility.
    !nk = 10000
    !n_int_ind = maxii
    !n_cont_ind = maxci

    call MPI_BCAST(nk, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
    call MPI_BCAST(n_int_ind, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
    call MPI_BCAST(n_cont_ind, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)

    allocate (kpath(3, nk))

    !randomize dimension specifiers for int_ind and cont_ind.
    allocate (rint_ind(n_int_ind), int_ind(n_int_ind), &
              rcont_ind(n_cont_ind), cont_ind(n_cont_ind), &
              cont_ind_start(n_cont_ind), cont_ind_end(n_cont_ind))

    cont_ind_start = 1.0_wp
    cont_ind_end = 2.0_wp

    call random_number(rint_ind)
    int_ind = nint(real(minii_range, wp) + real(maxii_range - minii_range, wp)*rint_ind)
    call random_number(rcont_ind)
    cont_ind = nint(real(minci_range, wp) + real(maxci_range - minci_range, wp)*rcont_ind)

    !DBG: Manually set variables for error reproducibility.
    !int_ind = [3, 3]
    !cont_ind = [50, 50]

    call MPI_BCAST(int_ind, size(int_ind), MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
    call MPI_BCAST(cont_ind, size(cont_ind), MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
    call MPI_BCAST(cont_ind_start, size(cont_ind_start), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
    call MPI_BCAST(cont_ind_end, size(cont_ind_end), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)

    if (rank == 0) write (error_unit, "(A)") "Info:"
    if (rank == 0) write (error_unit, "(A, i0, A)") "  Number of k-points = ", nk, "."
    if (rank == 0) write (error_unit, "(A)") "  Integer indices specifier:"
    if (rank == 0) write (error_unit, *) int_ind
    if (rank == 0) write (error_unit, "(A)") "  Continuous indices specifier:"
    if (rank == 0) write (error_unit, *) cont_ind

    do i = 1, nk
      kpath(:, i) = -[0.5_wp, 0.5_wp, 0.5_wp] + [1.0_wp, 1.0_wp, 1.0_wp]*real(i - 1, wp)/real(nk - 1, wp)
    enddo
    kpart = nint(real(nk, wp)**(1.0_wp/3.0_wp))

    call tsk%construct(name="test", &
                       int_ind=int_ind, &
                       cont_data_start=cont_ind_start, &
                       cont_data_end=cont_ind_end, &
                       cont_data_steps=cont_ind, &
                       exponent=10.0_wp, &
                       calculator=test_sampler)

    call tsk%sample(dummy, kpath, result, "MPI+OMP")

    re = real(result(1, 1, 1), wp)
    im = aimag(result(1, 1, 1))

    if ((abs(re - 1.0_wp) > tol*epsilon(1.0_wp)) .or. &
        (abs(im) > tol*epsilon(1.0_wp))) &
      allocate (error)

    call tsk%sample(dummy, kpart, result, "MPI+OMP")

    re = real(result(1, 1, 1), wp)
    im = aimag(result(1, 1, 1))

    if ((abs(re - 1.0_wp) > tol*epsilon(1.0_wp)) .or. &
        (abs(im) > tol*epsilon(1.0_wp))) &
      allocate (error)

    call tsk%construct(name="test", &
                       int_ind=int_ind, &
                       cont_data_start=cont_ind_start, &
                       cont_data_end=cont_ind_end, &
                       cont_data_steps=cont_ind, &
                       calculator=test_sampler)

    call tsk%sample(dummy, kpath, result, "MPI+OMP")

    re = real(result(1, 1, 1), wp)
    im = aimag(result(1, 1, 1))

    if ((abs(re - 1.0_wp) > tol*epsilon(1.0_wp)) .or. &
        (abs(im) > tol*epsilon(1.0_wp))) &
      allocate (error)

    call tsk%sample(dummy, kpart, result, "MPI+OMP")

    re = real(result(1, 1, 1), wp)
    im = aimag(result(1, 1, 1))

    if ((abs(re - 1.0_wp) > tol*epsilon(1.0_wp)) .or. &
        (abs(im) > tol*epsilon(1.0_wp))) &
      allocate (error)

    deallocate (int_ind, cont_ind, cont_ind_start, cont_ind_end, rint_ind, rcont_ind, kpath)

  end subroutine sample_rnd_no_sum

  subroutine sample_rnd_sum(error)
    type(error_type), allocatable, intent(out) :: error

    type(crystal) :: dummy
    type(task_specifier) :: tsk
    complex(wp), allocatable :: result(:, :)

    integer :: nk, n_int_ind, n_cont_ind
    integer, allocatable :: int_ind(:), cont_ind(:)
    real(wp), allocatable :: cont_ind_start(:), cont_ind_end(:)

    real(wp) :: rnk, rn_int_ind, rn_cont_ind
    real(wp), allocatable :: rint_ind(:), rcont_ind(:)

    real(wp), allocatable :: kpath(:, :)
    integer :: i, kpart(3)

    real(wp) :: re, im

    integer :: ierror, rank

    !Params for sampling and summing tests.
    integer, parameter :: minnk = 10000, maxnk = 100000, &
                          minii = 1, maxii = 2, &
                          minci = 1, maxci = 2, &
                          minii_range = 1, minci_range = 3, &
                          maxii_range = 2, maxci_range = 100

    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)

    call random_seed()

    !randomize nk, number of int_ind and cont_ind.
    call random_number(rnk)
    nk = nint(real(minnk, wp) + real(maxnk - minnk, wp)*rnk)
    call random_number(rn_int_ind)
    n_int_ind = nint(real(minii, wp) + real(maxii - minii, wp)*rn_int_ind)
    call random_number(rn_cont_ind)
    n_cont_ind = nint(real(minci, wp) + real(maxci - minci, wp)*rn_cont_ind)

    call MPI_BCAST(nk, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
    call MPI_BCAST(n_int_ind, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
    call MPI_BCAST(n_cont_ind, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)

    !DBG: Manually set variables for error reproducibility.
    !nk = maxnk
    !n_int_ind = maxii
    !n_cont_ind = maxci

    allocate (kpath(3, nk))

    !randomize dimension specifiers for int_ind and cont_ind.
    allocate (rint_ind(n_int_ind), int_ind(n_int_ind), &
              rcont_ind(n_cont_ind), cont_ind(n_cont_ind), &
              cont_ind_start(n_cont_ind), cont_ind_end(n_cont_ind))

    cont_ind_start = 1.0_wp
    cont_ind_end = 2.0_wp

    call random_number(rint_ind)
    int_ind = nint(real(minii_range, wp) + real(maxii_range - minii_range, wp)*rint_ind)
    call random_number(rcont_ind)
    cont_ind = nint(real(minci_range, wp) + real(maxci_range - minci_range, wp)*rcont_ind)

    !DBG: Manually set variables for error reproducibility.
    !int_ind = [3, 3]
    !cont_ind = [100, 100]

    call MPI_BCAST(int_ind, size(int_ind), MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
    call MPI_BCAST(cont_ind, size(cont_ind), MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
    call MPI_BCAST(cont_ind_start, size(cont_ind_start), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
    call MPI_BCAST(cont_ind_end, size(cont_ind_end), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)

    if (rank == 0) write (error_unit, "(A)") "Info:"
    if (rank == 0) write (error_unit, "(A, i0, A)") "  Number of k-points = ", nk, "."
    if (rank == 0) write (error_unit, "(A)") "  Integer indices specifier:"
    if (rank == 0) write (error_unit, *) int_ind
    if (rank == 0) write (error_unit, "(A)") "  Continuous indices specifier:"
    if (rank == 0) write (error_unit, *) cont_ind

    do i = 1, nk
      kpath(:, i) = -[0.5_wp, 0.5_wp, 0.5_wp] + [1.0_wp, 1.0_wp, 1.0_wp]*real(i - 1, wp)/real(nk - 1, wp)
    enddo
    kpart = nint(real(nk, wp)**(1.0_wp/3.0_wp))

    call tsk%construct(name="test", &
                       int_ind=int_ind, &
                       cont_data_start=cont_ind_start, &
                       cont_data_end=cont_ind_end, &
                       cont_data_steps=cont_ind, &
                       exponent=10.0_wp, &
                       calculator=test_sampler)

    call tsk%sample(dummy, kpath, result, "MPI+OMP")

    re = real(result(1, 1), wp)
    im = aimag(result(1, 1))

    if ((abs(re - real(nk, wp)) > tol*epsilon(1.0_wp)) .or. &
        (abs(im) > tol*epsilon(1.0_wp))) &
      allocate (error)

    call tsk%sample(dummy, kpart, result, "MPI+OMP")

    re = real(result(1, 1), wp)
    im = aimag(result(1, 1))

    if ((abs(re - real(nint(real(nk, wp)**(1.0_wp/3.0_wp)), wp)**3) > tol*epsilon(1.0_wp)) .or. &
        (abs(im) > tol*epsilon(1.0_wp))) &
      allocate (error)

    call tsk%construct(name="test", &
                       int_ind=int_ind, &
                       cont_data_start=cont_ind_start, &
                       cont_data_end=cont_ind_end, &
                       cont_data_steps=cont_ind, &
                       calculator=test_sampler)

    call tsk%sample(dummy, kpath, result, "MPI+OMP")

    re = real(result(1, 1), wp)
    im = aimag(result(1, 1))

    if ((abs(re - real(nk, wp)) > tol*epsilon(1.0_wp)) .or. &
        (abs(im) > tol*epsilon(1.0_wp))) &
      allocate (error)

    call tsk%sample(dummy, kpart, result, "MPI+OMP")

    re = real(result(1, 1), wp)
    im = aimag(result(1, 1))

    if ((abs(re - real(nint(real(nk, wp)**(1.0_wp/3.0_wp)), wp)**3) > tol*epsilon(1.0_wp)) .or. &
        (abs(im) > tol*epsilon(1.0_wp))) &
      allocate (error)

    deallocate (int_ind, cont_ind, cont_ind_start, cont_ind_end, rint_ind, rcont_ind, kpath)

  end subroutine sample_rnd_sum

  function test_sampler(self, crys, k, other)
    class(task_specifier), intent(in) :: self
    class(crystal), intent(in) :: crys
    real(wp), intent(in) :: k(3)
    class(*), optional, intent(in) :: other

    logical :: dbg1, dbg2, dbg3

    complex(wp) :: test_sampler(self%idims%size(), self%cdims%size())

    dbg1 = crys%initialized()
    dbg2 = present(other)
    dbg3 = (nint(k(1)) == 0)

    test_sampler = cmplx_1

  end function test_sampler

end module Randomized_Suite
