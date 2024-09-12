#include "macros.inc"
module SsTC_driver

  use MPI
  use OMP_LIB

  use SsTC_driver_kinds, only: wp => dp
  use SsTC_driver_definitions, only: cmplx_0
  use SsTC_driver_comms, only: is_mpi_initialized, is_mpi_finalized, &
    rank, nProcs, ierror, &
    get_MPI_task_partition
  use SsTC_driver_utils, only: kpath, kslice, crys_to_cart, cart_to_crys, kpath_length
  use MAC, only: container_specifier, container
  use WannInt, only: crystal, diagonalize, inverse, &
    dirac_delta, &
    deg_list, schur, &
    SVD, expsh, logu

  implicit none
  private

  public :: kpath
  public :: kslice
  public :: crys_to_cart
  public :: cart_to_crys
  public :: kpath_length
  public :: container_specifier, container
  public :: crystal
  public :: diagonalize
  public :: inverse
  public :: dirac_delta
  public :: deg_list
  public :: schur
  public :: SVD
  public :: expsh
  public :: logu

  type, public :: task_specifier
    private
    character(len=120) :: nm
    type(container_specifier), public :: idims
    type(container_specifier), public :: cdims
    procedure(cs), pointer, pass(self) :: calc => null()
    type(container), allocatable :: cont_dt(:)
    logical :: task_initialized = .false.
  contains
    private
    procedure, public, pass(self) :: construct => cons
    procedure, pass(self) :: sample_pre_no_sum
    procedure, pass(self) :: sample_part_no_sum
    procedure, pass(self) :: sample_pre_sum
    procedure, pass(self) :: sample_part_sum
    generic, public :: sample => sample_pre_no_sum, sample_part_no_sum, &
      sample_pre_sum, sample_part_sum
    procedure, public, pass(self) :: name
    procedure, public, pass(self) :: cdt
    procedure, public, pass(self) :: initialized => is_initialized
  end type

  abstract interface
    function cs(self, crys, k, other)
      import wp, task_specifier, crystal
      class(task_specifier), intent(in) :: self
      class(crystal), intent(in) :: crys
      real(wp), intent(in) :: k(3)
      class(*), optional, intent(in) :: other

      complex(wp) :: cs(self%idims%size(), self%cdims%size())
    end function cs
  end interface

contains

  subroutine cons(self, name, &
                  int_ind, &
                  cont_data_start, cont_data_end, cont_data_steps, exponent, &
                  calculator)

    class(task_specifier), intent(out) :: self
    character(len=*), intent(in) :: name

    integer, optional, intent(in) :: int_ind(:)
    real(wp), optional, intent(in) :: cont_data_start(:), cont_data_end(:)
    integer, optional, intent(in) :: cont_data_steps(:)
    real(wp), optional, intent(in) :: exponent

    procedure(cs) :: calculator

    integer :: i, j
    real(wp) :: discretization

    character(len=1024) :: errormsg
    integer :: istat

    self%nm = name

    if (present(int_ind)) then
      call self%idims%specify(dimension_specifier=int_ind)
    else
      call self%idims%specify(dimension_specifier=[1])
    endif

    if (present(cont_data_start) .and. &
        present(cont_data_end) .and. &
        present(cont_data_steps)) then

      do i = 1, size(cont_data_steps)
        if (cont_data_steps(i) < 1) error stop &
          "SsTC_driver: Error #1: the elements of cont_data_steps must be positive integers."
      enddo

      call self%cdims%specify(dimension_specifier=cont_data_steps)

      if ((size(cont_data_start) /= self%cdims%rank()) .or. &
          (size(cont_data_end) /= self%cdims%rank())) then

        error stop "SsTC_driver: Error #1: size of cont_data_start, cont_data_end must match size of cont_data_steps."

      else

        allocate (self%cont_dt(self%cdims%rank()), stat=istat)
        if (istat /= 0) then
          write (errormsg, "(i20)") istat
          errormsg = "SsTC_driver: Error #2: failure allocating cont_dt. stat = "//trim(adjustl(errormsg))//"."
          error stop trim(errormsg)
        endif

        if (present(exponent)) then

          do i = 1, self%cdims%rank() !For each dimension,
            !initialize 1D container holding data,
            call self%cont_dt(i)%construct(container_type="real_dp", dimension_specifier=[cont_data_steps(i)])
            do j = 1, cont_data_steps(i) !and for each data step,
              !compute data discretization,
              if (cont_data_steps(i) == 1) then
                discretization = cont_data_start(i)
              else
                discretization = cont_data_start(i) + &
                                 (cont_data_end(i) - cont_data_start(i))*real(j - 1, wp)/real(cont_data_steps(i) - 1, wp)
              endif
              call self%cont_dt(i)%set(val=exponent**discretization, at=j) !and store.
            enddo
          enddo

        else

          do i = 1, self%cdims%rank() !For each dimension,
            !initialize 1D container holding data,
            call self%cont_dt(i)%construct(container_type="real_dp", dimension_specifier=[cont_data_steps(i)])
            do j = 1, cont_data_steps(i) !and for each data step,
              !compute data discretization,
              if (cont_data_steps(i) == 1) then
                discretization = cont_data_start(i)
              else
                discretization = cont_data_start(i) + &
                                 (cont_data_end(i) - cont_data_start(i))*real(j - 1, wp)/real(cont_data_steps(i) - 1, wp)
              endif
              call self%cont_dt(i)%set(val=discretization, at=j) !and store.
            enddo
          enddo

        endif

      endif

    else

      call self%cdims%specify(dimension_specifier=[1])

      allocate (self%cont_dt(1), stat=istat)
      if (istat /= 0) then
        write (errormsg, "(i20)") istat
        errormsg = "SsTC_driver: Error #2: failure allocating cont_dt. stat = "//trim(adjustl(errormsg))//"."
        error stop trim(errormsg)
      endif

      call self%cont_dt(1)%construct(container_type="real_dp", dimension_specifier=[1])
      self%cont_dt(1)%rdp_storage = 1.0_wp

    endif

    self%calc => calculator

    self%task_initialized = .true.

  end subroutine cons

  subroutine sample_pre_no_sum(self, crys, klist, store_at, parallelization, other)
    class(task_specifier), intent(in) :: self
    class(crystal), intent(in) :: crys
    real(wp), intent(in) :: klist(:, :)
    complex(wp), allocatable, intent(out) :: store_at(:, :, :)
    character(len=*), optional, intent(in) :: parallelization
    class(*), optional, intent(in) :: other

    logical :: using_mpi = .false., &
               using_omp = .true.
    logical :: is_mpi_initialized, is_mpi_finalized
    integer, allocatable :: counts(:), displs(:)

    complex(wp), allocatable :: local_data_k(:, :, :)

    integer :: i, r

    integer :: nk, ik
    real(wp) :: k(3)

    character(len=1024) :: errormsg
    integer :: istat

    if (.not. (self%task_initialized)) error stop &
      "SsTC_driver: Error #4: task_specifier is not initialized."

    nk = size(klist(1, :))
    if (size(klist(:, 1)) /= 3) error stop "SsTC_driver: Error #1: size of klist(:, 1) is not 3."

    if (present(parallelization)) then
      select case (parallelization)
      case ("MPI+OMP", "OMP+MPI")
        using_mpi = .true.
        using_omp = .true.
      case ("MPI")
        using_mpi = .true.
        using_omp = .false.
      case ("OMP")
        using_mpi = .false.
        using_omp = .true.
      case ("none")
        using_mpi = .false.
        using_omp = .false.
      case default
        error stop "SsTC_driver: Error #1: requested parallelization not recognized."
      end select
    endif

    allocate (store_at(nk, self%idims%size(), self%cdims%size()), stat=istat)
    if (istat /= 0) then
      write (errormsg, "(i20)") istat
      errormsg = "SsTC_driver: Error #2: failure allocating store_at. stat = "//trim(adjustl(errormsg))//"."
      error stop trim(errormsg)
    endif
    store_at = cmplx_0

    if (using_mpi) then

      call MPI_INITIALIZED(is_mpi_initialized, ierror)
      call MPI_FINALIZED(is_mpi_finalized, ierror)
      if (ierror /= 0) error stop "SsTC_driver: Error #3: MPI error when running MPI_INITIALIZED, MPI_FINALIZED."

      if (.not. ((is_mpi_initialized) .and. ((.not. is_mpi_finalized)))) then
        error stop "SsTC_driver: Error #3: MPI has not been initialized or has been finalized."
      endif

      call MPI_COMM_SIZE(MPI_COMM_WORLD, nProcs, ierror)
      call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)
      if (ierror /= 0) error stop "SsTC_driver: Error #3: MPI error when running MPI_COMM_SIZE, MPI_COMM_RANK."

      allocate (counts(0:nProcs - 1), displs(0:nProcs - 1), stat=istat)
      if (istat /= 0) then
        write (errormsg, "(i20)") istat
        errormsg = "SsTC_driver: Error #2: failure allocating counts, displs. stat = "//trim(adjustl(errormsg))//"."
        error stop trim(errormsg)
      endif
      call get_MPI_task_partition(nk, nProcs, counts, displs)

      allocate (local_data_k(displs(rank) + 1:displs(rank) + counts(rank), self%idims%size(), self%cdims%size()), &
                stat=istat)
      if (istat /= 0) then
        write (errormsg, "(i20)") istat
        errormsg = "SsTC_driver: Error #2: failure allocating local_data_k. stat = "//trim(adjustl(errormsg))//"."
        error stop trim(errormsg)
      endif
      local_data_k = cmplx_0

      if (using_omp) then !MPI+OMP.

        !_OMPTGT_(PARALLEL DO &)
        !_OMPTGT_(SHARED(self, crys, displs, counts, rank, local_data_k) &)
        !_OMPTGT_(PRIVATE(ik, k))
        do ik = displs(rank) + 1, displs(rank) + counts(rank)
          k = klist(:, ik)
          if (present(other)) then
            local_data_k(ik, :, :) = self%calc(crys, k, other)
          else
            local_data_k(ik, :, :) = self%calc(crys, k)
          endif
        enddo

      else !MPI

        do ik = displs(rank) + 1, displs(rank) + counts(rank)
          k = klist(:, ik)
          if (present(other)) then
            local_data_k(ik, :, :) = self%calc(crys, k, other)
          else
            local_data_k(ik, :, :) = self%calc(crys, k)
          endif
        enddo

      endif

      do r = 1, self%cdims%size() !For each integer index.
        do i = 1, self%idims%size() !For each continuous index.
          call MPI_ALLGATHERV(local_data_k(:, i, r), &
                              size(local_data_k(:, i, r)), &
                              MPI_COMPLEX16, &
                              store_at(:, i, r), &
                              counts, &
                              displs, &
                              MPI_COMPLEX16, &
                              MPI_COMM_WORLD, &
                              ierror)
          if (ierror /= 0) error stop "SsTC_driver: Error #3: MPI error when running MPI_ALLGATHERV."
        enddo
      enddo

    else

      if (using_omp) then !OMP.

        !_OMPTGT_(PARALLEL DO &)
        !_OMPTGT_(SHARED(self, crys, displs, counts, rank, store_at) &)
        !_OMPTGT_(PRIVATE(ik, k))
        do ik = 1, nk
          k = klist(:, ik)
          if (present(other)) then
            store_at(ik, :, :) = self%calc(crys, k, other)
          else
            store_at(ik, :, :) = self%calc(crys, k)
          endif
        enddo

      else !No parallelization.

        do ik = 1, nk
          k = klist(:, ik)
          if (present(other)) then
            store_at(ik, :, :) = self%calc(crys, k, other)
          else
            store_at(ik, :, :) = self%calc(crys, k)
          endif
        enddo

      endif

    endif

  end subroutine sample_pre_no_sum

  subroutine sample_part_no_sum(self, crys, kpart, store_at, parallelization, other)
    class(task_specifier), intent(in) :: self
    class(crystal), intent(in) :: crys
    integer, intent(in) :: kpart(3)
    complex(wp), allocatable, intent(out) :: store_at(:, :, :)
    character(len=*), optional, intent(in) :: parallelization
    class(*), optional, intent(in) :: other

    logical :: using_mpi = .false., &
               using_omp = .true.
    logical :: is_mpi_initialized, is_mpi_finalized
    integer, allocatable :: counts(:), displs(:)

    complex(wp), allocatable :: local_data_k(:, :, :)

    integer :: i, r, l

    integer :: nk, ik, k_arr(3)
    real(wp) :: k(3)
    type(container_specifier) :: k_handle

    character(len=1024) :: errormsg
    integer :: istat

    if (.not. (self%task_initialized)) error stop &
      "SsTC_driver: Error #4: task_specifier is not initialized."

    do l = 1, 3
      if (kpart(l) < 1) error stop "SsTC_driver: Error #1: the elements of kpart must be positive integers."
    enddo
    nk = product(kpart)

    if (present(parallelization)) then
      select case (parallelization)
      case ("MPI+OMP", "OMP+MPI")
        using_mpi = .true.
        using_omp = .true.
      case ("MPI")
        using_mpi = .true.
        using_omp = .false.
      case ("OMP")
        using_mpi = .false.
        using_omp = .true.
      case ("none")
        using_mpi = .false.
        using_omp = .false.
      case default
        error stop "SsTC_driver: Error #1: requested parallelization not recognized."
      end select
    endif

    allocate (store_at(nk, self%idims%size(), self%cdims%size()), stat=istat)
    if (istat /= 0) then
      write (errormsg, "(i20)") istat
      errormsg = "SsTC_driver: Error #2: failure allocating store_at. stat = "//trim(adjustl(errormsg))//"."
      error stop trim(errormsg)
    endif
    store_at = cmplx_0

    if (using_mpi) then

      call MPI_INITIALIZED(is_mpi_initialized, ierror)
      call MPI_FINALIZED(is_mpi_finalized, ierror)
      if (ierror /= 0) error stop "SsTC_driver: Error #3: MPI error when running MPI_INITIALIZED, MPI_FINALIZED."

      if (.not. ((is_mpi_initialized) .and. ((.not. is_mpi_finalized)))) then
        error stop "SsTC_driver: Error #3: MPI has not been initialized or has been finalized."
      endif

      call MPI_COMM_SIZE(MPI_COMM_WORLD, nProcs, ierror)
      call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)
      if (ierror /= 0) error stop "SsTC_driver: Error #3: MPI error when running MPI_COMM_SIZE, MPI_COMM_RANK."

      allocate (counts(0:nProcs - 1), displs(0:nProcs - 1), stat=istat)
      if (istat /= 0) then
        write (errormsg, "(i20)") istat
        errormsg = "SsTC_driver: Error #2: failure allocating counts, displs. stat = "//trim(adjustl(errormsg))//"."
        error stop trim(errormsg)
      endif
      call get_MPI_task_partition(nk, nProcs, counts, displs)

      allocate (local_data_k(displs(rank) + 1:displs(rank) + counts(rank), self%idims%size(), self%cdims%size()), &
                stat=istat)
      if (istat /= 0) then
        write (errormsg, "(i20)") istat
        errormsg = "SsTC_driver: Error #2: failure allocating local_data_k. stat = "//trim(adjustl(errormsg))//"."
        error stop trim(errormsg)
      endif
      local_data_k = cmplx_0

      if (using_omp) then !MPI+OMP.

        call k_handle%specify(kpart)
        !_OMPTGT_(PARALLEL DO &)
        !_OMPTGT_(SHARED(self, crys, displs, counts, rank, local_data_k, k_handle) &)
        !_OMPTGT_(PRIVATE(ik, k_arr, k))
        do ik = displs(rank) + 1, displs(rank) + counts(rank)
          k_arr = k_handle%ind(ik)
          if (kpart(1) == 1) then
            k(1) = 0.0_wp
          else
            k(1) = -0.5_wp + real(k_arr(1) - 1, wp)/real(kpart(1) - 1, wp)
          endif
          if (kpart(2) == 1) then
            k(2) = 0.0_wp
          else
            k(2) = -0.5_wp + real(k_arr(2) - 1, wp)/real(kpart(2) - 1, wp)
          endif
          if (kpart(3) == 1) then
            k(3) = 0.0_wp
          else
            k(3) = -0.5_wp + real(k_arr(3) - 1, wp)/real(kpart(3) - 1, wp)
          endif
          if (present(other)) then
            local_data_k(ik, :, :) = self%calc(crys, k, other)
          else
            local_data_k(ik, :, :) = self%calc(crys, k)
          endif
        enddo

      else !MPI

        call k_handle%specify(kpart)
        do ik = displs(rank) + 1, displs(rank) + counts(rank)
          k_arr = k_handle%ind(ik)
          if (kpart(1) == 1) then
            k(1) = 0.0_wp
          else
            k(1) = -0.5_wp + real(k_arr(1) - 1, wp)/real(kpart(1) - 1, wp)
          endif
          if (kpart(2) == 1) then
            k(2) = 0.0_wp
          else
            k(2) = -0.5_wp + real(k_arr(2) - 1, wp)/real(kpart(2) - 1, wp)
          endif
          if (kpart(3) == 1) then
            k(3) = 0.0_wp
          else
            k(3) = -0.5_wp + real(k_arr(3) - 1, wp)/real(kpart(3) - 1, wp)
          endif
          if (present(other)) then
            local_data_k(ik, :, :) = self%calc(crys, k, other)
          else
            local_data_k(ik, :, :) = self%calc(crys, k)
          endif
        enddo
      endif

      do r = 1, self%cdims%size() !For each integer index.
        do i = 1, self%idims%size() !For each continuous index.
          call MPI_ALLGATHERV(local_data_k(:, i, r), &
                              size(local_data_k(:, i, r)), &
                              MPI_COMPLEX16, &
                              store_at(:, i, r), &
                              counts, &
                              displs, &
                              MPI_COMPLEX16, &
                              MPI_COMM_WORLD, &
                              ierror)
          if (ierror /= 0) error stop "SsTC_driver: Error #3: MPI error when running MPI_ALLGATHERV."
        enddo
      enddo

    else

      if (using_omp) then !OMP.

        call k_handle%specify(kpart)
        !_OMPTGT_(PARALLEL DO &)
        !_OMPTGT_(SHARED(self, crys, displs, counts, rank, store_at, k_handle) &)
        !_OMPTGT_(PRIVATE(ik, k_arr, k))
        do ik = 1, nk
          k_arr = k_handle%ind(ik)
          if (kpart(1) == 1) then
            k(1) = 0.0_wp
          else
            k(1) = -0.5_wp + real(k_arr(1) - 1, wp)/real(kpart(1) - 1, wp)
          endif
          if (kpart(2) == 1) then
            k(2) = 0.0_wp
          else
            k(2) = -0.5_wp + real(k_arr(2) - 1, wp)/real(kpart(2) - 1, wp)
          endif
          if (kpart(3) == 1) then
            k(3) = 0.0_wp
          else
            k(3) = -0.5_wp + real(k_arr(3) - 1, wp)/real(kpart(3) - 1, wp)
          endif
          if (present(other)) then
            store_at(ik, :, :) = self%calc(crys, k, other)
          else
            store_at(ik, :, :) = self%calc(crys, k)
          endif
        enddo

      else !No parallelization.

        call k_handle%specify(kpart)
        do ik = 1, nk
          k_arr = k_handle%ind(ik)
          if (kpart(1) == 1) then
            k(1) = 0.0_wp
          else
            k(1) = -0.5_wp + real(k_arr(1) - 1, wp)/real(kpart(1) - 1, wp)
          endif
          if (kpart(2) == 1) then
            k(2) = 0.0_wp
          else
            k(2) = -0.5_wp + real(k_arr(2) - 1, wp)/real(kpart(2) - 1, wp)
          endif
          if (kpart(3) == 1) then
            k(3) = 0.0_wp
          else
            k(3) = -0.5_wp + real(k_arr(3) - 1, wp)/real(kpart(3) - 1, wp)
          endif
          if (present(other)) then
            store_at(ik, :, :) = self%calc(crys, k, other)
          else
            store_at(ik, :, :) = self%calc(crys, k)
          endif
        enddo

      endif

    endif

  end subroutine sample_part_no_sum

  subroutine sample_pre_sum(self, crys, klist, store_at, parallelization, other)
    class(task_specifier), intent(in) :: self
    class(crystal), intent(in) :: crys
    real(wp), intent(in) :: klist(:, :)
    complex(wp), allocatable, intent(out) :: store_at(:, :)
    character(len=*), optional, intent(in) :: parallelization
    class(*), optional, intent(in) :: other

    logical :: using_mpi = .false., &
               using_omp = .true.
    logical :: is_mpi_initialized, is_mpi_finalized
    integer, allocatable :: counts(:), displs(:)

    complex(wp), allocatable :: local_data_k(:, :)

    integer :: nk, ik
    real(wp) :: k(3)

    character(len=1024) :: errormsg
    integer :: istat

    if (.not. (self%task_initialized)) error stop &
      "SsTC_driver: Error #4: task_specifier is not initialized."

    nk = size(klist(1, :))
    if (size(klist(:, 1)) /= 3) error stop "SsTC_driver: Error #1: size of klist(:, 1) is not 3."

    if (present(parallelization)) then
      select case (parallelization)
      case ("MPI+OMP", "OMP+MPI")
        using_mpi = .true.
        using_omp = .true.
      case ("MPI")
        using_mpi = .true.
        using_omp = .false.
      case ("OMP")
        using_mpi = .false.
        using_omp = .true.
      case ("none")
        using_mpi = .false.
        using_omp = .false.
      case default
        error stop "SsTC_driver: Error #1: requested parallelization not recognized."
      end select
    endif

    allocate (store_at(self%idims%size(), self%cdims%size()), stat=istat)
    if (istat /= 0) then
      write (errormsg, "(i20)") istat
      errormsg = "SsTC_driver: Error #2: failure allocating store_at. stat = "//trim(adjustl(errormsg))//"."
      error stop trim(errormsg)
    endif
    store_at = cmplx_0

    if (using_mpi) then

      call MPI_INITIALIZED(is_mpi_initialized, ierror)
      call MPI_FINALIZED(is_mpi_finalized, ierror)
      if (ierror /= 0) error stop "SsTC_driver: Error #3: MPI error when running MPI_INITIALIZED, MPI_FINALIZED."

      if (.not. ((is_mpi_initialized) .and. ((.not. is_mpi_finalized)))) then
        error stop "SsTC_driver: Error #3: MPI has not been initialized or has been finalized."
      endif

      call MPI_COMM_SIZE(MPI_COMM_WORLD, nProcs, ierror)
      call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)
      if (ierror /= 0) error stop "SsTC_driver: Error #3: MPI error when running MPI_COMM_SIZE, MPI_COMM_RANK."

      allocate (counts(0:nProcs - 1), displs(0:nProcs - 1), stat=istat)
      if (istat /= 0) then
        write (errormsg, "(i20)") istat
        errormsg = "SsTC_driver: Error #2: failure allocating counts, displs. stat = "//trim(adjustl(errormsg))//"."
        error stop trim(errormsg)
      endif
      call get_MPI_task_partition(nk, nProcs, counts, displs)

      allocate (local_data_k(self%idims%size(), self%cdims%size()), &
                stat=istat)
      if (istat /= 0) then
        write (errormsg, "(i20)") istat
        errormsg = "SsTC_driver: Error #2: failure allocating local_data_k. stat = "//trim(adjustl(errormsg))//"."
        error stop trim(errormsg)
      endif
      local_data_k = cmplx_0

      if (using_omp) then !MPI+OMP.

        !_OMPTGT_(PARALLEL DO REDUCTION (+: local_data_k) &)
        !_OMPTGT_(SHARED(self, crys, displs, counts, rank) &)
        !_OMPTGT_(PRIVATE(ik, k))
        do ik = displs(rank) + 1, displs(rank) + counts(rank)
          k = klist(:, ik)
          if (present(other)) then
            local_data_k = local_data_k + self%calc(crys, k, other)
          else
            local_data_k = local_data_k + self%calc(crys, k)
          endif
        enddo

      else !MPI

        do ik = displs(rank) + 1, displs(rank) + counts(rank)
          k = klist(:, ik)
          if (present(other)) then
            local_data_k = local_data_k + self%calc(crys, k, other)
          else
            local_data_k = local_data_k + self%calc(crys, k)
          endif
        enddo

      endif

      call MPI_ALLREDUCE(local_data_k, store_at, size(local_data_k), MPI_COMPLEX16, MPI_SUM, MPI_COMM_WORLD, ierror)
      if (ierror /= 0) error stop "SsTC_driver: Error #3: MPI error when running MPI_ALLREDUCE."

    else

      if (using_omp) then !OMP.

        !_OMPTGT_(PARALLEL DO REDUCTION (+: store_at) &)
        !_OMPTGT_(SHARED(self, crys, displs, counts, rank) &)
        !_OMPTGT_(PRIVATE(ik, k))
        do ik = 1, nk
          k = klist(:, ik)
          if (present(other)) then
            store_at = store_at + self%calc(crys, k, other)
          else
            store_at = store_at + self%calc(crys, k)
          endif
        enddo

      else !No parallelization.

        do ik = 1, nk
          k = klist(:, ik)
          if (present(other)) then
            store_at = store_at + self%calc(crys, k, other)
          else
            store_at = store_at + self%calc(crys, k)
          endif
        enddo

      endif

    endif

  end subroutine sample_pre_sum

  subroutine sample_part_sum(self, crys, kpart, store_at, parallelization, other)
    class(task_specifier), intent(in) :: self
    class(crystal), intent(in) :: crys
    integer, intent(in) :: kpart(3)
    complex(wp), allocatable, intent(out) :: store_at(:, :)
    character(len=*), optional, intent(in) :: parallelization
    class(*), optional, intent(in) :: other

    logical :: using_mpi = .false., &
               using_omp = .true.
    logical :: is_mpi_initialized, is_mpi_finalized
    integer, allocatable :: counts(:), displs(:)

    complex(wp), allocatable :: local_data_k(:, :)

    integer :: l

    integer :: nk, ik, k_arr(3)
    real(wp) :: k(3)
    type(container_specifier) :: k_handle

    character(len=1024) :: errormsg
    integer :: istat

    if (.not. (self%task_initialized)) error stop &
      "SsTC_driver: Error #4: task_specifier is not initialized."

    do l = 1, 3
      if (kpart(l) < 1) error stop "SsTC_driver: Error #1: the elements of kpart must be positive integers."
    enddo
    nk = product(kpart)

    if (present(parallelization)) then
      select case (parallelization)
      case ("MPI+OMP", "OMP+MPI")
        using_mpi = .true.
        using_omp = .true.
      case ("MPI")
        using_mpi = .true.
        using_omp = .false.
      case ("OMP")
        using_mpi = .false.
        using_omp = .true.
      case ("none")
        using_mpi = .false.
        using_omp = .false.
      case default
        error stop "SsTC_driver: Error #1: requested parallelization not recognized."
      end select
    endif

    allocate (store_at(self%idims%size(), self%cdims%size()), stat=istat)
    if (istat /= 0) then
      write (errormsg, "(i20)") istat
      errormsg = "SsTC_driver: Error #2: failure allocating store_at. stat = "//trim(adjustl(errormsg))//"."
      error stop trim(errormsg)
    endif
    store_at = cmplx_0

    if (using_mpi) then

      call MPI_INITIALIZED(is_mpi_initialized, ierror)
      call MPI_FINALIZED(is_mpi_finalized, ierror)
      if (ierror /= 0) error stop "SsTC_driver: Error #3: MPI error when running MPI_INITIALIZED, MPI_FINALIZED."

      if (.not. ((is_mpi_initialized) .and. ((.not. is_mpi_finalized)))) then
        error stop "SsTC_driver: Error #3: MPI has not been initialized or has been finalized."
      endif

      call MPI_COMM_SIZE(MPI_COMM_WORLD, nProcs, ierror)
      call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)
      if (ierror /= 0) error stop "SsTC_driver: Error #3: MPI error when running MPI_COMM_SIZE, MPI_COMM_RANK."

      allocate (counts(0:nProcs - 1), displs(0:nProcs - 1), stat=istat)
      if (istat /= 0) then
        write (errormsg, "(i20)") istat
        errormsg = "SsTC_driver: Error #2: failure allocating counts, displs. stat = "//trim(adjustl(errormsg))//"."
        error stop trim(errormsg)
      endif
      call get_MPI_task_partition(nk, nProcs, counts, displs)

      allocate (local_data_k(self%idims%size(), self%cdims%size()), &
                stat=istat)
      if (istat /= 0) then
        write (errormsg, "(i20)") istat
        errormsg = "SsTC_driver: Error #2: failure allocating local_data_k. stat = "//trim(adjustl(errormsg))//"."
        error stop trim(errormsg)
      endif
      local_data_k = cmplx_0

      if (using_omp) then !MPI+OMP.

        call k_handle%specify(kpart)
        !_OMPTGT_(PARALLEL DO REDUCTION (+: local_data_k) &)
        !_OMPTGT_(SHARED(self, crys, displs, counts, rank, k_handle) &)
        !_OMPTGT_(PRIVATE(ik, k_arr, k))
        do ik = displs(rank) + 1, displs(rank) + counts(rank)
          k_arr = k_handle%ind(ik)
          if (kpart(1) == 1) then
            k(1) = 0.0_wp
          else
            k(1) = -0.5_wp + real(k_arr(1) - 1, wp)/real(kpart(1) - 1, wp)
          endif
          if (kpart(2) == 1) then
            k(2) = 0.0_wp
          else
            k(2) = -0.5_wp + real(k_arr(2) - 1, wp)/real(kpart(2) - 1, wp)
          endif
          if (kpart(3) == 1) then
            k(3) = 0.0_wp
          else
            k(3) = -0.5_wp + real(k_arr(3) - 1, wp)/real(kpart(3) - 1, wp)
          endif
          if (present(other)) then
            local_data_k = local_data_k + self%calc(crys, k, other)
          else
            local_data_k = local_data_k + self%calc(crys, k)
          endif
        enddo

      else !MPI

        call k_handle%specify(kpart)
        do ik = displs(rank) + 1, displs(rank) + counts(rank)
          k_arr = k_handle%ind(ik)
          if (kpart(1) == 1) then
            k(1) = 0.0_wp
          else
            k(1) = -0.5_wp + real(k_arr(1) - 1, wp)/real(kpart(1) - 1, wp)
          endif
          if (kpart(2) == 1) then
            k(2) = 0.0_wp
          else
            k(2) = -0.5_wp + real(k_arr(2) - 1, wp)/real(kpart(2) - 1, wp)
          endif
          if (kpart(3) == 1) then
            k(3) = 0.0_wp
          else
            k(3) = -0.5_wp + real(k_arr(3) - 1, wp)/real(kpart(3) - 1, wp)
          endif
          if (present(other)) then
            local_data_k = local_data_k + self%calc(crys, k, other)
          else
            local_data_k = local_data_k + self%calc(crys, k)
          endif
        enddo
      endif

      call MPI_ALLREDUCE(local_data_k, store_at, size(local_data_k), MPI_COMPLEX16, MPI_SUM, MPI_COMM_WORLD, ierror)
      if (ierror /= 0) error stop "SsTC_driver: Error #3: MPI error when running MPI_ALLREDUCE."

    else

      if (using_omp) then !OMP.

        call k_handle%specify(kpart)
        !_OMPTGT_(PARALLEL DO REDUCTION (+: store_at) &)
        !_OMPTGT_(SHARED(self, crys, displs, counts, rank, k_handle) &)
        !_OMPTGT_(PRIVATE(ik, k_arr, k))
        do ik = 1, nk
          k_arr = k_handle%ind(ik)
          if (kpart(1) == 1) then
            k(1) = 0.0_wp
          else
            k(1) = -0.5_wp + real(k_arr(1) - 1, wp)/real(kpart(1) - 1, wp)
          endif
          if (kpart(2) == 1) then
            k(2) = 0.0_wp
          else
            k(2) = -0.5_wp + real(k_arr(2) - 1, wp)/real(kpart(2) - 1, wp)
          endif
          if (kpart(3) == 1) then
            k(3) = 0.0_wp
          else
            k(3) = -0.5_wp + real(k_arr(3) - 1, wp)/real(kpart(3) - 1, wp)
          endif
          if (present(other)) then
            store_at = store_at + self%calc(crys, k, other)
          else
            store_at = store_at + self%calc(crys, k)
          endif
        enddo

      else !No parallelization.

        call k_handle%specify(kpart)
        do ik = 1, nk
          k_arr = k_handle%ind(ik)
          if (kpart(1) == 1) then
            k(1) = 0.0_wp
          else
            k(1) = -0.5_wp + real(k_arr(1) - 1, wp)/real(kpart(1) - 1, wp)
          endif
          if (kpart(2) == 1) then
            k(2) = 0.0_wp
          else
            k(2) = -0.5_wp + real(k_arr(2) - 1, wp)/real(kpart(2) - 1, wp)
          endif
          if (kpart(3) == 1) then
            k(3) = 0.0_wp
          else
            k(3) = -0.5_wp + real(k_arr(3) - 1, wp)/real(kpart(3) - 1, wp)
          endif
          if (present(other)) then
            store_at = store_at + self%calc(crys, k, other)
          else
            store_at = store_at + self%calc(crys, k)
          endif
        enddo

      endif

    endif

  end subroutine sample_part_sum

  pure elemental function name(self)
    class(task_specifier), intent(in) :: self
    character(len=120) :: name
    if (.not. (self%task_initialized)) error stop &
      "SsTC_driver: Error #4: task_specifier is not initialized."
    name = trim(self%nm)
  end function name

  pure elemental function cdt(self, var, step)
    class(task_specifier), intent(in) :: self
    integer, intent(in) :: var, step
    real(wp) :: cdt
    if (.not. (self%task_initialized)) error stop &
      "SsTC_driver: Error #4: task_specifier is not initialized."
    cdt = self%cont_dt(var)%rdp_storage(step)
  end function cdt

  pure elemental logical function is_initialized(self)
    class(task_specifier), intent(in) :: self
    is_initialized = self%task_initialized
  end function is_initialized

end module SsTC_driver
