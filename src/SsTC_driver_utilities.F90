module SsTC_driver_utilities

  use SsTC_driver_kinds, only: wp => dp
  use MAC, only: container_specifier

  implicit none

  private

  public :: kpath
  public :: kslice

contains

  function kpath(vecs, nkpts)
    real(wp), intent(in) :: vecs(:, :)
    integer, intent(in)  :: nkpts(:)

    real(wp), allocatable :: kpath(:, :)

    integer :: nvec, nk
    integer :: i, ivec, isampling, &
               count, icount

    character(len=1024) :: errormsg
    integer :: istat

    nvec = size(vecs(:, 1))
    if (nvec < 2) error stop "SsTC_driver: Error #1: size of vecs(:, 1) must be greater than 1."
    if (size(vecs(1, :)) /= 3) error stop "SsTC_driver: Error #1: size of vecs(1, :) must be 3."
    if (size(nkpts) /= nvec - 1) error stop "SsTC_driver: Error #1: size of nkpts must be equal to size of vecs(:, 1) - 1."

    do i = 1, nvec - 1
      if (nkpts(i) < 2) &
        error stop "SsTC_driver: Error #1: components of nkpts must be positive integers greater than 1."
    enddo

    nk = sum(nkpts) - (nvec - 2)
    allocate (kpath(3, nk), stat=istat)
    if (istat /= 0) then
      write (errormsg, "(i20)") istat
      errormsg = "SsTC_driver: Error #2: failure allocating kpath. stat = "//trim(adjustl(errormsg))//"."
      error stop trim(errormsg)
    endif

    ivec = 1
    isampling = 1

    do i = 1, nk - 1

      count = 0
      do icount = 1, nvec - 1
        count = count + nkpts(icount) - 1
        if (count >= i) then
          ivec = icount
          isampling = i - (count - (nkpts(icount) - 1))
          exit
        endif
      enddo

      kpath(:, i) = vecs(ivec, :) + &
                    (vecs(ivec + 1, :) - vecs(ivec, :))*real(isampling - 1, wp) &
                    /real(nkpts(ivec) - 1, wp)

    enddo
    ivec = nvec - 1
    isampling = nkpts(ivec)

    kpath(:, nk) = vecs(ivec, :) + &
                   (vecs(ivec + 1, :) - vecs(ivec, :))*real(isampling - 1, wp) &
                   /real(nkpts(ivec) - 1, wp)

  end function kpath

  function kslice(corner, vec_a, vec_b, part)
    real(wp), intent(in) :: corner(3), vec_a(3), vec_b(3)
    integer, intent(in) :: part(2)

    real(wp), allocatable :: kslice(:, :)

    integer :: i, ik(2)
    type(container_specifier) :: handle

    character(len=1024) :: errormsg
    integer :: istat

    do i = 1, 2
      if (part(i) < 2) &
        error stop "SsTC_driver: Error #1: components of part must be positive integers greater than 1."
    enddo
    call handle%specify(part)

    allocate (kslice(3, product(part)), stat=istat)
    if (istat /= 0) then
      write (errormsg, "(i20)") istat
      errormsg = "SsTC_driver: Error #2: failure allocating kslice. stat = "//trim(adjustl(errormsg))//"."
      error stop trim(errormsg)
    endif

    do i = 1, handle%size()
      ik = handle%ind(i)

      kslice(:, i) = corner + &
                     vec_a*real(ik(1) - 1, wp)/real(part(1) - 1, wp) + &
                     vec_b*real(ik(2) - 1, wp)/real(part(2) - 1, wp)

    enddo

  end function kslice

end module SsTC_driver_utilities
