module Utilities_Suite
  use, intrinsic :: iso_fortran_env, only: error_unit

  use SsTC_driver_kinds, only: wp => dp
  use SsTC_driver_utils, only: kpath, kslice

  use testdrive, only: error_type

  public :: test_kpath, test_kslice

  real(wp) :: tol = 1.0E1_wp

contains

  subroutine test_kpath(error)
    type(error_type), allocatable, intent(out) :: error

    integer, parameter :: nvec = 5
    real(wp) :: vecs(nvec, 3)
    real(wp), allocatable :: path(:, :)

    vecs(1, :) = [0.0_wp, 0.0_wp, 0.0_wp]
    vecs(2, :) = [0.5_wp, 0.0_wp, 0.0_wp]
    vecs(3, :) = [0.5_wp, 0.5_wp, 0.0_wp]
    vecs(4, :) = [0.0_wp, 0.5_wp, 0.0_wp]
    vecs(5, :) = [0.0_wp, 0.0_wp, 0.0_wp]

    path = kpath(vecs, [10, 10, 10, 10])

    if (size(path(1, :)) /= 40 - 5 + 2) allocate (error)

    if (abs(path(1, 3) - 1.0_wp/9.0_wp) > tol*epsilon(1.0_wp)) allocate (error)

  end subroutine test_kpath

  subroutine test_kslice(error)
    type(error_type), allocatable, intent(out) :: error

    real(wp), allocatable :: slice(:, :)

    integer :: i

    slice = kslice([-0.5_wp, -0.5_wp, 0.0_wp], [1.0_wp, 0.0_wp, 0.0_wp], &
                   [0.0_wp, 1.0_wp, 0.0_wp], [10, 10])

    if (size(slice(1, :)) /= 100) allocate (error)

    if (abs(slice(1, 2) - (-0.5_wp + 1.0_wp/9.0_wp)) > tol*epsilon(1.0_wp)) allocate (error)

  end subroutine test_kslice

end module Utilities_Suite
