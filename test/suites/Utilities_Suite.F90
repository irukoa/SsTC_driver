module Utilities_Suite
  use, intrinsic :: iso_fortran_env, only: error_unit

  use SsTC_driver_kinds, only: wp => dp
  use SsTC_driver_utils, only: kpath, kpath_length, kslice, &
    crys_to_cart, cart_to_crys
  use WannInt, only: crystal

  use testdrive, only: error_type

  implicit none
  private

  public :: test_kpath, test_kslice, test_kpath_length, test_crys_to_cart_cart_to_crys

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

  subroutine test_kpath_length(error)
    type(error_type), allocatable, intent(out) :: error

    type(crystal) :: GaAs
    integer, parameter :: nvec = 5
    real(wp) :: vecs(nvec, 3), total_length
    real(wp), allocatable :: path_length(:)

    call GaAs%construct(name="GaAs", &
                        from_file="./material_data/GaAs_tb.dat", &
                        fermi_energy=7.7414_wp)

    vecs(1, :) = [0.0_wp, 0.0_wp, 0.0_wp]
    vecs(2, :) = [0.5_wp, 0.0_wp, 0.0_wp]
    vecs(3, :) = [0.5_wp, 0.5_wp, 0.0_wp]
    vecs(4, :) = [0.0_wp, 0.5_wp, 0.0_wp]
    vecs(5, :) = [0.0_wp, 0.0_wp, 0.0_wp]

    path_length = kpath_length(vecs, [10, 10, 10, 10], GaAs)

    total_length = norm2(crys_to_cart(vecs(2, :) - vecs(1, :), GaAs)) + &
                   norm2(crys_to_cart(vecs(3, :) - vecs(2, :), GaAs)) + &
                   norm2(crys_to_cart(vecs(4, :) - vecs(3, :), GaAs)) + &
                   norm2(crys_to_cart(vecs(5, :) - vecs(4, :), GaAs))
    if (abs(total_length - path_length(size(path_length))) > tol*epsilon(1.0_wp)) allocate (error)

  end subroutine test_kpath_length

  subroutine test_kslice(error)
    type(error_type), allocatable, intent(out) :: error

    real(wp), allocatable :: slice(:, :)

    integer :: i

    slice = kslice([-0.5_wp, -0.5_wp, 0.0_wp], [1.0_wp, 0.0_wp, 0.0_wp], &
                   [0.0_wp, 1.0_wp, 0.0_wp], [10, 10])

    if (size(slice(1, :)) /= 100) allocate (error)

    if (abs(slice(1, 2) - (-0.5_wp + 1.0_wp/9.0_wp)) > tol*epsilon(1.0_wp)) allocate (error)

  end subroutine test_kslice

  subroutine test_crys_to_cart_cart_to_crys(error)
    type(error_type), allocatable, intent(out) :: error

    type(crystal) :: GaAs
    real(wp) :: point_cart(3), point_crys(3), diff(3)

    call random_seed()

    call GaAs%construct(name="GaAs", &
                        from_file="./material_data/GaAs_tb.dat", &
                        fermi_energy=7.7414_wp)

    call random_number(point_cart)
    point_cart = 100.0_wp*(point_cart - 0.5_wp)
    point_crys = cart_to_crys(point_cart, GaAs)

    diff = point_cart - crys_to_cart(point_crys, GaAs)
    if (norm2(diff) > 100*tol*epsilon(1.0_wp)) allocate (error)

  end subroutine test_crys_to_cart_cart_to_crys

end module Utilities_Suite
