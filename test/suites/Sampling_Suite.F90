module Sampling_Suite

  use MPI
  use OMP_LIB

  use SsTC_driver_kinds, only: wp => dp
  use SsTC_driver_definitions, only: cmplx_i
  use SsTC_driver, only: task_specifier
  use WannInt, only: crystal

  use testdrive, only: error_type

  implicit none
  private

  real(wp) :: tol = 1.0E4_wp

  public :: sample_no_sum_par_1, sample_no_sum_par_2, &
            sample_no_sum_par_3, sample_no_sum_par_4, &
            sample_sum_par_1, sample_sum_par_2, &
            sample_sum_par_3, sample_sum_par_4, &
            sample_no_sum_s_par_1, sample_no_sum_s_par_2, &
            sample_no_sum_s_par_3, sample_no_sum_s_par_4, &
            sample_sum_s_par_1, sample_sum_s_par_2, &
            sample_sum_s_par_3, sample_sum_s_par_4

contains

  subroutine sample_no_sum_par_1(error)
    type(error_type), allocatable, intent(out) :: error

    type(crystal) :: dummy
    type(task_specifier) :: tsk
    complex(wp), allocatable :: result(:, :, :)

    integer, parameter :: nk = 500

    real(wp) :: kpath(3, nk)
    integer :: i, kpart(3)
    real(wp) :: re, im

    do i = 1, nk
      kpath(:, i) = -[0.5_wp, 0.5_wp, 0.5_wp] + [1.0_wp, 1.0_wp, 1.0_wp]*real(i - 1, wp)/real(nk - 1, wp)
    enddo
    kpart = nint(real(nk, wp)**(1.0_wp/3.0_wp))

    call tsk%construct(name="test", &
                       int_ind=[3], &
                       cont_data_start=[0.0_wp], cont_data_end=[1.0_wp], cont_data_steps=[10], exponent=10.0_wp, &
                       calculator=test_sampler)

    call tsk%sample(dummy, kpath, result, "MPI+OMP")

    !re + i*im = 1*e^(-i*(10^0)*sqrt(3)/2).
    re = real(result(1, 1, 1), wp)
    im = aimag(result(1, 1, 1))

    if ((abs(re - cos(sqrt(3.0_wp)/2.0_wp)) > tol*epsilon(1.0_wp)) .or. &
        (abs(im + sin(sqrt(3.0_wp)/2.0_wp)) > tol*epsilon(1.0_wp))) &
      allocate (error)

    call tsk%sample(dummy, kpart, result, "MPI+OMP")

    !re + i*im = 1*e^(-i*(10^0)*sqrt(3)/2).
    re = real(result(1, 1, 1), wp)
    im = aimag(result(1, 1, 1))

    if ((abs(re - cos(sqrt(3.0_wp)/2.0_wp)) > tol*epsilon(1.0_wp)) .or. &
        (abs(im + sin(sqrt(3.0_wp)/2.0_wp)) > tol*epsilon(1.0_wp))) &
      allocate (error)

    call tsk%construct(name="test", &
                       int_ind=[3], &
                       cont_data_start=[0.0_wp], cont_data_end=[1.0_wp], cont_data_steps=[10], &
                       calculator=test_sampler)

    call tsk%sample(dummy, kpath, result, "MPI+OMP")

    !re + i*im = 1*e^(-i*0*sqrt(3)/2).
    re = real(result(1, 1, 1), wp)
    im = aimag(result(1, 1, 1))

    if ((abs(re - 1.0_wp) > tol*epsilon(1.0_wp)) .or. &
        (abs(im) > tol*epsilon(1.0_wp))) &
      allocate (error)

    call tsk%sample(dummy, kpart, result, "MPI+OMP")

    !re + i*im = 1*e^(-i*0*sqrt(3)/2).
    re = real(result(1, 1, 1), wp)
    im = aimag(result(1, 1, 1))

    if ((abs(re - 1.0_wp) > tol*epsilon(1.0_wp)) .or. &
        (abs(im) > tol*epsilon(1.0_wp))) &
      allocate (error)

  end subroutine sample_no_sum_par_1

  subroutine sample_no_sum_par_2(error)
    type(error_type), allocatable, intent(out) :: error

    type(crystal) :: dummy
    type(task_specifier) :: tsk
    complex(wp), allocatable :: result(:, :, :)

    integer, parameter :: nk = 500

    real(wp) :: kpath(3, nk)
    integer :: i, kpart(3)
    real(wp) :: re, im

    do i = 1, nk
      kpath(:, i) = -[0.5_wp, 0.5_wp, 0.5_wp] + [1.0_wp, 1.0_wp, 1.0_wp]*real(i - 1, wp)/real(nk - 1, wp)
    enddo
    kpart = nint(real(nk, wp)**(1.0_wp/3.0_wp))

    call tsk%construct(name="test", &
                       int_ind=[3], &
                       cont_data_start=[0.0_wp], cont_data_end=[1.0_wp], cont_data_steps=[10], exponent=10.0_wp, &
                       calculator=test_sampler)

    call tsk%sample(dummy, kpath, result, "MPI")

    !re + i*im = 1*e^(-i*(10^0)*sqrt(3)/2).
    re = real(result(1, 1, 1), wp)
    im = aimag(result(1, 1, 1))

    if ((abs(re - cos(sqrt(3.0_wp)/2.0_wp)) > tol*epsilon(1.0_wp)) .or. &
        (abs(im + sin(sqrt(3.0_wp)/2.0_wp)) > tol*epsilon(1.0_wp))) &
      allocate (error)

    call tsk%sample(dummy, kpart, result, "MPI")

    !re + i*im = 1*e^(-i*(10^0)*sqrt(3)/2).
    re = real(result(1, 1, 1), wp)
    im = aimag(result(1, 1, 1))

    if ((abs(re - cos(sqrt(3.0_wp)/2.0_wp)) > tol*epsilon(1.0_wp)) .or. &
        (abs(im + sin(sqrt(3.0_wp)/2.0_wp)) > tol*epsilon(1.0_wp))) &
      allocate (error)

    call tsk%construct(name="test", &
                       int_ind=[3], &
                       cont_data_start=[0.0_wp], cont_data_end=[1.0_wp], cont_data_steps=[10], &
                       calculator=test_sampler)

    call tsk%sample(dummy, kpath, result, "MPI")

    !re + i*im = 1*e^(-i*0*sqrt(3)/2).
    re = real(result(1, 1, 1), wp)
    im = aimag(result(1, 1, 1))

    if ((abs(re - 1.0_wp) > tol*epsilon(1.0_wp)) .or. &
        (abs(im) > tol*epsilon(1.0_wp))) &
      allocate (error)

    call tsk%sample(dummy, kpart, result, "MPI")

    !re + i*im = 1*e^(-i*0*sqrt(3)/2).
    re = real(result(1, 1, 1), wp)
    im = aimag(result(1, 1, 1))

    if ((abs(re - 1.0_wp) > tol*epsilon(1.0_wp)) .or. &
        (abs(im) > tol*epsilon(1.0_wp))) &
      allocate (error)

  end subroutine sample_no_sum_par_2

  subroutine sample_no_sum_par_3(error)
    type(error_type), allocatable, intent(out) :: error

    type(crystal) :: dummy
    type(task_specifier) :: tsk
    complex(wp), allocatable :: result(:, :, :)

    integer, parameter :: nk = 500

    real(wp) :: kpath(3, nk)
    integer :: i, kpart(3)
    real(wp) :: re, im

    do i = 1, nk
      kpath(:, i) = -[0.5_wp, 0.5_wp, 0.5_wp] + [1.0_wp, 1.0_wp, 1.0_wp]*real(i - 1, wp)/real(nk - 1, wp)
    enddo
    kpart = nint(real(nk, wp)**(1.0_wp/3.0_wp))

    call tsk%construct(name="test", &
                       int_ind=[3], &
                       cont_data_start=[0.0_wp], cont_data_end=[1.0_wp], cont_data_steps=[10], exponent=10.0_wp, &
                       calculator=test_sampler)

    call tsk%sample(dummy, kpath, result, "OMP")

    !re + i*im = 1*e^(-i*(10^0)*sqrt(3)/2).
    re = real(result(1, 1, 1), wp)
    im = aimag(result(1, 1, 1))

    if ((abs(re - cos(sqrt(3.0_wp)/2.0_wp)) > tol*epsilon(1.0_wp)) .or. &
        (abs(im + sin(sqrt(3.0_wp)/2.0_wp)) > tol*epsilon(1.0_wp))) &
      allocate (error)

    call tsk%sample(dummy, kpart, result, "OMP")

    !re + i*im = 1*e^(-i*(10^0)*sqrt(3)/2).
    re = real(result(1, 1, 1), wp)
    im = aimag(result(1, 1, 1))

    if ((abs(re - cos(sqrt(3.0_wp)/2.0_wp)) > tol*epsilon(1.0_wp)) .or. &
        (abs(im + sin(sqrt(3.0_wp)/2.0_wp)) > tol*epsilon(1.0_wp))) &
      allocate (error)

    call tsk%construct(name="test", &
                       int_ind=[3], &
                       cont_data_start=[0.0_wp], cont_data_end=[1.0_wp], cont_data_steps=[10], &
                       calculator=test_sampler)

    call tsk%sample(dummy, kpath, result, "OMP")

    !re + i*im = 1*e^(-i*0*sqrt(3)/2).
    re = real(result(1, 1, 1), wp)
    im = aimag(result(1, 1, 1))

    if ((abs(re - 1.0_wp) > tol*epsilon(1.0_wp)) .or. &
        (abs(im) > tol*epsilon(1.0_wp))) &
      allocate (error)

    call tsk%sample(dummy, kpart, result, "OMP")

    !re + i*im = 1*e^(-i*0*sqrt(3)/2).
    re = real(result(1, 1, 1), wp)
    im = aimag(result(1, 1, 1))

    if ((abs(re - 1.0_wp) > tol*epsilon(1.0_wp)) .or. &
        (abs(im) > tol*epsilon(1.0_wp))) &
      allocate (error)

  end subroutine sample_no_sum_par_3

  subroutine sample_no_sum_par_4(error)
    type(error_type), allocatable, intent(out) :: error

    type(crystal) :: dummy
    type(task_specifier) :: tsk
    complex(wp), allocatable :: result(:, :, :)

    integer, parameter :: nk = 500

    real(wp) :: kpath(3, nk)
    integer :: i, kpart(3)
    real(wp) :: re, im

    do i = 1, nk
      kpath(:, i) = -[0.5_wp, 0.5_wp, 0.5_wp] + [1.0_wp, 1.0_wp, 1.0_wp]*real(i - 1, wp)/real(nk - 1, wp)
    enddo
    kpart = nint(real(nk, wp)**(1.0_wp/3.0_wp))

    call tsk%construct(name="test", &
                       int_ind=[3], &
                       cont_data_start=[0.0_wp], cont_data_end=[1.0_wp], cont_data_steps=[10], exponent=10.0_wp, &
                       calculator=test_sampler)

    call tsk%sample(dummy, kpath, result, "none")

    !re + i*im = 1*e^(-i*(10^0)*sqrt(3)/2).
    re = real(result(1, 1, 1), wp)
    im = aimag(result(1, 1, 1))

    if ((abs(re - cos(sqrt(3.0_wp)/2.0_wp)) > tol*epsilon(1.0_wp)) .or. &
        (abs(im + sin(sqrt(3.0_wp)/2.0_wp)) > tol*epsilon(1.0_wp))) &
      allocate (error)

    call tsk%sample(dummy, kpart, result, "none")

    !re + i*im = 1*e^(-i*(10^0)*sqrt(3)/2).
    re = real(result(1, 1, 1), wp)
    im = aimag(result(1, 1, 1))

    if ((abs(re - cos(sqrt(3.0_wp)/2.0_wp)) > tol*epsilon(1.0_wp)) .or. &
        (abs(im + sin(sqrt(3.0_wp)/2.0_wp)) > tol*epsilon(1.0_wp))) &
      allocate (error)

    call tsk%construct(name="test", &
                       int_ind=[3], &
                       cont_data_start=[0.0_wp], cont_data_end=[1.0_wp], cont_data_steps=[10], &
                       calculator=test_sampler)

    call tsk%sample(dummy, kpath, result, "none")

    !re + i*im = 1*e^(-i*0*sqrt(3)/2).
    re = real(result(1, 1, 1), wp)
    im = aimag(result(1, 1, 1))

    if ((abs(re - 1.0_wp) > tol*epsilon(1.0_wp)) .or. &
        (abs(im) > tol*epsilon(1.0_wp))) &
      allocate (error)

    call tsk%sample(dummy, kpart, result, "none")

    !re + i*im = 1*e^(-i*0*sqrt(3)/2).
    re = real(result(1, 1, 1), wp)
    im = aimag(result(1, 1, 1))

    if ((abs(re - 1.0_wp) > tol*epsilon(1.0_wp)) .or. &
        (abs(im) > tol*epsilon(1.0_wp))) &
      allocate (error)

  end subroutine sample_no_sum_par_4

  subroutine sample_sum_par_1(error)
    type(error_type), allocatable, intent(out) :: error

    type(crystal) :: dummy
    type(task_specifier) :: tsk
    complex(wp), allocatable :: result(:, :)

    integer, parameter :: nk = 500

    real(wp) :: kpath(3, nk)
    integer :: i, kpart(3)
    real(wp) :: re, im

    do i = 1, nk
      kpath(:, i) = -[0.5_wp, 0.5_wp, 0.5_wp] + [1.0_wp, 1.0_wp, 1.0_wp]*real(i - 1, wp)/real(nk - 1, wp)
    enddo
    kpart = nint(real(nk, wp)**(1.0_wp/3.0_wp))

    call tsk%construct(name="test", &
                       int_ind=[3], &
                       cont_data_start=[0.0_wp], cont_data_end=[1.0_wp], cont_data_steps=[10], exponent=10.0_wp, &
                       calculator=test_sampler)

    call tsk%sample(dummy, kpath, result, "MPI+OMP")

    re = real(result(1, 1), wp)
    im = aimag(result(1, 1))

    if ((abs(re - 439.57014432898791_wp) > tol*epsilon(1.0_wp)) .or. &
        (abs(im + 203.66430906454599_wp) > tol*epsilon(1.0_wp))) &
      allocate (error)

    call tsk%sample(dummy, kpart, result, "MPI+OMP")

    re = real(result(1, 1), wp)
    im = aimag(result(1, 1))

    if ((abs(re - 432.43603843687964_wp) > tol*epsilon(1.0_wp)) .or. &
        (abs(im + 262.53892640642482_wp) > tol*epsilon(1.0_wp))) &
      allocate (error)

    call tsk%construct(name="test", &
                       int_ind=[3], &
                       cont_data_start=[0.0_wp], cont_data_end=[1.0_wp], cont_data_steps=[10], &
                       calculator=test_sampler)

    call tsk%sample(dummy, kpath, result, "MPI+OMP")

    re = real(result(1, 1), wp)
    im = aimag(result(1, 1))

    if ((abs(re - 500.0_wp) > tol*epsilon(1.0_wp)) .or. &
        (abs(im) > tol*epsilon(1.0_wp))) &
      allocate (error)

    call tsk%sample(dummy, kpart, result, "MPI+OMP")

    re = real(result(1, 1), wp)
    im = aimag(result(1, 1))

    if ((abs(re - real(nint(real(nk, wp)**(1.0_wp/3.0_wp)), wp)**3) > tol*epsilon(1.0_wp)) .or. &
        (abs(im) > tol*epsilon(1.0_wp))) &
      allocate (error)

  end subroutine sample_sum_par_1

  subroutine sample_sum_par_2(error)
    type(error_type), allocatable, intent(out) :: error

    type(crystal) :: dummy
    type(task_specifier) :: tsk
    complex(wp), allocatable :: result(:, :)

    integer, parameter :: nk = 500

    real(wp) :: kpath(3, nk)
    integer :: i, kpart(3)
    real(wp) :: re, im

    do i = 1, nk
      kpath(:, i) = -[0.5_wp, 0.5_wp, 0.5_wp] + [1.0_wp, 1.0_wp, 1.0_wp]*real(i - 1, wp)/real(nk - 1, wp)
    enddo
    kpart = nint(real(nk, wp)**(1.0_wp/3.0_wp))

    call tsk%construct(name="test", &
                       int_ind=[3], &
                       cont_data_start=[0.0_wp], cont_data_end=[1.0_wp], cont_data_steps=[10], exponent=10.0_wp, &
                       calculator=test_sampler)

    call tsk%sample(dummy, kpath, result, "MPI")

    re = real(result(1, 1), wp)
    im = aimag(result(1, 1))

    if ((abs(re - 439.57014432898791_wp) > tol*epsilon(1.0_wp)) .or. &
        (abs(im + 203.66430906454599_wp) > tol*epsilon(1.0_wp))) &
      allocate (error)

    call tsk%sample(dummy, kpart, result, "MPI")

    re = real(result(1, 1), wp)
    im = aimag(result(1, 1))

    if ((abs(re - 432.43603843687964_wp) > tol*epsilon(1.0_wp)) .or. &
        (abs(im + 262.53892640642482_wp) > tol*epsilon(1.0_wp))) &
      allocate (error)

    call tsk%construct(name="test", &
                       int_ind=[3], &
                       cont_data_start=[0.0_wp], cont_data_end=[1.0_wp], cont_data_steps=[10], &
                       calculator=test_sampler)

    call tsk%sample(dummy, kpath, result, "MPI")

    re = real(result(1, 1), wp)
    im = aimag(result(1, 1))

    if ((abs(re - 500.0_wp) > tol*epsilon(1.0_wp)) .or. &
        (abs(im) > tol*epsilon(1.0_wp))) &
      allocate (error)

    call tsk%sample(dummy, kpart, result, "MPI")

    re = real(result(1, 1), wp)
    im = aimag(result(1, 1))

    if ((abs(re - real(nint(real(nk, wp)**(1.0_wp/3.0_wp)), wp)**3) > tol*epsilon(1.0_wp)) .or. &
        (abs(im) > tol*epsilon(1.0_wp))) &
      allocate (error)

  end subroutine sample_sum_par_2

  subroutine sample_sum_par_3(error)
    type(error_type), allocatable, intent(out) :: error

    type(crystal) :: dummy
    type(task_specifier) :: tsk
    complex(wp), allocatable :: result(:, :)

    integer, parameter :: nk = 500

    real(wp) :: kpath(3, nk)
    integer :: i, kpart(3)
    real(wp) :: re, im

    do i = 1, nk
      kpath(:, i) = -[0.5_wp, 0.5_wp, 0.5_wp] + [1.0_wp, 1.0_wp, 1.0_wp]*real(i - 1, wp)/real(nk - 1, wp)
    enddo
    kpart = nint(real(nk, wp)**(1.0_wp/3.0_wp))

    call tsk%construct(name="test", &
                       int_ind=[3], &
                       cont_data_start=[0.0_wp], cont_data_end=[1.0_wp], cont_data_steps=[10], exponent=10.0_wp, &
                       calculator=test_sampler)

    call tsk%sample(dummy, kpath, result, "OMP")

    re = real(result(1, 1), wp)
    im = aimag(result(1, 1))

    if ((abs(re - 439.57014432898791_wp) > tol*epsilon(1.0_wp)) .or. &
        (abs(im + 203.66430906454599_wp) > tol*epsilon(1.0_wp))) &
      allocate (error)

    call tsk%sample(dummy, kpart, result, "OMP")

    re = real(result(1, 1), wp)
    im = aimag(result(1, 1))

    if ((abs(re - 432.43603843687964_wp) > tol*epsilon(1.0_wp)) .or. &
        (abs(im + 262.53892640642482_wp) > tol*epsilon(1.0_wp))) &
      allocate (error)

    call tsk%construct(name="test", &
                       int_ind=[3], &
                       cont_data_start=[0.0_wp], cont_data_end=[1.0_wp], cont_data_steps=[10], &
                       calculator=test_sampler)

    call tsk%sample(dummy, kpath, result, "OMP")

    re = real(result(1, 1), wp)
    im = aimag(result(1, 1))

    if ((abs(re - 500.0_wp) > tol*epsilon(1.0_wp)) .or. &
        (abs(im) > tol*epsilon(1.0_wp))) &
      allocate (error)

    call tsk%sample(dummy, kpart, result, "OMP")

    re = real(result(1, 1), wp)
    im = aimag(result(1, 1))

    if ((abs(re - real(nint(real(nk, wp)**(1.0_wp/3.0_wp)), wp)**3) > tol*epsilon(1.0_wp)) .or. &
        (abs(im) > tol*epsilon(1.0_wp))) &
      allocate (error)

  end subroutine sample_sum_par_3

  subroutine sample_sum_par_4(error)
    type(error_type), allocatable, intent(out) :: error

    type(crystal) :: dummy
    type(task_specifier) :: tsk
    complex(wp), allocatable :: result(:, :)

    integer, parameter :: nk = 500

    real(wp) :: kpath(3, nk)
    integer :: i, kpart(3)
    real(wp) :: re, im

    do i = 1, nk
      kpath(:, i) = -[0.5_wp, 0.5_wp, 0.5_wp] + [1.0_wp, 1.0_wp, 1.0_wp]*real(i - 1, wp)/real(nk - 1, wp)
    enddo
    kpart = nint(real(nk, wp)**(1.0_wp/3.0_wp))

    call tsk%construct(name="test", &
                       int_ind=[3], &
                       cont_data_start=[0.0_wp], cont_data_end=[1.0_wp], cont_data_steps=[10], exponent=10.0_wp, &
                       calculator=test_sampler)

    call tsk%sample(dummy, kpath, result, "none")

    re = real(result(1, 1), wp)
    im = aimag(result(1, 1))

    if ((abs(re - 439.57014432898791_wp) > tol*epsilon(1.0_wp)) .or. &
        (abs(im + 203.66430906454599_wp) > tol*epsilon(1.0_wp))) &
      allocate (error)
    return

    call tsk%sample(dummy, kpart, result, "none")

    re = real(result(1, 1), wp)
    im = aimag(result(1, 1))

    if ((abs(re - 432.43603843687964_wp) > tol*epsilon(1.0_wp)) .or. &
        (abs(im + 262.53892640642482_wp) > tol*epsilon(1.0_wp))) &
      allocate (error)
    return

    call tsk%construct(name="test", &
                       int_ind=[3], &
                       cont_data_start=[0.0_wp], cont_data_end=[1.0_wp], cont_data_steps=[10], &
                       calculator=test_sampler)

    call tsk%sample(dummy, kpath, result, "none")

    re = real(result(1, 1), wp)
    im = aimag(result(1, 1))

    if ((abs(re - 500.0_wp) > tol*epsilon(1.0_wp)) .or. &
        (abs(im) > tol*epsilon(1.0_wp))) &
      allocate (error)
    return

    call tsk%sample(dummy, kpart, result, "none")

    re = real(result(1, 1), wp)
    im = aimag(result(1, 1))

    if ((abs(re - real(nint(real(nk, wp)**(1.0_wp/3.0_wp)), wp)**3) > tol*epsilon(1.0_wp)) .or. &
        (abs(im) > tol*epsilon(1.0_wp))) &
      allocate (error)
    return

  end subroutine sample_sum_par_4

  subroutine sample_no_sum_s_par_1(error)
    type(error_type), allocatable, intent(out) :: error

    type(crystal) :: dummy
    type(task_specifier) :: tsk
    complex(wp), allocatable :: result(:, :, :)

    integer, parameter :: nk = 500

    real(wp) :: kpath(3, nk)
    integer :: i, kpart(3)
    real(wp) :: re, im

    do i = 1, nk
      kpath(:, i) = -[0.5_wp, 0.5_wp, 0.5_wp] + [1.0_wp, 1.0_wp, 1.0_wp]*real(i - 1, wp)/real(nk - 1, wp)
    enddo
    kpart = nint(real(nk, wp)**(1.0_wp/3.0_wp))

    call tsk%construct(name="test", &
                       calculator=test_sampler_simple)

    call tsk%sample(dummy, kpath, result, "MPI+OMP")

    !re + i*im = 1*e^(-i*(10^0)*sqrt(3)/2).
    re = real(result(1, 1, 1), wp)
    im = aimag(result(1, 1, 1))

    if ((abs(re - cos(sqrt(3.0_wp)/2.0_wp)) > tol*epsilon(1.0_wp)) .or. &
        (abs(im + sin(sqrt(3.0_wp)/2.0_wp)) > tol*epsilon(1.0_wp))) &
      allocate (error)

    call tsk%sample(dummy, kpart, result, "MPI+OMP")

    !re + i*im = 1*e^(-i*(10^0)*sqrt(3)/2).
    re = real(result(1, 1, 1), wp)
    im = aimag(result(1, 1, 1))

    if ((abs(re - cos(sqrt(3.0_wp)/2.0_wp)) > tol*epsilon(1.0_wp)) .or. &
        (abs(im + sin(sqrt(3.0_wp)/2.0_wp)) > tol*epsilon(1.0_wp))) &
      allocate (error)

  end subroutine sample_no_sum_s_par_1

  subroutine sample_no_sum_s_par_2(error)
    type(error_type), allocatable, intent(out) :: error

    type(crystal) :: dummy
    type(task_specifier) :: tsk
    complex(wp), allocatable :: result(:, :, :)

    integer, parameter :: nk = 500

    real(wp) :: kpath(3, nk)
    integer :: i, kpart(3)
    real(wp) :: re, im

    do i = 1, nk
      kpath(:, i) = -[0.5_wp, 0.5_wp, 0.5_wp] + [1.0_wp, 1.0_wp, 1.0_wp]*real(i - 1, wp)/real(nk - 1, wp)
    enddo
    kpart = nint(real(nk, wp)**(1.0_wp/3.0_wp))

    call tsk%construct(name="test", &
                       calculator=test_sampler_simple)

    call tsk%sample(dummy, kpath, result, "MPI")

    !re + i*im = 1*e^(-i*(10^0)*sqrt(3)/2).
    re = real(result(1, 1, 1), wp)
    im = aimag(result(1, 1, 1))

    if ((abs(re - cos(sqrt(3.0_wp)/2.0_wp)) > tol*epsilon(1.0_wp)) .or. &
        (abs(im + sin(sqrt(3.0_wp)/2.0_wp)) > tol*epsilon(1.0_wp))) &
      allocate (error)

    call tsk%sample(dummy, kpart, result, "MPI")

    !re + i*im = 1*e^(-i*(10^0)*sqrt(3)/2).
    re = real(result(1, 1, 1), wp)
    im = aimag(result(1, 1, 1))

    if ((abs(re - cos(sqrt(3.0_wp)/2.0_wp)) > tol*epsilon(1.0_wp)) .or. &
        (abs(im + sin(sqrt(3.0_wp)/2.0_wp)) > tol*epsilon(1.0_wp))) &
      allocate (error)

  end subroutine sample_no_sum_s_par_2

  subroutine sample_no_sum_s_par_3(error)
    type(error_type), allocatable, intent(out) :: error

    type(crystal) :: dummy
    type(task_specifier) :: tsk
    complex(wp), allocatable :: result(:, :, :)

    integer, parameter :: nk = 500

    real(wp) :: kpath(3, nk)
    integer :: i, kpart(3)
    real(wp) :: re, im

    do i = 1, nk
      kpath(:, i) = -[0.5_wp, 0.5_wp, 0.5_wp] + [1.0_wp, 1.0_wp, 1.0_wp]*real(i - 1, wp)/real(nk - 1, wp)
    enddo
    kpart = nint(real(nk, wp)**(1.0_wp/3.0_wp))

    call tsk%construct(name="test", &
                       calculator=test_sampler_simple)

    call tsk%sample(dummy, kpath, result, "OMP")

    !re + i*im = 1*e^(-i*(10^0)*sqrt(3)/2).
    re = real(result(1, 1, 1), wp)
    im = aimag(result(1, 1, 1))

    if ((abs(re - cos(sqrt(3.0_wp)/2.0_wp)) > tol*epsilon(1.0_wp)) .or. &
        (abs(im + sin(sqrt(3.0_wp)/2.0_wp)) > tol*epsilon(1.0_wp))) &
      allocate (error)

    call tsk%sample(dummy, kpart, result, "OMP")

    !re + i*im = 1*e^(-i*(10^0)*sqrt(3)/2).
    re = real(result(1, 1, 1), wp)
    im = aimag(result(1, 1, 1))

    if ((abs(re - cos(sqrt(3.0_wp)/2.0_wp)) > tol*epsilon(1.0_wp)) .or. &
        (abs(im + sin(sqrt(3.0_wp)/2.0_wp)) > tol*epsilon(1.0_wp))) &
      allocate (error)

  end subroutine sample_no_sum_s_par_3

  subroutine sample_no_sum_s_par_4(error)
    type(error_type), allocatable, intent(out) :: error

    type(crystal) :: dummy
    type(task_specifier) :: tsk
    complex(wp), allocatable :: result(:, :, :)

    integer, parameter :: nk = 500

    real(wp) :: kpath(3, nk)
    integer :: i, kpart(3)
    real(wp) :: re, im

    do i = 1, nk
      kpath(:, i) = -[0.5_wp, 0.5_wp, 0.5_wp] + [1.0_wp, 1.0_wp, 1.0_wp]*real(i - 1, wp)/real(nk - 1, wp)
    enddo
    kpart = nint(real(nk, wp)**(1.0_wp/3.0_wp))

    call tsk%construct(name="test", &
                       calculator=test_sampler_simple)

    call tsk%sample(dummy, kpath, result, "none")

    !re + i*im = 1*e^(-i*(10^0)*sqrt(3)/2).
    re = real(result(1, 1, 1), wp)
    im = aimag(result(1, 1, 1))

    if ((abs(re - cos(sqrt(3.0_wp)/2.0_wp)) > tol*epsilon(1.0_wp)) .or. &
        (abs(im + sin(sqrt(3.0_wp)/2.0_wp)) > tol*epsilon(1.0_wp))) &
      allocate (error)

    call tsk%sample(dummy, kpart, result, "none")

    !re + i*im = 1*e^(-i*(10^0)*sqrt(3)/2).
    re = real(result(1, 1, 1), wp)
    im = aimag(result(1, 1, 1))

    if ((abs(re - cos(sqrt(3.0_wp)/2.0_wp)) > tol*epsilon(1.0_wp)) .or. &
        (abs(im + sin(sqrt(3.0_wp)/2.0_wp)) > tol*epsilon(1.0_wp))) &
      allocate (error)

  end subroutine sample_no_sum_s_par_4

  subroutine sample_sum_s_par_1(error)
    type(error_type), allocatable, intent(out) :: error

    type(crystal) :: dummy
    type(task_specifier) :: tsk
    complex(wp), allocatable :: result(:, :)

    integer, parameter :: nk = 500

    real(wp) :: kpath(3, nk)
    integer :: i, kpart(3)
    real(wp) :: re, im

    do i = 1, nk
      kpath(:, i) = -[0.5_wp, 0.5_wp, 0.5_wp] + [1.0_wp, 1.0_wp, 1.0_wp]*real(i - 1, wp)/real(nk - 1, wp)
    enddo
    kpart = nint(real(nk, wp)**(1.0_wp/3.0_wp))

    call tsk%construct(name="test", &
                       calculator=test_sampler_simple)

    call tsk%sample(dummy, kpath, result, "MPI+OMP")

    re = real(result(1, 1), wp)
    im = aimag(result(1, 1))

    if ((abs(re - 439.57014432898791_wp) > tol*epsilon(1.0_wp)) .or. &
        (abs(im + 203.66430906454599_wp) > tol*epsilon(1.0_wp))) &
      allocate (error)

    call tsk%sample(dummy, kpart, result, "MPI+OMP")

    re = real(result(1, 1), wp)
    im = aimag(result(1, 1))

    if ((abs(re - 432.43603843687964_wp) > tol*epsilon(1.0_wp)) .or. &
        (abs(im + 262.53892640642482_wp) > tol*epsilon(1.0_wp))) &
      allocate (error)

  end subroutine sample_sum_s_par_1

  subroutine sample_sum_s_par_2(error)
    type(error_type), allocatable, intent(out) :: error

    type(crystal) :: dummy
    type(task_specifier) :: tsk
    complex(wp), allocatable :: result(:, :)

    integer, parameter :: nk = 500

    real(wp) :: kpath(3, nk)
    integer :: i, kpart(3)
    real(wp) :: re, im

    do i = 1, nk
      kpath(:, i) = -[0.5_wp, 0.5_wp, 0.5_wp] + [1.0_wp, 1.0_wp, 1.0_wp]*real(i - 1, wp)/real(nk - 1, wp)
    enddo
    kpart = nint(real(nk, wp)**(1.0_wp/3.0_wp))

    call tsk%construct(name="test", &
                       calculator=test_sampler_simple)

    call tsk%sample(dummy, kpath, result, "MPI")

    re = real(result(1, 1), wp)
    im = aimag(result(1, 1))

    if ((abs(re - 439.57014432898791_wp) > tol*epsilon(1.0_wp)) .or. &
        (abs(im + 203.66430906454599_wp) > tol*epsilon(1.0_wp))) &
      allocate (error)

    call tsk%sample(dummy, kpart, result, "MPI")

    re = real(result(1, 1), wp)
    im = aimag(result(1, 1))

    if ((abs(re - 432.43603843687964_wp) > tol*epsilon(1.0_wp)) .or. &
        (abs(im + 262.53892640642482_wp) > tol*epsilon(1.0_wp))) &
      allocate (error)

  end subroutine sample_sum_s_par_2

  subroutine sample_sum_s_par_3(error)
    type(error_type), allocatable, intent(out) :: error

    type(crystal) :: dummy
    type(task_specifier) :: tsk
    complex(wp), allocatable :: result(:, :)

    integer, parameter :: nk = 500

    real(wp) :: kpath(3, nk)
    integer :: i, kpart(3)
    real(wp) :: re, im

    do i = 1, nk
      kpath(:, i) = -[0.5_wp, 0.5_wp, 0.5_wp] + [1.0_wp, 1.0_wp, 1.0_wp]*real(i - 1, wp)/real(nk - 1, wp)
    enddo
    kpart = nint(real(nk, wp)**(1.0_wp/3.0_wp))

    call tsk%construct(name="test", &
                       calculator=test_sampler_simple)

    call tsk%sample(dummy, kpath, result, "OMP")

    re = real(result(1, 1), wp)
    im = aimag(result(1, 1))

    if ((abs(re - 439.57014432898791_wp) > tol*epsilon(1.0_wp)) .or. &
        (abs(im + 203.66430906454599_wp) > tol*epsilon(1.0_wp))) &
      allocate (error)

    call tsk%sample(dummy, kpart, result, "OMP")

    re = real(result(1, 1), wp)
    im = aimag(result(1, 1))

    if ((abs(re - 432.43603843687964_wp) > tol*epsilon(1.0_wp)) .or. &
        (abs(im + 262.53892640642482_wp) > tol*epsilon(1.0_wp))) &
      allocate (error)

  end subroutine sample_sum_s_par_3

  subroutine sample_sum_s_par_4(error)
    type(error_type), allocatable, intent(out) :: error

    type(crystal) :: dummy
    type(task_specifier) :: tsk
    complex(wp), allocatable :: result(:, :)

    integer, parameter :: nk = 500

    real(wp) :: kpath(3, nk)
    integer :: i, kpart(3)
    real(wp) :: re, im

    do i = 1, nk
      kpath(:, i) = -[0.5_wp, 0.5_wp, 0.5_wp] + [1.0_wp, 1.0_wp, 1.0_wp]*real(i - 1, wp)/real(nk - 1, wp)
    enddo
    kpart = nint(real(nk, wp)**(1.0_wp/3.0_wp))

    call tsk%construct(name="test", &
                       calculator=test_sampler_simple)

    call tsk%sample(dummy, kpath, result, "none")

    re = real(result(1, 1), wp)
    im = aimag(result(1, 1))

    if ((abs(re - 439.57014432898791_wp) > tol*epsilon(1.0_wp)) .or. &
        (abs(im + 203.66430906454599_wp) > tol*epsilon(1.0_wp))) &
      allocate (error)
    return

    call tsk%sample(dummy, kpart, result, "none")

    re = real(result(1, 1), wp)
    im = aimag(result(1, 1))

    if ((abs(re - 432.43603843687964_wp) > tol*epsilon(1.0_wp)) .or. &
        (abs(im + 262.53892640642482_wp) > tol*epsilon(1.0_wp))) &
      allocate (error)
    return

  end subroutine sample_sum_s_par_4

  function test_sampler(self, crys, k, other)
    class(task_specifier), intent(in) :: self
    class(crystal), intent(in) :: crys
    real(wp), intent(in) :: k(3)
    class(*), optional, intent(in) :: other

    logical :: dbg1, dbg2

    complex(wp) :: test_sampler(self%idims%size(), self%cdims%size())

    integer :: i, r

    dbg1 = crys%initialized()
    dbg2 = present(other)

    do i = 1, self%idims%size()
      do r = 1, self%cdims%size()
        test_sampler(i, r) = &
          real(i, wp)*exp(-cmplx_i*self%cdt(var=1, step=r)*norm2(k))
      enddo
    enddo

  end function test_sampler

  function test_sampler_simple(self, crys, k, other)
    class(task_specifier), intent(in) :: self
    class(crystal), intent(in) :: crys
    real(wp), intent(in) :: k(3)
    class(*), optional, intent(in) :: other

    logical :: dbg1, dbg2

    complex(wp) :: test_sampler_simple(self%idims%size(), self%cdims%size())

    integer :: i, r

    dbg1 = crys%initialized()
    dbg2 = present(other)

    do i = 1, self%idims%size()
      do r = 1, self%cdims%size()
        test_sampler_simple(i, r) = &
          exp(-cmplx_i*norm2(k))
      enddo
    enddo

  end function test_sampler_simple

end module Sampling_Suite
