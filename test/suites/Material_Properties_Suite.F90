module Material_Properties_Suite
  use, intrinsic :: iso_fortran_env, only: error_unit
  use MPI
  use OMP_LIB

  use SsTC_driver_kinds, only: wp => dp
  use SsTC_driver_definitions, only: cmplx_0, cmplx_i
  use SsTC_driver, only: task_specifier
  use MAC, only: container
  use WannInt, only: crystal, diagonalize, dirac_delta

  use testdrive, only: error_type

  implicit none
  private

  real(wp) :: tol = 1.0E10_wp

  public :: sample_and_sum_GaAs_optcond_in_peak

contains

  subroutine sample_and_sum_GaAs_optcond_in_peak(error)
    type(error_type), allocatable, intent(out) :: error

    type(crystal) :: GaAs
    type(task_specifier) :: optcond
    complex(wp), allocatable :: result(:, :)

    real(wp), parameter :: omega_start = 420.0_wp/99, &
                           omega_end = 420.0_wp/99
    integer, parameter  :: omega_steps = 1

    integer, dimension(3), parameter :: partition = [25, 25, 25]

    call GaAs%construct(name="GaAs", &
                        from_file="./material_data/GaAs_tb.dat", &
                        fermi_energy=7.7414_wp)

    call optcond%construct(name="optcond", &
                           int_ind=[3, 3], &
                           cont_data_start=[omega_start], &
                           cont_data_end=[omega_end], &
                           cont_data_steps=[omega_steps], &
                           calculator=optcond_calculator)

    call optcond%sample(crys=GaAs, kpart=partition, store_at=result, &
                        parallelization="MPI+OMP")
    result = result/product(partition)

    if ((abs(real(result(1, 1), wp) + 6.073940790512_wp) > tol*epsilon(1.0_wp))) &
      allocate (error)

  end subroutine sample_and_sum_GaAs_optcond_in_peak

  function optcond_calculator(self, crys, k, other)
    class(task_specifier), intent(in) :: self
    class(crystal), intent(in) :: crys
    real(wp), intent(in) :: k(3)
    class(*), optional, intent(in) :: other

    complex(wp) :: optcond_calculator(self%idims%size(), self%cdims%size())

    complex(wp) :: HW(crys%num_bands(), crys%num_bands()), &
                   DHW(crys%num_bands(), crys%num_bands(), 3), &
                   HH(crys%num_bands(), crys%num_bands()), &
                   AW(crys%num_bands(), crys%num_bands(), 3), &
                   AH(crys%num_bands(), crys%num_bands(), 3), &
                   NAD(crys%num_bands(), crys%num_bands(), 3), &
                   rot(crys%num_bands(), crys%num_bands()), &
                   optcs(self%cdims%size(), 3, 3)

    type(container) :: interpolator

    real(wp) :: eig(crys%num_bands()), omega, delta, smr
    integer :: occ(crys%num_bands())

    integer :: i, j, r, nbnd, n, m, iarr(2)

    real(wp), parameter :: deg_thr = 0.01_wp, &
                           deg_offset = 0.04_wp

    nbnd = crys%num_bands()
    !Smearing for the dirac delta in eV.
    smr = 0.1_wp

    !Obtain H in Hamiltonian basis.
    HW = crys%hamiltonian(kpt=k)
    call diagonalize(matrix=HW, P=rot, D=HH, eig=eig)

    !Obtain dH/dk and A in Wannier basis.
    interpolator = crys%hamiltonian(kpt=k, derivative=1)
    DHW = reshape(interpolator%cdp_storage, [nbnd, nbnd, 3])
    AW = crys%berry_connection(kpt=k)

    !Set up occupation lis.
    occ = 0
    do n = 1, nbnd
      if (eig(n) <= crys%fermi_energy()) occ(n) = 1
    enddo

    !Obtain NAD to transform A from Wannier to Hamiltonian
    !basis.
    do i = 1, 3
      NAD(:, :, i) = matmul(matmul(transpose(conjg(rot)), DHW(:, :, i)), rot)

      do n = 1, nbnd
        do m = 1, nbnd
          if (abs(eig(n) - eig(m)) < deg_thr) then
            NAD(n, m, i) = cmplx_0
          else
            NAD(n, m, i) = NAD(n, m, i)*((eig(m) - eig(n))/((eig(m) - eig(n))**2 + (deg_offset)**2))
          endif
        enddo
      enddo

    enddo

    !Transform connection to Hamiltonian basis.
    do i = 1, 3
      AH(:, :, i) = matmul(matmul(transpose(conjg(rot)), AW(:, :, i)), rot) + &
                    cmplx_i*NAD(:, :, i)
    enddo

    !Compute complex optical conductivity tensor.
    optcs = cmplx_0
    do j = 1, 3
      do i = 1, 3
        do r = 1, self%cdims%size()

          do n = 1, nbnd
            do m = 1, nbnd
              if (n == m) cycle
              omega = self%cdt(var=1, step=r)
              delta = dirac_delta(x=eig(m) - eig(n) - omega, smr=smr)
              optcs(r, i, j) = optcs(r, i, j) + &
                               (occ(n) - occ(m))*(eig(n) - eig(m))*AH(m, n, i)*AH(n, m, j)*delta
            enddo
          enddo

          !Parse to result.
          iarr = [i, j]
          optcond_calculator(self%idims%ind(iarr), self%cdims%ind(r)) = optcs(r, i, j)
        enddo
      enddo
    enddo

  end function optcond_calculator

end module Material_Properties_Suite
