program GaAs_Complex_Optical_Conductivity
  use, intrinsic :: iso_fortran_env, only: output_unit
  use MPI
  use OMP_LIB

  use SsTC_driver_kinds, only: wp => dp
  use SsTC_driver_definitions, only: cmplx_0, cmplx_i
  use SsTC_driver, only: task_specifier
  use MAC, only: container
  use WannInt, only: crystal, diagonalize, dirac_delta

  !In this example we calculate the complex optical conductivity
  !tensor \eps^{ij}(omega) of GaAs.
  !Using symmetry arguments, it can be shown that
  !for the space group of GaAs, only the \eps^{11}(omega) =
  !\eps^{22}(omega) = \eps^{33}(omega) are nonzero.

  implicit none

  type(crystal) :: GaAs
  type(task_specifier) :: optcond
  complex(wp), allocatable :: result(:, :)

  integer :: i, j, iarr(2), r, out
  character(len=120) :: filename, ci, cj

  integer :: ierror, rank

  !Tell the program which frequency
  !range we will be sampling (in eV).
  real(wp), parameter :: omega_start = 0.0_wp, &
                         omega_end = 5.0_wp
  integer, parameter  :: omega_steps = 100

  integer, dimension(3), parameter :: partition = [100, 100, 100]

  !MPI things.
  call MPI_INIT(ierror)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)

  !Load data for GaAs from file using WannInt.
  call GaAs%construct(name="GaAs", &
                      from_file="./material_data/GaAs_tb.dat", &
                      fermi_energy=7.7414_wp)

  !Here we describe the functional dependence of the optical conductivity:
  !2 integer indices (i, j) ranging from each from 1 to 3,
  !and a continuous index (omega) starting at omega_start, ending at
  !omega_end and discretized in omega_steps steps.
  !We also point the task where the calculator is to be found.
  call optcond%construct(name="Optical_Conductivity", &
                         int_ind=[3, 3], &
                         cont_data_start=[omega_start], &
                         cont_data_end=[omega_end], &
                         cont_data_steps=[omega_steps], &
                         calculator=optcond_calculator)

  !Here we do the actual sampling. We tell the program to partition the BZ
  !and sample. Notice that since result is a rank-2 array the contributions
  !from all k will be implicitly summed over. Thats why we then divide.
  call optcond%sample(crys=GaAs, kpart=partition, store_at=result, &
                      parallelization="MPI+OMP")
  result = result/product(partition)
  !Now result holds the optical conductivity tensor for GaAs in A^{-1}.
  !We still need to multiply by constants, but for the present example
  !it suffices to print results already.

  !Lastly, we print to files.
  !You can visualize the spectra using gnuplot: gnuplot> p 'GaAs-optcond_11.dat' w l
  do j = 1, 3
    write (cj, "(i0)") j
    do i = 1, 3
      write (ci, "(i0)") i
      filename = trim(adjustl("GaAs-optcond_")//trim(adjustl(cj))//trim(adjustl(ci))//".dat")
      iarr = [i, j]
      if (rank == 0) open (newunit=out, action="write", file=filename)
      do r = 1, omega_steps
        if (rank == 0) write (unit=out, fmt="(3(1xE15.8))") optcond%cdt(var=1, step=r), &
          real(result(optcond%idims%ind(iarr), r), wp), aimag(result(optcond%idims%ind(iarr), r))
      enddo
      if (rank == 0) close (unit=out)
    enddo
  enddo

  if (rank == 0) write (unit=output_unit, fmt="(A)") "Optical conductivity of GaAs calculated."

  call MPI_FINALIZE(ierror)

contains

  !This function computes the optical conductivity
  !contribution to point k. It suposses 2 integer
  !indices i, j and a continuous index omega.
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

end program
