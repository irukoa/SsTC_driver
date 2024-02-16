program Si_Band_Structure
  use, intrinsic :: iso_fortran_env, only: output_unit
  use MPI
  use OMP_LIB

  use SsTC_driver_kinds, only: wp => dp
  use SsTC_driver_utils, only: kpath
  use SsTC_driver, only: task_specifier
  use WannInt, only: crystal, diagonalize

  !In this example we calculate the band structure of Si in the
  !L - G - X - K - G path.

  implicit none

  real(wp), allocatable :: path(:, :)
  real(wp) :: coords(5, 3)

  type(crystal) :: Si
  type(task_specifier) :: band_structure
  complex(wp), allocatable :: result(:, :, :)

  integer :: out
  integer :: ik, ibnd

  integer :: ierror, rank

  !MPI things.
  call MPI_INIT(ierror)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)

  coords(1, :) = [0.5_wp, 0.5_wp, 0.5_wp] !L
  coords(2, :) = [0.0_wp, 0.0_wp, 0.0_wp] !G
  coords(3, :) = [0.5_wp, 0.0_wp, 0.5_wp] !X
  coords(4, :) = [0.375_wp, -0.375_wp, 0.0_wp] !K
  coords(5, :) = [0.0_wp, 0.0_wp, 0.0_wp] !G

  !This will create a kpoint dictionary traversing the path.
  !nkpts determines the number of points between each pair
  !of points.
  path = kpath(vecs=coords, nkpts=[100, 100, 100, 100])

  !Load data for Silicon from file using WannInt.
  call Si%construct(name="Si", &
                    from_file="./material_data/Si_tb.dat", &
                    fermi_energy=6.3869_wp)

  !Here, we define the task to sample. It has a single integer-like
  !index ranging from 1 to Si%num_bands() (which in this case is crystal dependent)
  !and the way to calculate it is found in the function band_structure_calculator.
  !Take a look to the function in the contains section.
  call band_structure%construct(name="band_structure", &
                                int_ind=[Si%num_bands()], &
                                calculator=band_structure_calculator)

  !Here we do the actual sampling. The program is told to load electronic structure
  !data and sample along the given path. It will store the results in the array 'results'.
  !The sampling will involve a hybrid MPI+OpenMP parallelization scheme.
  call band_structure%sample(crys=Si, klist=path, store_at=result, &
                             parallelization="MPI+OMP")

  !Lastly, we print to files.
  !You can visualize the bands using gnuplot: gnuplot> p 'Si-bands.dat' u 1:5 w l
  if (rank == 0) open (newunit=out, action="write", file=trim(Si%name())//"-bands.dat")
  do ibnd = 1, Si%num_bands()
    do ik = 1, size(path(1, :))
      !Notice, data is always stored in the next order: kpt index, integer-like index and
      !continuous-like index.
      if (rank == 0) write (unit=out, fmt="(i0, 1x, 4(1xE15.8))") ik, path(:, ik), real(result(ik, ibnd, 1), wp)
    enddo
    if (rank == 0) write (unit=out, fmt=*) ""
  enddo
  if (rank == 0) close (unit=out)

  if (rank == 0) write (unit=output_unit, fmt="(A)") "Band structure of Si calculated."

  !You may have noticed that this workflow can be generalized. By using SsTC_driver and WannInt, you can
  !create an automatic band structure sampler which takes as input any crystal, a set of k-point coordinates,
  !and an array specifying the number of sampling points between pairs of k-point coordinates.
  !Go and try it out!

  call MPI_FINALIZE(ierror)

  deallocate (path, result)

contains

  !This function passes the energy eigenvales of a crystal as
  !integer indices of a task.
  function band_structure_calculator(self, crys, k, other)
    class(task_specifier), intent(in) :: self
    class(crystal), intent(in) :: crys
    real(wp), intent(in) :: k(3)
    class(*), optional, intent(in) :: other

    complex(wp) :: band_structure_calculator(self%idims%size(), self%cdims%size())

    complex(wp) :: H(crys%num_bands(), crys%num_bands()), &
                   P(crys%num_bands(), crys%num_bands())
    real(wp) :: eig(crys%num_bands())

    integer :: i, r

    !Get Hamiltonian.
    H = crys%hamiltonian(kpt=k)
    !Get eigenvalues.
    call diagonalize(matrix=H, P=P, eig=eig)

    !Set to result.
    do i = 1, self%idims%size()
      do r = 1, self%cdims%size()
        band_structure_calculator(i, r) = &
          cmplx(eig(i), 0.0_wp, wp)
      enddo
    enddo

  end function band_structure_calculator

end program Si_Band_Structure
