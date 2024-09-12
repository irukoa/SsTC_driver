module SsTC_driver_utils

  use SsTC_driver_kinds, only: wp => dp
  use MAC, only: container_specifier
  use WannInt, only: crystal, SVD

  implicit none

  private

  public :: kpath, kslice
  public :: cart_to_crys, crys_to_cart

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

  function crys_to_cart(k_crys, crys) result(k_cart)
    real(wp), intent(in) :: k_crys(3)
    class(crystal), intent(in) :: crys
    real(wp) :: k_cart(3)

    k_cart = matmul(transpose(crys%reciprocal_lattice_basis()), k_crys)

  end function crys_to_cart

  function cart_to_crys(k_cart, crys) result(k_crys)
    real(wp), intent(in) :: k_cart(3)
    class(crystal), intent(in) :: crys
    real(wp) :: k_crys(3)

    integer :: i
    complex(wp) :: recip_latt_basis(3, 3), &
                   U(3, 3), V(3, 3), SGM(3, 3)
    real(wp) :: inv_recip_latt_basis(3, 3)

    recip_latt_basis = cmplx(transpose(crys%reciprocal_lattice_basis()), 0.0_wp, wp)

    call SVD(matrix=recip_latt_basis, U=U, V=V, sigma=sgm)
    do i = 1, 3
      SGM(i, i) = cmplx(1.0_wp/real(SGM(i, i), wp), 0.0_wp, wp)
    enddo
    inv_recip_latt_basis = real(matmul(matmul(V, SGM), transpose(conjg(U))), wp)

    k_crys = matmul(inv_recip_latt_basis, k_cart)

  end function cart_to_crys

end module SsTC_driver_utils
