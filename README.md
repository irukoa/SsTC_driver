[![Language](https://img.shields.io/badge/-Fortran-734f96?logo=fortran&logoColor=white)](https://github.com/topics/fortran)
ZENODO PLACEHOLDER
[![Testing suite](https://github.com/irukoa/SsTC_driver/actions/workflows/CI.yml/badge.svg)](https://github.com/irukoa/SsTC_driver/actions/workflows/CI.yml)

# SsTC Driver

### Solid state task constructor driver

This is a modern Fortran library sister to [WannInt](https://github.com/irukoa/WannInt/). This library is meant to provide a backend for automation of sampling and integration workflows of functions defined on the Brillouin zone (BZ) of a crystal.

### Working principle

Many quantities in solid state physics (Berry curvature, optical responses, transport properties...) are expressed as a function, or as the integral of a function, defined in reciprocal space. The most general form of such function is

$$
C^{\alpha}(\textbf{k}; \beta),
$$

where $\textbf{k}$ is a vector in the BZ, $\alpha$ denotes a set of "integer indices", which take only integer values and $\beta$ denotes a set of "continuous variables", which take only real values. The idea behind the library is to automate the task of sampling a generic function $C$ in the BZ and provide the necessary utilities for a user to express the particular function $C$ he/she wishes to sample.

To abstract the idea of integer and continuous indices, [MAC](https://github.com/irukoa/MAC) is used.

# API

The derived type
``` fortran
type, public :: task_specifier
```
is defined.

## `type(task_specifier) :: tsk`
### Constructor.

A task is created by calling

```fortran
call tsk%construct(name, [&
                   int_ind, &
                   cont_data_start, &
                   cont_data_end, &
                   cont_data_steps, &
                   exponent, ]&
                   calculator)
```
where

- `character(len=*), intent(in) :: name` is the name of the task.
- `integer, optional, intent(in) :: int_ind(:)` is the dimension specifier of the integer indices $\alpha$ in $C^{\alpha}(\textbf{k}; \beta)$. That is, the size of the array is the number of integer indices, and the value of each element is the number of values that index can take. If `int_ind` is not present, $C^{\alpha}(\textbf{k}; \beta)$ does not depend on $\alpha$.
- `integer, optional, intent(in) :: cont_data_steps(:)` is the dimension specifier of the continuous indices $\beta$ in $C^{\alpha}(\textbf{k}; \beta)$. That is, the size of the array is the number of continuous variables, and the value of each element is the number of steps that the variable has been discretized into. If `cont_data_steps` is not present, $C^{\alpha}(\textbf{k}; \beta)$ does not depend on $\beta$.
- `real(dp), optional, intent(in) :: cont_data_start(:)` is a real array where each array element specifies the starting sampling point for each of the continuous variables.
- `real(dp), optional, intent(in) :: cont_data_end(:)` is a real array where each array element specifies the ending sampling point for each of the continuous variables.
- `real(dp), optional, intent(in) :: exponent` is a real number $e$, which if not present, will sample the variable $\beta_i$ according to $\beta_i(j)=$ `cont_data_start(i) + (cont_data_end(i) - cont_data_start(i))` $\times (j-1)/$ `(cont_data_steps(i) - 1)`, and if present as $\beta_i(j) \rightarrow e^{\beta_i(j)}$. If `cont_data_steps(i) = 1`, then $\beta_i(j)=$ `cont_data_start(i)`.
- `procedure(cs), pass(self) :: calculator` is a function with interface `cs` which corresponds the implementation of $C^{\alpha}(\textbf{k}; \beta)$. See [interface](#interface).

### Component `type(container_specifier) :: tsk%idims`

Is a [MAC](https://github.com/irukoa/MAC) container specifier initialized to the dimension specifier given by `int_ind`.

### Component `type(container_specifier) :: tsk%cdims`

Is a [MAC](https://github.com/irukoa/MAC) container specifier initialized to the dimension specifier given by `cont_data_steps`.

### Name handle

Is called as,

```fortran
name = tsk%name()
```

where `character(len=120) :: name`.

### Initialization query

Is called as,

```fortran
initialized = tsk%initialized()
```

where `logical :: initialized` is `.true.` if the task has been initialized and `.false.` otherwise.

### Continuous variable retriever

This utility serves to retrieve the value of $\beta_i(j)$. It is called as,

```fortran
beta = tsk%cdt(var, step)
```

where `real(dp) :: beta` has the value $\beta_i(j)$ when `var=`$i$, `step=`$j$.

<a id="interface"></a>
### Interface `cs`

The interface
```fortran
use MAC, only: container_specifier, container
use WannInt, only: crystal
abstract interface
  function cs(self, crys, k, other)
    import dp, task_specifier, crystal
    class(task_specifier), intent(in) :: self
    class(crystal), intent(in) :: crys
    real(dp), intent(in) :: k(3)
    class(*), optional, intent(in) :: other

    complex(dp) :: cs(self%idims%size(), self%cdims%size())
  end function cs
end interface
```
describes the shape of any calculator. The input values for these types of functions are always an object of class `task_specifier`, a [WannInt](https://github.com/irukoa/WannInt/) crystal class object `crys`, a `real(dp), dimension(3)` variable `k` representing a vector in the BZ of the crystal and an optional unlimited polymorphic object `other`. The output value is a `complex(dp)` rank-2 array which holds, in each dimension, the memory layout index corresponding to a particular permutation of the integer and continuous variables, respectively.

#### Best practices

- If you need to pass further external information to the calculator (a value for a Dirac delta function smearing for integration tasks, an exception, or an error case...), it is best to extend either `task_specifier` or `crystal`, depending on what information needs to be passed.
- It is a clever idea to traverse the whole output array with the aid of [MAC](https://github.com/irukoa/MAC)'s indexing when defining calculators.
```fortran
do r = 1, tsk%cdims%size()
  r_arr = tsk%cdims%ind(r)
  cv1 = tsk%cdt(var = 1, step = r_arr(1)) !external variable 1,
  cv2 = tsk%cdt(var = 2, step = r_arr(2)) !external variable 2,
  ![...]
  do i = 1, tsk%idims%size()
    i_arr = tsk%idims%ind(i)
    iv1 = i_arr(1) !integer variable 1,
    iv2 = i_arr(2) !integer variable 2,
    cs(i, r) = ... !your implementation in terms of iv1, cv1...
    ![...]
  enddo
enddo
```
- Most of the times (implementing optical responses, transport properties...), the number of integer and continuous values is known. In these cases it is a clever idea to loop though the variables in a transparent way and then indexing.
```fortran
cvshape = tsk%cdims%shape()
ivshape = tsk%idims%shape()
do w1 = 1, cvshape(1) !external variable 1,
  cv1 = tsk%cdt(var = 1, step = w1)
  do w2 = 1, cvshape(2) !external variable 2,
    cv2 = tsk%cdt(var = 2, step = w2)
    ![...]
    r_arr = [w1, w2, ...]
    do i1 = 1, ivshape(1) !integer variable 1,
      do i2 = 1, ivshape(2) !integer variable 2,
      ![...]
      i_arr = [i1, i2, ...]
      !your implementation in terms of i1, cv1...
      cs(tsk%idims%ind(i_arr), tsk%rdims%ind(r_arr)) = ...
      enddo
    enddo
    ![...]
  enddo
enddo
```

### Sampler

An initialized task is sampled by calling (1st way),
```fortran
call tsk%sample(crys, kpart, store_at [, parallelization, other])
```
where

- `class(crystal), intent(in) :: crys` is a [WannInt](https://github.com/irukoa/WannInt/) crystal.
- `integer, intent(in) :: kpart(3)` is a integer array where each component specifies the number of points to partition each dimension of the BZ into. Each component `kpart(i)` discretizes the reciprocal lattice vector $\mathbf{b}_i$.
- `complex(dp), allocatable, intent(out) :: store_at(:, :, :)/store_at(:, :)` is a rank-3 or rank-2 array. If rank-3:
  - The first dimension specifies the memory layout representation of $\textbf{k}$ in the order that will be laid out by a [MAC](https://github.com/irukoa/MAC) container specifier with `dimension_specifier = kpart`.
  - The second dimension is the memory layout representation of a particular permutation of integer indices $\alpha$ in the order laid out by tsk%idims.
  - The third dimension is the memory layout representation of a particular permutation of continuous variables $\beta$ in the order laid out by tsk%cdims.
- If rank-2:
  - The first dimension is the memory layout representation of a particular permutation of integer indices $\alpha$ in the order laid out by tsk%idims.
  - The second dimension is the memory layout representation of a particular permutation of continuous variables $\beta$ in the order laid out by tsk%cdims.
  - `store_at` holds an unnormalized sum over all the `product(kpart)` points given as input.
- `character(len=*), optional, intent(in) :: parallelization` is an specification of the parallelization method to employ to distribute the sampling. Options are `MPI+OMP`, `MPI`, `OMP`, `none`. Default is `OMP`.
- `class(*), optional, intent(in) :: other` is an unlimited polymorphic object to be passed to the calculator during sampling. See [interface](#interface).

Or by calling (2nd way),
An initialized task is sampled by calling (1st way),
```fortran
call tsk%sample(crys, klist, store_at [, parallelization, other])
```
where all present parameters are the same as described in the previous calling way except for
- `real(dp), intent(in) :: klist(3, nk)`, which is a list of `nk` points where each `klist(:, j)` represents a BZ point $\textbf{k}^j=(k^j_1, k^j_2, k^j_3)$, where each $k^j_i$ is given in coordinates relative to the reciprocal lattice vectors $\mathbf{b}_i$ (crystal coordinates).

## Utilities

The library defines the utilities `kpath`, `kslice`. These utilities are of use to create a list of $\textbf{k}$-points of the form `klist(3, nk)` which traverse a path and a slice of the BZ, respectively.

### Kpath

The utility is used as
```fortran
path = kpath(vecs, nkpts)
```
where,

- `real(dp), intent(in) :: vecs(nv, 3)` are the coordinates of `nk` points where each `vecs(:, j)` represents a point $\textbf{k}^j=(k^j_1, k^j_2, k^j_3)$, where each $k^j_i$ is given in coordinates relative to the reciprocal lattice vectors $\mathbf{b}_i$ (crystal coordinates).
- `integer, intent(in)  :: nkpts(nv - 1)` each element `nkpts(i)` is the number of points in a straight line to consider between points `vecs(i, :)`, `vecs(i + 1, :)`.
- `real(dp), allocatable :: path(:, :)` is the return value. On output, the shape is `[3, nk]`, where `nk = sum(nkpts) - (nvec - 2)`. Each element `path(:, j)` represents a point $\textbf{k}^j=(k^j_1, k^j_2, k^j_3)$, along the ordered path, where each $k^j_i$ is given in coordinates relative to the reciprocal lattice vectors $\mathbf{b}_i$ (crystal coordinates).

### Kslice

The utility is used as
```fortran
slice = kslice(corner, vec_a, vec_b, part)
```
where,

- `real(dp), intent(in) :: corner(3)` are the coordinates of a point $\textbf{k}^j=(k^j_1, k^j_2, k^j_3)$, where each $k^j_i$ is given in coordinates relative to the reciprocal lattice vectors $\mathbf{b}_i$ (crystal coordinates). This parameter represents the bottom-left point of the slice.
- `real(dp), intent(in) :: vec_a(3), vec_b(3)` are the coordinates of vectors $\textbf{k}^1$ and $\textbf{k}^2$, where each $k^j_i$ is given in coordinates relative to the reciprocal lattice vectors $\mathbf{b}_i$ (crystal coordinates), representing two vectors defining the slice.
- `integer, intent(in) :: part(2)` is a integer array where each component specifies the number of points to partition each dimension of the BZ into. Each component `part(i)` discretizes the reciprocal lattice vector $\mathbf{b}_i$.
- `real(dp), allocatable :: slice(:, :)` is the return value. On output, the shape is `[3, nk]`, where `nk = product(part)`. Each element `slice(:, j)` represents a point $\textbf{k}^j=(k^j_1, k^j_2, k^j_3)$, where each $k^j_i$ is given in coordinates relative to the reciprocal lattice vectors $\mathbf{b}_i$ (crystal coordinates). The ordering if $\textbf{k}$ points is given by the order that will be laid out by a [MAC](https://github.com/irukoa/MAC) container specifier with `dimension_specifier = part`.

# Build

An automated build is available for [Fortran Package Manager](https://fpm.fortran-lang.org/) users. This is the recommended way to build and use WannInt in your projects. You can add WannInt to your project dependencies by including

```
[dependencies]
SsTC_driver = { git="https://github.com/irukoa/SsTC_driver.git" }
```
to the `fpm.toml` file.

[MAC](https://github.com/irukoa/MAC)'s objects
``` fortran
type, public :: container_specifier
type, extends(container_specifier), public :: container
```
and [WannInt](https://github.com/irukoa/WannInt/)'s objects and utilities
``` fortran
type, public :: crystal
public :: diagonalize
public :: dirac_delta
public :: deg_list
public :: schur
public :: SVD
public :: expsh
public :: logu
```
are made public by SsTC_driver.

# Limitations

The mayor limitation lies in memory management when sampling with parallelization. The array `store_at` needs to be dynamically allocated. Its size is given by $($ `nk` $) \times$ `product(int_ind)` $\times$ `product(cont_data_steps)`. Before making any calculation it is suggested to check the value of its size and choose a parallelization scheme accordingly. It is recommended to use `none` or `OMP` parallelization in PCs and `MPI` or `MPI+OMP` in multinode clusters. If the size of `store_at` is too big, [segmentation faults](https://stackoverflow.com/a/19797710/22403953) may occur. We recommend using only `store_at` in its rank-3 version in sampling calls if peeping the index of $\textbf{k}$ points is important (for plotting or for special integration schemes).

We recommend the compilers in the [Intel oneAPI](https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit.html) toolkit `ifort/mpiifort` and `ifx/mpiifx` or the GNU compilers `gfortran/mpifort`. If possible, we recommend using the `--heap-arrays` flag for Intel compilers and the `-fmax-stack-var-size=n` flag for GNU compilers.
