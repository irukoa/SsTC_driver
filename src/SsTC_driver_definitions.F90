module SsTC_driver_definitions

  use SsTC_driver_kinds, only: wp => dp

  implicit none

  private

  real(wp), parameter, public :: pi = acos(-1.0_wp)

  complex(wp), parameter, public :: cmplx_0 = cmplx(0.0_wp, 0.0_wp, wp), &
                                    cmplx_1 = cmplx(1.0_wp, 0.0_wp, wp), &
                                    cmplx_i = cmplx(0.0_wp, 1.0_wp, wp)

end module SsTC_driver_definitions
