module parameters

integer, parameter :: dp = selected_real_kind(15, 307)
integer, parameter :: qp = selected_real_kind(33, 4931)
integer, parameter :: i4 = selected_int_kind(9)
integer, parameter :: i8 = selected_int_kind(18)

complex( kind = dp ), parameter :: re = (1.0_dp, 0.0_dp)
complex( kind = dp ), parameter :: im = (0.0_dp, 1.0_dp)

integer :: stepsForJudge = 1000
logical :: quietON = .true.
real( kind = dp ) :: convCriterion = 1.0d-08
real( kind = dp ) :: chopCutoff = 1.0d-20
real( kind = dp ) :: dt = 0.001_dp

end module
