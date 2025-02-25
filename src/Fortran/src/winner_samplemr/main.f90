program generate_gaussian
use mkl_vsl
implicit none

type(VSL_STREAM_STATE) :: stream
integer :: status, error_code
  integer, parameter :: n = 1000
  real(8), dimension(n) :: r
  real(8) :: mean, sigma

  ! Set the parameters for the Gaussian distribution
  mean  = 0.0d0
  sigma = 1.0d0

! Create a new random stream using the Mersenne Twister generator with seed 777
error_code = vslNewStream(stream, VSL_BRNG_MT19937, 777)
if (error_code /= 0) then
print *, "Error initializing random stream: ", error_code
stop
endif

! Generate 'n' random numbers with a Gaussian distribution (Box-Muller method)
error_code = vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER, stream, n, r, mean, sigma)
if (error_code /= 0) then
print *, "Error generating random numbers: ", error_code
stop
endif

  ! Example: print the first 10 generated random numbers
  print *, r(1:10)

! Delete the random stream to free resources
error_code = vslDeleteStream(stream)
if (error_code /= 0) then
print *, "Error deleting stream: ", error_code
endif

end program generate_gaussian

