	subroutine initialize(nc, xc, q, alpha)
		use params_global
		implicit none
		integer, intent(in) :: nc
		real, dimension(nc), intent(in) :: xc
		real, dimension(nc, 3), intent(out) :: q
		real, intent(in) :: alpha
		integer :: i
		real :: A, dAdx

		do i=1, nc
			call sigma(xc(i), A, dAdx, alpha)
			q(i, 1) = rhoo*A
			q(i, 2) = rhoo*uinf*A
			q(i, 3) = (po/(gamma - 1) + 0.5*rhoo*uinf*uinf)*A
		end do

	end subroutine initialize