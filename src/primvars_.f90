	subroutine primvars_(nc, xc, q, rho, u, p, alpha)
		use params_global
		implicit none
		real, dimension(nc), intent(in) :: xc
		real, dimension(nc, 3), intent(in) :: q
		real, dimension(nc+2), intent(out) :: rho, u, p
		real, intent(in) :: alpha
		integer, intent(in) :: nc
		real :: A, dAdx
		integer :: i
		do i=2, nc+1
			call sigma(xc(i-1), A, dAdx, alpha)
			rho(i) = q(i-1,1)/A;
			u(i) = q(i-1,2)/rho(i)/A;
			p(i) = (q(i-1,3)/A - 0.5*rho(i)*u(i)*u(i))*(gamma-1);
		end do

	end subroutine primvars_