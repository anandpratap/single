	subroutine primvars(nc, xc, q, rho, u, p, alpha)
		use params_global
		implicit none
		real, dimension(nc), intent(in) :: xc
		real, dimension(nc), intent(out) :: rho, u, p
		real, dimension(nc, 3), intent(in) :: q
		integer, intent(in) :: nc
		real, intent(in) :: alpha

		real :: A, dAdx
		
		integer :: i
		do i=1, nc
			call sigma(xc(i), A, dAdx, alpha)
			rho(i) = q(i,1)/A;
			u(i) = q(i,2)/rho(i)/A;
			p(i) = (q(i,3)/A - 0.5*rho(i)*u(i)*u(i))*(gamma-1);
		end do

	end subroutine primvars