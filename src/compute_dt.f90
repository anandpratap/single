	subroutine compute_dt(nc, xc, dxc, q, dt, alpha)
		use params_global
		implicit none
		integer, intent(in) :: nc
		real, dimension(nc), intent(in) :: xc, dxc
		real, dimension(nc, 3), intent(in) :: q
		real, intent(out) :: dt
		real, intent(in) :: alpha
		real, dimension(nc) :: rho, u, p, c, lam
		

		call primvars(nc, xc, q, rho, u, p, alpha)
		c = sqrt(gamma*p/rho);
		lam = c + abs(u);
		dt = minval(dxc(:))/maxval(lam(:))*cfl;

	end subroutine compute_dt