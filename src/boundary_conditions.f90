	subroutine boundary_conditions(nc, idx, xc, q, rho, u, p)
		use params_global
		implicit none
		real, intent(in) :: xc
		real, dimension(3), intent(out) :: rho, u, p
		real, dimension(3, 3), intent(in) :: q
		integer, intent(in) :: nc
		
		real :: To, cp, rat, T
		
		call primvars_(nc, xc, q, rho, u, p)
		
		u(1) = u(2)
		To = po/rhoo/R
		cp = gamma*R/(gamma-1)
		T = To - u(1)*u(1)/2.0/cp
		rat = T/To
		rho(1) = rhoo*rat**(1/(gamma-1))
		p(1) = po*rat**(gamma/(gamma-1))
	
		rho(nc+2) = rho(nc+1);
		u(nc+2) = u(nc+1);
		p(nc+2) = pout + 0.5*rho(nc+2)*u(nc+2)**2;

	end subroutine boundary_conditions