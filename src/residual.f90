	subroutine residual(nc, idx, x, xc, dxc, q, res, alpha)
		use params_global
		implicit none
		integer, intent(in) :: nc, idx
		real, dimension(2), intent(in) :: x
		real, dimension(3), intent(in) :: xc, dxc
		real, dimension(3, 3), intent(in) :: q
		real, dimension(1, 3), intent(out) :: res
		real, intent(in) :: alpha
		real, dimension(3) :: rho, u, p
		real, dimension(2) :: rhol, rhor, ul, ur, pl, pr
		real, dimension(3) :: fl, fr
		integer :: i, rdx
		real :: A, dAdx
		real :: To, cp, rat, T

		if(idx .eq. 1) then
			u(1) = u(2)
			To = po/rhoo/R
			cp = gamma*R/(gamma-1)
			T = To - u(1)*u(1)/2.0/cp
			rat = T/To
			rho(1) = rhoo*rat**(1/(gamma-1))
			p(1) = po*rat**(gamma/(gamma-1))

			do i=1, 2
				rdx = i
				call sigma(xc(rdx), A, dAdx, alpha)
				rho(rdx+1) = q(rdx,1)/A;
				u(rdx+1) = q(rdx,2)/rho(rdx+1)/A;
				p(rdx+1) = (q(rdx,3)/A - 0.5*rho(rdx+1)*u(rdx+1)*u(rdx+1))*(gamma-1);
			end do
		else if(idx .eq. nc) then
			rho(3) = rho(2);
			u(3) = u(2);
			p(3) = pout + 0.5*rho(3)*u(3)**2;

			do i=1, 2
				rdx = i+1
				call sigma(xc(rdx), A, dAdx, alpha)
				rho(i) = q(rdx,1)/A;
				u(i) = q(rdx,2)/rho(rdx-1)/A;
				p(i) = (q(rdx,3)/A - 0.5*rho(rdx-1)*u(rdx-1)*u(rdx-1))*(gamma-1);
			end do
		else
			do i=1, 3
				call sigma(xc(i), A, dAdx, alpha)
				rho(i) = q(i,1)/A;
				u(i) = q(i,2)/rho(i)/A;
				p(i) = (q(i,3)/A - 0.5*rho(i)*u(i)*u(i))*(gamma-1);
			end do
		end if

		rhol(1) = rho(1)
		rhor(1) = rho(2)
		ul(1) = u(1)
		ur(1) = u(2)
		pl(1) = p(1)
		pr(1) = p(2)

		rhol(2) = rho(2)
		rhor(2) = rho(3)
		ul(2) = u(2)
		ur(2) = u(3)
		pl(2) = p(2)
		pr(2) = p(3)

		call roeflux(rhol(1), rhor(1), ul(1), ur(1), pl(1), pr(1), fl(:))
		call roeflux(rhol(2), rhor(2), ul(2), ur(2), pl(2), pr(2), fr(:))
		call sigma(xc(2), A, dAdx, alpha)

		res(1,:) =  -(fr(:) - fl(:))/dxc(2)*A
		res(1,1) = res(1,1) - rho(2)*u(2)*(dAdx);
		res(1,2) = res(1,2) - rho(2)*u(2)*u(2)*(dAdx);
		res(1,3) = res(1,3) - (rho(2)*u(2)*u(2)/2.0 + p(2)*gamma/(gamma-1))*u(2)*(dAdx);

	end subroutine residual