subroutine reconstruct(nc, rho, u, p, rhol, rhor, ul, ur, pl, pr)
	use params_global
	implicit none
	real, dimension(nc+2), intent(in) :: rho, u, p
	real, dimension(nc+1), intent(out) :: rhol, rhor, ul, ur, pl, pr
	integer, intent(in) :: nc

	integer :: n
	integer :: i
	n = nc + 1

	do i=1, nc+1
		rhol(i) = rho(i);
		rhor(i) = rho(i+1);

		ul(i) = u(i);
		ur(i) = u(i+1);

		pl(i) = p(i);
		pr(i) = p(i+1);
	end do
end subroutine reconstruct