subroutine write_data(nc, xc, q, alpha)
	use params_global
	implicit none
	integer, intent(in) :: nc
	real, dimension(nc), intent(in) :: xc
	real, dimension(nc, 3), intent(in) :: q
	real, intent(in) :: alpha
 
	real, dimension(nc) :: rho, u, p
	integer :: i

	call primvars(nc, xc, q, rho, u, p, alpha)
	open(unit=log_data, file='data.dat', status='unknown', form='formatted')
	write(log_data, *) "TITLE = ""Single Solution FIle"""
	write(log_data, *) "VARIABLES = ""X"" ""rho"" ""u"" ""p"" ""M"""
	write(log_data, *) "ZONE T=""Only Zone"", I = ", nc, ", F = Point"
	do i = 1, nc
		write(log_data, *) xc(i), rho(i), u(i), p(i), u(i)/sqrt(gamma*p(i)/rho(i))
	end do
	close(log_data)
          
end subroutine write_data