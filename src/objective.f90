	subroutine objective(nc, xc, q, obj, alpha)
		use params_global
		implicit none
		integer, intent(in) :: nc
		real, dimension(nc), intent(in) :: xc
		real, dimension(nc, 3), intent(in) :: q
		real, intent(out) :: obj
		real, intent(in) :: alpha
		real, dimension(nc) :: rho, u, p
		integer :: i
		call primvars(nc, xc, q, rho, u, p, alpha)

		obj = 0.0

		do i=1,nc
			obj = obj + p(i)**2
		end do
	end subroutine objective