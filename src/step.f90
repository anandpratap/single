	subroutine step(nc, x, xc, dxc, q, dt, l2norm)
		use params_global
		implicit none
		integer, intent(in) :: nc
		real, dimension(nc+1), intent(in) :: x
		real, dimension(nc), intent(in) :: xc, dxc
		real, dimension(nc, 3), intent(inout) :: q
		real, intent(out) :: dt, l2norm

		real, dimension(nc, 3) :: res, dq
		real, dimension(nc*3, nc*3) :: jac
		real, dimension(nc*3) :: b
		real :: alpha, alphab
		integer :: i, j, ii, jj, pivot(3*nc), ok, i_idx, j_idx
		
		alpha = 1.0
		alphab = 0.0

		call compute_dt(nc, xc, dxc, q, dt, alpha)
		call calc_residual(nc, x, xc, dxc, q, res, alpha)
		call calc_residual_jacobian(nc, x, xc, dxc, q, res, jac, alpha, alphab)


		! regularization
		do i = 1, 3*nc
			jac(i, i) = 1/dt + jac(i, i)
		end do
		
		! RHS calculation
		do j=1, 3
			do i=1, nc
				b(3*(i-1) +j) = res(i, j)
			end do
		end do
		
		! solve system of equation
		call dgesv(3*nc, 1, jac, 3*nc, pivot, b, 3*nc, ok)

		! delta q
		do j=1, 3
			do i=1, nc
				dq(i,j) = b(3*(i-1) +j)
			end do
		end do

		! increment
		q(:, :) = q(:, :) + dq(:, :)

		l2norm = 0.0
		do i=1, nc
			l2norm = l2norm + res(i,1)**2
		!	print *, i, l2norm, res(i,1)
		end do
		l2norm = sqrt(l2norm)/nc
	end subroutine step