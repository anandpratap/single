	subroutine sensitivity(nc, x, xc, dxc, q, alpha, psi)
		use params_global
		implicit none
		integer, intent(in) :: nc
		real, dimension(nc), intent(in) :: xc, dxc
		real, dimension(nc+1), intent(in) :: x
		real, dimension(3*nc), intent(in) :: psi
		real, dimension(nc, 3), intent(inout) :: q
		real, intent(in) :: alpha

		real :: alphab, obj, objb
		real, dimension(nc, 3) :: res
		real, dimension(nc, 3) :: qb, q1
		real, dimension(3*nc) :: dRdalpha

		
		objb = 1.0
		qb(:, :) = 0.0
		alphab = 1.0
		q1(:, :) = q(:, :)
		call calc_residual_alpha(nc, x, xc, dxc, q, res, dRdalpha, alpha, alphab)
		q(:, :) = q1(:, :)

		objb = 1.0
		qb(:,:) = 0.0
		alphab = 1.0
		call objective(nc, xc, q, obj, alpha)
		call objective_bq(nc, xc, q, qb, obj, objb, alpha, alphab)
		write(11, *) obj
		write(11, *) alphab + sum(psi*dRdalpha)
		print *,'dJdAlpha: ', alphab + sum(psi*dRdalpha)

	end subroutine sensitivity