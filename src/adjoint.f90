	subroutine adjoint(nc, x , xc, dxc, q, alpha, psi)
		use params_global
		implicit none
		integer, intent(in) :: nc
		real, dimension(nc), intent(in) :: xc, dxc
		real, dimension(nc+1), intent(in) :: x
		real, dimension(nc, 3), intent(inout) :: q
		real, dimension(3*nc), intent(out) :: psi
		real, intent(in) :: alpha
			
		real :: alphab, obj, objb
		real, dimension(3*nc, 3*nc) :: jac
		real, dimension(nc, 3) :: q1, qb
		real, dimension(nc, 3) :: res
		real, dimension(3*nc) :: dJdq
		integer, dimension(3*nc) :: pivot
		integer :: i, j, ok

		alphab = 0.0
		call calc_residual(nc, x, xc, dxc, q, res, alpha)
		call calc_residual_jacobian(nc, x, xc, dxc, q, res, jac, alpha, alphab)

		jac = transpose(jac)

		call objective(nc, xc, q, obj, alpha)

		alphab = 0.0
		objb = 1.0
		q1(:,:) = q(:,:)
		call objective_bq(nc, xc, q, qb, obj, objb, alpha, alphab)
		q(:,:) = q1(:,:)


		! RHS calculation
		do j=1, 3
			do i=1, nc
				dJdq(3*(i-1) +j) = -qb(i, j)
			end do
		end do
		    
		! solve system of equation
		call dgesv(3*nc, 1, jac, 3*nc, pivot, dJdq, 3*nc, ok)
		print *, 'Adjoint solver done:', ok
		print *, 'Objective function value: ', obj
		psi(:) = dJdq(:)

	end subroutine adjoint