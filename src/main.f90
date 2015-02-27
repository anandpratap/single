	program single
		use params_global
		implicit none

		integer :: n, nc
		real:: l
		real, dimension(:), allocatable :: x
		real, dimension(:), allocatable :: xc, dxc, b, dRdalpha, psi
		integer, dimension(:), allocatable :: pivot
		real, dimension(:, :), allocatable :: q, jac, jact, res, qb, q1

		real :: t, dt, l2norm, A, dAdx, obj, objb
		integer :: i, j, counter, ok
		real :: alpha, alphab
		
		alpha = 1.0+1e-6
		n = 1001
		nc = n - 1
		l = 1.0
		allocate(x(n), xc(nc), dxc(nc))
		allocate(q(nc,3))


		do i=1, n
			x(i) = (i-1)*l/(n-1)
		end do

		do i=1, nc
			xc(i) = (x(i+1) + x(i))/2.0
			dxc(i) = (x(i+1) - x(i))
		end do
		do i=1, nc
			call sigma(xc(i), A, dAdx, alpha)
			q(i, 1) = rhoo*A
			q(i, 2) = rhoo*uinf*A
			q(i, 3) = (po/(gamma - 1) + 0.5*rhoo*uinf*uinf)*A
		end do

		counter = 0
		do 
			call step(nc, x, xc, dxc, q, dt, l2norm)
			counter = counter + 1
			if(mod(counter, 1) .eq. 0) then
				print *, 'step: ', counter, 'l2norm: ', l2norm, 'dt', dt
				call write_data(nc, xc, q)
			end if
			
			if(l2norm .le. 1e-6) then
				call write_data(nc, xc, q)
				exit
			end if
			!stop;
		end do

		allocate(jac(3*nc, 3*nc), jact(3*nc, 3*nc), res(nc, 3), qb(nc, 3), b(3*nc), pivot(3*nc), psi(3*nc), dRdalpha(3*nc), q1(nc,3))
		
		alphab = 0.0
		call calc_residual(nc, x, xc, dxc, q, res, alpha)

		call calc_residual_jacobian(nc, x, xc, dxc, q, res, jac, alpha, alphab)
		
		do i=1,3*nc
			do j=1,3*nc
				jact(i, j) = jac(j, i)
			end do
		end do
	

		call calc_residual_alpha(nc, x, xc, dxc, q, res, dRdalpha, alpha, alphab)
		
		call objective(nc, xc, q, obj, alpha)
		
		alphab = 0.0
		objb = 1.0
		
		q1(:,:) = q(:,:)
		call objective_bq(nc, xc, q, qb, obj, objb, alpha, alphab)
		print *, alphab
		q(:,:) = q1(:,:)


		! RHS calculation
		do i=1, nc
			do j=1, 3
				b(3*(i-1) +j) = -qb(i, j)
			end do
		end do
		! solve system of equation
		
		call dgesv(3*nc, 1, jact, 3*nc, pivot, b, 3*nc, ok)
		print *, 'normb', norm2(b)
		print *, 'ok:', ok
		print *, 'Objective function: ', obj
		psi(:) = b(:)

		objb = 1.0
		qb(:,:) = 0.0
		alphab = 1.0
		alpha = 1.0
		call objective(nc, xc, q, obj, alpha)
		call objective_bq(nc, xc, q, qb, obj, objb, alpha, alphab)
		print *, alphab, objb
		print *, alphab + sum(psi*dRdalpha)
		deallocate(dRdalpha, psi, pivot, b, qb, res, jac)
		deallocate(q, q1)
		deallocate(x, xc, dxc)
	end program single