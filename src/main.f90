	program single
		use omp_lib
		use params_global
		implicit none

		integer :: n, nc, i, j, counter
		real:: l, t, dt, l2norm, alpha
		real, dimension(:), allocatable :: x, xc, dxc, psi
		real, dimension(:, :), allocatable :: q

		alpha = 1.0+1e-6
		n = 1001
		nc = n - 1
		l = 1.0
		allocate(x(n), xc(nc), dxc(nc))
		allocate(q(nc,3))

		call metrics(nc, x, xc, dxc, l)
		call initialize(nc, xc, q, alpha)

		counter = 0
		do
			call step(nc, x, xc, dxc, q, dt, l2norm)

			if(mod(counter, 1) .eq. 0) then
				print *, 'step: ', counter, 'l2norm: ', l2norm, 'dt', dt
				call write_data(nc, xc, q, alpha)
			end if
			
			if(l2norm .le. 1e-6) then
				call write_data(nc, xc, q, alpha)
				exit
			end if
			counter = counter + 1
		end do

		if(ifadjoint) then
			! calculate adjoint
			allocate(psi(3*nc))
			call adjoint(nc, x , xc, dxc, q, alpha, psi)
			! calculate sensitivity
			call sensitivity(nc, x, xc, dxc, q, alpha, psi)
			deallocate(psi)
		end if

		deallocate(q)
		deallocate(x, xc, dxc)
	end program single