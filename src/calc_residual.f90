	subroutine calc_residual(nc, x, xc, dxc, q, res, alpha)
		use params_global
		implicit none
		integer, intent(in) :: nc
		real, intent(in) :: alpha
		real, dimension(nc+1), intent(in) :: x
		real, dimension(nc), intent(in) :: xc, dxc
		real, dimension(nc, 3), intent(in) :: q
		real, dimension(nc, 3), intent(out) :: res
		
		integer :: i
		res(:,:) = 0.0

		do i=1, nc
			if (i .eq. 1) then
				call residual(nc, i, x(i:i+1), xc(1:3), dxc(1:3), q(1:3,:), res(i,:), alpha)
			else if(i .eq. nc) then
				call residual(nc, i, x(i:i+1), xc(nc-2:nc), dxc(nc-2:nc), q(nc-2:nc,:), res(i,:), alpha)
			else 
				call residual(nc, i, x(i:i+1), xc(i-1:i+1), dxc(i-1:i+1), q(i-1:i+1,:), res(i,:), alpha)
			end if
		end do
		

	end subroutine calc_residual

	subroutine calc_residual_jacobian(nc, x, xc, dxc, q, res, jac, alpha, alphab)
		use params_global
		implicit none
		integer, intent(in) :: nc
		real, dimension(nc+1), intent(in) :: x
		real, dimension(nc), intent(in) :: xc, dxc
		real, dimension(nc, 3), intent(inout) :: q
		real, dimension(nc, 3), intent(inout) :: res
		real, dimension(3*nc, 3*nc), intent(out) :: jac
		real, intent(inout) :: alpha, alphab
		real, dimension(nc, 3) :: q1, res1
		real, dimension(3, 3) :: qb
		real, dimension(1, 3) :: resb
		integer :: i, j, ii, jj, i_idx, j_idx

		real :: a1, a1b
		resb(:, :) = 0.0
		jac(:,:) = 0.0
		
		res1(:,:) = res(:,:)
		q1(:,:) = q(:,:)

		do i=1, nc
			do j=1, 3
				resb(1,j) = 1.0
				if (i .eq. 1) then
					res1(:,:) = res(:,:)
					q1(:,:) = q(:,:)
					a1 = alpha
					a1b = alphab
					call residual_bq(nc, i, x(i:i+1), xc(1:3), dxc(1:3), q(1:3,:), qb(:,:), res(i,:), resb(1,:), alpha, alphab)
					res(:,:) = res1(:,:)
					q(:,:) = q1(:,:)
					alpha = a1
					alphab = a1b
					do ii = i-1, i
						do jj = 1, 3
							i_idx = 3*(i-1) + j
							j_idx = 3*(ii-1+1) + jj
							jac(i_idx, j_idx) = -qb(ii-i+2, jj) + jac(i_idx, j_idx)
						end do
					end do


				else if(i .eq. nc) then
					res1(:,:) = res(:,:)
					q1(:,:) = q(:,:)
					a1 = alpha
					a1b = alphab
					
					call residual_bq(nc, i, x(i:i+1), xc(nc-2:nc), dxc(nc-2:nc), q(nc-2:nc,:),  qb(:,:) , res(i,:), resb(1,:), alpha, alphab)
					res(:,:) = res1(:,:)
					q(:,:) = q1(:,:)
					alpha = a1
					alphab = a1b
					
					do ii = i, i+1
						do jj = 1, 3
							i_idx = 3*(i-1) + j
							j_idx = 3*(ii-1-1) + jj
							jac(i_idx, j_idx) = -qb(ii-i+2, jj) + jac(i_idx, j_idx)
						end do
					end do


				else 
					res1(:,:) = res(:,:)
					q1(:,:) = q(:,:)
					a1 = alpha
					a1b = alphab
					
					call residual_bq(nc, i, x(i:i+1), xc(i-1:i+1), dxc(i-1:i+1), q(i-1:i+1,:), qb(:,:), res(i,:), resb(1,:), alpha, alphab)
					res(:,:) = res1(:,:)
					q(:,:) = q1(:,:)
					alpha = a1
					alphab = a1b
					
					do ii = i-1, i+1
						do jj = 1, 3
							i_idx = 3*(i-1) + j
							j_idx = 3*(ii-1) + jj
							jac(i_idx, j_idx) = -qb(ii-i+2, jj) + jac(i_idx, j_idx)
						end do
					end do
				end if
				resb(:,:) = 0.0
				qb(:,:) = 0.0
			end do
		end do


	end subroutine calc_residual_jacobian


	subroutine calc_residual_alpha(nc, x, xc, dxc, q, res, jac, alpha, alphab)
		use params_global
		implicit none
		integer, intent(in) :: nc
		real, dimension(nc+1), intent(in) :: x
		real, dimension(nc), intent(in) :: xc, dxc
		real, dimension(nc, 3), intent(inout) :: q
		real, dimension(nc, 3), intent(inout) :: res
		real, dimension(3*nc), intent(out) :: jac
		real, intent(inout) :: alpha, alphab
		real, dimension(nc, 3) :: q1, res1
		real, dimension(3, 3) :: qb
		real, dimension(1, 3) :: resb
		integer :: i, j, ii, jj, i_idx, j_idx

		real :: a1, a1b
		resb(:, :) = 0.0
		jac(:) = 0.0
		
		res1(:,:) = res(:,:)
		q1(:,:) = q(:,:)

		do i=1, nc
			do j=1, 3
				resb(1,j) = 1.0
				if (i .eq. 1) then
					res1(:,:) = res(:,:)
					q1(:,:) = q(:,:)
					a1 = alpha
					a1b = alphab
					qb(:,:) = 0.0
					call residual_bq(nc, i, x(i:i+1), xc(1:3), dxc(1:3), q(1:3,:), qb(:,:), res(i,:), resb(1,:), alpha, alphab)
					jac(3*(i-1) + j) = alphab
					res(:,:) = res1(:,:)
					q(:,:) = q1(:,:)

					alpha = a1
					alphab = a1b
				


				else if(i .eq. nc) then
					res1(:,:) = res(:,:)
					q1(:,:) = q(:,:)
					a1 = alpha
					a1b = alphab
					qb(:,:) = 0.0
					call residual_bq(nc, i, x(i:i+1), xc(nc-2:nc), dxc(nc-2:nc), q(nc-2:nc,:),  qb(:,:) , res(i,:), resb(1,:), alpha, alphab)
					jac(3*(i-1) + j) = alphab
					res(:,:) = res1(:,:)
					q(:,:) = q1(:,:)
					alpha = a1
					alphab = a1b
					
				

				else 
					res1(:,:) = res(:,:)
					q1(:,:) = q(:,:)
					a1 = alpha
					a1b = alphab
					qb(:,:) = 0.0
					call residual_bq(nc, i, x(i:i+1), xc(i-1:i+1), dxc(i-1:i+1), q(i-1:i+1,:), qb(:,:), res(i,:), resb(1,:), alpha, alphab)
					jac(3*(i-1) + j) = alphab
					res(:,:) = res1(:,:)
					q(:,:) = q1(:,:)
					alpha = a1
					alphab = a1b
					
					
				end if
				resb(:,:) = 0.0
				qb(:,:) = 0.0
			end do
		end do


	end subroutine calc_residual_alpha