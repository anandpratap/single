	subroutine metrics(nc, x, xc, dxc, l)
		use params_global
		implicit none
		integer :: nc
		real, dimension(nc+1), intent(out) :: x
		real, dimension(nc), intent(out) :: xc, dxc
		real, intent(in) :: l
		integer :: i
		
		do i=1, nc+1
			x(i) = (i-1)*l/(nc)
		end do

		do i=1, nc
			xc(i) = (x(i+1) + x(i))/2.0
			dxc(i) = (x(i+1) - x(i))
		end do

	end subroutine metrics