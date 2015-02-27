	subroutine sigma(x, A, dAdx, alpha)
		use params_global
		implicit none
		real, intent(in) :: x, alpha
		real, intent(out) :: A, dAdx
		A =  (1.0 - 0.8*x*(1-x))*alpha;
		dAdx = (1.6*x - 0.8)*alpha;
		
	end subroutine sigma