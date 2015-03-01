	module params_global
		real :: gamma, cfl
		real :: rhoinf, uinf, pinf, rhoo, po, R, pout
		
		data rhoinf, uinf, pinf/1.0, 0.0, 1.0/
		data rhoo, po, pout/1.0, 1.0, 0.7/
		data gamma, cfl, R/1.4, 1000.0, 287.0/
		
		integer :: log_data
		data log_data/10/
	
		logical :: ifadjoint
		data ifadjoint/.true./

	end module params_global