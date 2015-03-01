	subroutine roeflux(rhol, rhor, ul, ur, pl, pr, f)
		use params_global
		implicit none
		real, intent(in) :: rhol, rhor, ul, ur, pl, pr
		real, dimension(3), intent(out) :: f
		
		real:: el, er, hl, hr, den, uavg, havg, cavg, d1, d2, d3, alpha_1, alpha_2, &
		& alpha_3, lambda_1, lambda_2, lambda_3, f1, f2, f3, cl, cr, da
		
		

		el = pl/(gamma-1) + 0.5*rhol*ul*ul;
		er = pr/(gamma-1) + 0.5*rhor*ur*ur;
		hl = (el + pl)/rhol;
		hr = (er + pr)/rhor;
		cl = sqrt((gamma-1)*(hl - 0.5*ul*ul));
		cr = sqrt((gamma-1)*(hr - 0.5*ur*ur));



		den = sqrt(rhol) + sqrt(rhor);
		uavg = (sqrt(rhol)*ul + sqrt(rhor)*ur)/den;
		havg = (sqrt(rhol)*hl + sqrt(rhor)*hr)/den;
		cavg = sqrt((gamma-1)*(havg - 0.5*uavg*uavg));


		d1 = rhor - rhol;
		d2 = rhor*ur - rhol*ul;
		d3 = er - el;

		alpha_2 = (gamma-1)*((havg - uavg*uavg)*d1 + uavg*d2 - d3)/(cavg*cavg);
		alpha_3 = (d2 + (cavg - uavg)*d1 - cavg*alpha_2)/(2*cavg);
		alpha_1 = d1 - alpha_2 - alpha_3;

		lambda_1 =  abs(uavg - cavg);
		lambda_2 =  abs(uavg);
		lambda_3 =  abs(uavg + cavg);


 		! entropy fix
		da = max(0.0, 4.0*((ur - cr) - (ul - cl)))
		if(lambda_1 < 0.5*da) then
			lambda_1 = lambda_1*lambda_1/da + 0.25*da
		end if


		f1 = lambda_1*alpha_1 + lambda_2*alpha_2 + lambda_3*alpha_3;
		f2 = lambda_1*alpha_1*(uavg-cavg) + lambda_2*alpha_2*uavg + lambda_3*alpha_3*(uavg+cavg);
		f3 = lambda_1*alpha_1*(havg-cavg*uavg) + lambda_2*alpha_2*uavg*uavg/2.0 + lambda_3*alpha_3*(havg+cavg*uavg);
		f(1) = 0.5*((rhol*ul + rhor*ur) - f1);
		f(2) = 0.5*((rhol*ul*ul + pl + rhor*ur*ur + pr) - f2);
		f(3) = 0.5*(ul*hl*rhol + ur*hr*rhor - f3);
	end subroutine roeflux
