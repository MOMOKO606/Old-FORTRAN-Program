****************************************************************
*  函数功能：计算sin(kπ/2+x)与cos(kπ/2+x)。
*  当l=0时计算sin(kπ/2+x)，当l=1时计算cos(kπ/2+x)。
*  x为双精度实型变量。
****************************************************************
	double precision function ff(l,k,x)
	double precision x
	if(l.eq.0)then
	  select case(mod(k,4))
	  case(0)
	    ff=sin(x)
	  case(1)
		ff=cos(x)
	  case(2)
		ff=-sin(x)
	  case(3)
		ff=-cos(x)
	  end select
	else
	  select case(mod(k,4))
	  case(0)
	    ff=cos(x)
	  case(1)
		ff=-sin(x)
	  case(2)
		ff=-cos(x)
	  case(3)
		ff=sin(x)
	  end select
	endif
	end function

	  
