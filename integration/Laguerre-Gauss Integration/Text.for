*************************************************************
*  拉盖尔-高斯求积法子程序测试程序
*************************************************************	
	external f
	double precision a,f,r
	a=0.0
	call LGI(f,r)
	write(*,*) r
	end
	
	subroutine LGI(f,r)
	dimension t(6),coef(6)
	double precision f,r,x
	data t/0.2228466041,1.1889321016,2.9927363260,5.7751435691,
     &       9.8374674183,15.9828739806/
	data coef/0.4589646793,0.4170008307,0.1133733820,0.0103991975,
     &          0.0002610172,0.0000008985/
	r=0.0
	do i=1,6
	  x=t(i)
	  r=r+coef(i)*exp(x)*f(x)
	enddo
	end subroutine

	function f(x)
	double precision x
	f=x*exp(-x)
	end function