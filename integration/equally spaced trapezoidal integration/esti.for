*******************************************************************
*esti=equally spaced trapezoidal integration（等间距梯形积分） 	
*b,a分别为积分上下限
*******************************************************************
	real function esti(a,b,n,f)
	integer n
	real a,b,delta
	real f(0:n)
	delta=(b-a)/n
	esti=(f(0)+f(n))/2.0
	do i=1,n-1,1
	  esti=esti+f(i)
	enddo
	esti=esti*delta
	end function
