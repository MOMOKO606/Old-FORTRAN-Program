*******************************************************************
*uesti=unequally spaced trapezoidal integration（不等间距梯形积分） 	
*******************************************************************
	real function uesti(n,x,f)
	real x(0:n),f(0:n) 
	uesti=0.0
	do i=1,n
	  uesti=uesti+(f(i)+f(i-1))*(x(i)-x(i-1))
	enddo
	  uesti=uesti/2.0
	end function