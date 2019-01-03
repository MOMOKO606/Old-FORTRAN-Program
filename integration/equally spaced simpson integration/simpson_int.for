*************************************************************
*函数：计算等间距的辛普森积分
*************************************************************
	real function simpson_int(a,b,n,f)
	real a,b,h,s1,s2
	real f(0:n)
	if(mod(n,2).ne.0)then
	  m=n-1
	else
	  m=n
	endif
	h=(b-a)/n
	s1=0.0
	do i=1,m/2-1,1
	  s1=s1+f(2*i)
	enddo
	s2=0.0
	do i=1,m/2,1
	  s2=s2+f(2*i-1)
	enddo
	simpson_int=(f(0)+f(n)+2.0*s1+4.0*s2)*(h/3.0)
	if(mod(n,2).ne.0)then
	  simpson_int=simpson_int+(f(n-1)+f(n))*(h/2.0)
	endif
	end function

