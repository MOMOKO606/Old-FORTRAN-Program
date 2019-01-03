**********************************************************
*  切比雪夫积分子程序测试程序
**********************************************************	
	external f
	double precision a,b,f,s
	a=2.5
	b=8.4
	eps=0.000001
	call Cheb_I(a,b,eps,f,s)	
	write(*,*) s
	end


	subroutine Cheb_I(a,b,eps,f,s)	
	dimension t(5)
	double precision a,b,aa,bb,t,h,r,x,temp,f,s
	data t/-0.8324975,-0.3745414,0.0,0.3745414,0.8324975/
	n=1
10	h=(b-a)/n
	s=0.0
	do k=1,n
	  aa=a+(k-1)*h
	  bb=a+k*h
	  r=0.0
	  do i=1,5
		x=((bb-aa)*t(i)+(bb+aa))/2.0
		r=r+f(x)
	  enddo
	  r=r*h/5.0
	  s=s+r
	enddo
	if(abs(s-temp).ge.eps)then
	  temp=s
	  n=n+1
	  goto 10
	endif
	end subroutine

	function f(x)
	double precision x
	f=x*x+sin(x)
	end function
