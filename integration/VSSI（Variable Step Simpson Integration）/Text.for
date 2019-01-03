********************************************************************
* 变步长辛普森积分子程序测试程序	
********************************************************************
	external f
	double precision f,a,b,s
	a=0.0
	b=1.0
	eps=0.000001
	call VSSI(a,b,eps,f,s)
	write(*,*) s
	end


	subroutine VSSI(a,b,eps,f,s)
	double precision a,b,f,s,temp,x,t1,t2,s1
*写法1：
*	s1=(b-a)*(f(a)+4*f((b+a)/2.0)+f(b))/6.0
*	n=2
*	h=(b-a)/n
*	t1=h*(f(a)+2*f((b+a)/2.0)+f(b))/4.0
*写法2：
	n=1
	h=b-a
	t1=h*(f(a)+f(b))/2.0
	s1=t1
10	temp=0.0
	do k=0,n-1,1
	  x=a+(k+0.5)*h
	  temp=temp+f(x)
	enddo
	t2=(t1+temp*h)/2.0
	s=(4*t2-t1)/3.0
	if(abs(s-s1).ge.eps)then
	  t1=t2
	  s1=s
	  n=n+n
	  h=h/2.0
	  goto 10
	endif
	end subroutine

	double precision function f(x)
	double precision x
	f=log(1+x)/(1+x*x)
	end function
		

