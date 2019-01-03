********************************************************************
* 变步长梯形积分子程序测试程序	
********************************************************************
	external f
	double precision f,a,b,t
	a=0.0
	b=1.0
	eps=0.000001
	call VSTI(a,b,eps,f,t)
	write(*,*) t
	end


	subroutine VSTI(a,b,eps,f,t)
	double precision a,b,h,t1,t,s,x,f
	n=1
	h=b-a
	t1=(f(a)+f(b))*h/2.0
10	s=0.0
	do k=0,n-1,1
	  x=a+(k+0.5)*h
	  s=s+f(x)
	end do
	t=(t1+s*h)/2.0
	if (abs(t-t1).ge.eps)then
	  t1=t
	  n=2*n
	  h=h/2.0
	  goto 10
	endif
	end subroutine

	double precision function f(x)
	double precision x
	f=exp(-x*x)
	end function
		

