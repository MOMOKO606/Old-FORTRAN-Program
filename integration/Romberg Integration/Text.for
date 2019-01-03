********************************************************************
*龙贝格求积法子程序测试程序	
********************************************************************
	external f
	double precision f,a,b,r
	a=0.0
	b=1.0
	eps=0.000001
	call Romberg_I(a,b,eps,f,r)
	write(*,20) r
20	format(1x,'T=',d15.6)
	end


	subroutine Romberg_I(a,b,eps,f,r)
	dimension t(10)
	double precision a,b,f,r,t,t2,h,s,q,x,temp
	m=1
	n=1
	h=b-a
	t(1)=h*(f(a)+f(b))/2.0
10	t2=0.0
	do k=0,n-1
	  x=a+(k+0.5)*h
	  t2=t2+f(x)
	enddo
	t2=(t(1)+t2*h)/2.0
	s=1.0
	do k=1,m
	  s=s*4.0
	  q=(s*t2-t(k))/(s-1)
	  temp=t(k)
	  t(k)=t2
	  t2=q
	enddo
	if((abs(q-temp).ge.eps).and.(m.le.9))then
	  m=m+1
	  t(m)=q
	  h=h/2.0
	  n=n+n
	  goto 10
	endif
	  r=q
	end subroutine

	double precision function f(x)
	double precision x
	f=x/(4.0+x*x)
	end function
		

