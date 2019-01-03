********************************************************************
*一维积分连分式法子程序测试程序	
********************************************************************
	external f
	double precision f,a,b,t
	a=0.0
	b=4.3
	eps=0.000001
	call CFI_1D(a,b,eps,f,t)
	write(*,20) t
20	format(1x,'T=',d15.6)
	end


	subroutine CFI_1D(a,b,eps,f,t)
	dimension bb(0:10),s(0:10)
	double precision a,b,f,t,x,bb,s,u,temp,t1,t2,g
	n=1
	m=1
	h=b-a
	t1=h*(f(a)+f(b))/2.0
	bb(0)=t1
	s(0)=t1
10	t2=0.0
	hh=b-a
	do k=0,n-1
	  x=a+(k+0.5)*h
	  t2=t2+f(x)
	enddo
	t2=(t1+h*t2)/2.0
	u=t2
	do k=1,m
	  if(abs(u-bb(k-1).eq.0.0))then
	    u=sign(1.0,(h/2.0-hh))
	  else
	    u=(h/2.0-hh)/(u-bb(k-1))
	  endif
	    hh=hh/2.0
	enddo  
	h1=h
	bb(m)=u
	g=bb(m)
	do k=m,1,-1
	  if(g.eq.0.0)then
	    temp=bb(k-1)-h1
	  else
	    temp=bb(k-1)-h1/g
	  endif
	  h1=h1*2
	  g=temp
	enddo
	s(m)=temp
	if((abs(s(m)-s(m-1)).ge.eps).and.(m.lt.10))then
	  m=m+1
	  t1=t2
	  h=h/2.0
	  n=n+n
	  goto 10
	endif
	t=s(m)
	end subroutine

	double precision function f(x)
	double precision x
	f=exp(-x*x)
	end function
		

