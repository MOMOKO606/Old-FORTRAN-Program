***************************************************************
*  二重积分连分式法子程序测试程序
****************************************************************
	external fs,f
	double precision a,b,f,s
	a=0.0
	b=1.0
	eps=0.00005
	call CFI_2D(a,b,fs,f,eps,s)
	write(*,*) s
	end
	
	subroutine CFI_2D(a,b,fs,f,eps,s)
	external f,fs
	dimension h(10),c(10)
	double precision a,b,x,f,g,s,h1,h,s0,s1,s2,c,u,ga,gb
	n=1
	m=1
	call getgx(a,fs,f,eps,ga)
	call getgx(b,fs,f,eps,gb)
	h(1)=b-a
	h1=(b-a)*0.5
	s1=h1*(ga+gb)/2.0
	s0=s1
	c(1)=s1
10	s2=0.0
	do i=1,n
	  x=a+(2*i-1)*h1
	  call getgx(x,fs,f,eps,g)
	  s2=s2+g
	enddo
	s2=0.5*s1+s2*h1
	m=m+1
	h(m)=h(m-1)/2.0
	u=s2
	do k=1,m-1
	  if(abs(u-c(k))+1.0.eq.1.0)then
	    u=sign(1.0,(h(m)-h(k)))
	  else
	    u=(h(m)-h(k))/(u-c(k))
	  endif
	enddo
	c(m)=u
	s=c(m)
	do j=m,2,-1
	  if(abs(s)+1.0.eq.1.0)then
		s=c(j-1)-sign(1.0,(h(j-1)))
	  else
	    s=c(j-1)-h(j-1)/s
	  endif
	enddo
	if((abs(s-s0).ge.eps).and.(m.le.9))then
	  n=n+n
	  s1=s2
	  s0=s
	  h1=h1/2.0
	  goto 10
	endif
	end subroutine



	subroutine getgx(x,fs,f,eps,g)
	dimension hh(10),bb(10)
	double precision x,f,g,y1,y2,h2,hh,g0,g1,g2,bb,y,u
	n=1
	m=1
	call fs(x,y1,y2)
	hh(1)=y2-y1
	h2=(y2-y1)*0.5
	g1=h2*(f(x,y1)+f(x,y2))/2.0
	g0=g1
	bb(1)=g1
10	g2=0.0
	do i=1,n
	  y=y1+(2*i-1)*h2
	  g2=g2+f(x,y)
	enddo
	g2=0.5*g1+g2*h2
	m=m+1
	hh(m)=hh(m-1)/2.0
	u=g2
	do k=1,m-1
	  if(abs(u-bb(k))+1.0.eq.1.0)then
	    u=sign(1.0,(hh(m)-hh(k)))
	  else
	    u=(hh(m)-hh(k))/(u-bb(k))
	  endif
	enddo
	bb(m)=u
	g=bb(m)
	do j=m,2,-1
	  if(abs(g)+1.0.eq.1.0)then
		g=bb(j-1)-sign(1.0,(hh(j-1)))
	  else
	    g=bb(j-1)-hh(j-1)/g
	  endif
	enddo
	if((abs(g-g0).ge.eps).and.(m.le.9))then
	  n=n+n
	  g1=g2
	  g0=g
	  h2=h2/2.0
	  goto 10
	endif
	end subroutine

	subroutine fs(x,y1,y2)
	double precision x,y1,y2,q
	q=1-x*x
	y2=sqrt(q)
	y1=-sqrt(q)
	end subroutine

	double precision function f(x,y)
	double precision x,y,q
	q=x*x+y*y
	f=exp(q)
	end function

