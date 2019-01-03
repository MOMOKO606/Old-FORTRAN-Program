********************************************************************************
*  一维积分连分式法 (One-dimension continued fraction Integration)
*  b,a——分别为积分上下限。
*  eps——为精度要求。
*  f——双精度实型函数子程序名，输入参数，用于计算被积函数f(x)。
*  t——返回积分值。
*  写法1：
********************************************************************************	
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
*******************************************************************************
*  写法2：
*******************************************************************************
	subroutine CFI_1D(a,b,eps,f,t,l)
	dimension bb(10),h(10)
	double precision a,b,f,t,bb,h,hh,t1,t2,x,u,g,s1
	m=1
	n=1
	hh=b-a
	h(1)=b-a
	t1=hh*(f(a)+f(b))/2.0
	s1=t1
	bb(1)=t1
10	t2=0.0
	do k=0,n-1
	  x=a+(k+0.5)*hh
	  t2=t2+f(x)
	enddo
	t2=(t1+t2*hh)/2.0
	m=m+1
	h(m)=h(m-1)/2.0
	u=t2
	do k=1,m-1
	  g=u-bb(k)
	  if(abs(g)+1.0.eq.1.0)then
	    g=sign(1.0d+35,g)
	    u=g*sign(1.0d0,h(m)-h(k))
	  else
	    u=(h(m)-h(k))/g
	  endif
	enddo
	bb(m)=u
	t=bb(m)
	do k=m,2,-1
	  if(abs(t)+1.0.eq.1.0)then
	    t=sign(1.0d+35,t)
	    t=t*sign(1.0d0,h(k-1))
	    t=bb(k-1)-t
	  else
	    t=bb(k-1)-h(k-1)/t
	  endif
	enddo
	if((abs(t-s1).ge.eps).and.(m.le.9))then
	  n=n+n
	  t1=t2
	  s1=t
	  hh=hh/2.0
	  goto 10
	endif
	if(m.ge.10)then
	  l=0
	else
	  l=1
	endif
	end subroutine
	    

	
