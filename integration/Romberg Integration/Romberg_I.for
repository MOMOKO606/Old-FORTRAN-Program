*********************************************************************
*  龙贝格求积法子程序Romberg_I（Romberg Integration）
*  b,a——分别为积分上下限。
*  eps——为精度要求。
*  f——双精度实型函数子程序名，输入参数，用于计算被积函数f(x)。
*  r——返回积分值。
*********************************************************************
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

***********************************************************************
*do while写法：
***********************************************************************
	subroutine Romberg_I(a,b,eps,f,r)
	dimension t(10)
	double precision a,b,p,t,t2,f,s,temp,x,s1,r
	n=1
	m=1
	h=b-a
	p=1.0+eps
	t(1)=h*0.5*(f(a)+f(b))
	do while((p.gt.eps).and.(m.le.10))
	  t2=0.0
	  do k=0,n-1
	    x=(k+0.5)*h
	    t2=t2+f(x)
	  enddo
	  t2=0.5*(t(1)+h*t2)
	  s=1.0
	  do k=1,m
	    s=s*4.0
	    r=(s*t2-t(k))/(s-1)
		temp=t(k)
		t(k)=t2
		t2=r
	  enddo
	  p=abs(r-temp)
	  n=n+n
	  h=h/2.0
	  m=m+1
	  t(m)=r   	    
	enddo
	end subroutine