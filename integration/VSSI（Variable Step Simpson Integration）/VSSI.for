********************************************************************
*  变步长辛普森积分子程序VSSI（Variable Step Simpson Integration）
*  b,a——分别为积分上下限。
*  eps——为精度要求。
*  f——双精度实型函数子程序名，输入参数，用于计算被积函数f(x)。
*  s——返回积分值。
********************************************************************
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
		
**********************************************************************
* 使用do while的写法：
**********************************************************************
	subroutine VSSI(a,b,eps,f,s)
	double precision a,b,h,p,t1,t2,f,x,s1,s
	n=1
	h=b-a
	p=1.0+eps
	s1=0.0
	t1=h*0.5*(f(a)+f(b))
	do while(p.gt.eps)
	  t2=0.0
	  do k=0,n-1
	    x=(k+0.5)*h
	    t2=t2+f(x)
	  enddo
	  t2=0.5*(t1+h*t2)
	  s=(4.0*t2-t1)/3.0
	  p=abs(s-s1)
	  t1=t2
	  s1=s
	  n=n+n
	  h=h/2.0
	enddo
	end subroutine
