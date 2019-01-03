*********************************************************************
*  变步长梯形积分子程序VSTI（Variable Step Trapezoidal Integration）
*  b,a——分别为积分上下限。
*  eps——为精度要求。
*  f——双精度实型函数子程序名，输入参数，用于计算被积函数f(x)。
*  t——返回积分值。
*********************************************************************
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

**************************************************************************
*  另一种写法
**************************************************************************
	subroutine VSTI(a,b,eps,f,t)
	double precision a,b,h,t1,t,s,x,f
	n=1
	h=b-a
	t1=(f(a)+f(b))*h/2.0
	p=eps+1.0
	do while(p.ge.eps)
	  s=0.0
	  do k=0,n-1,1
	    x=a+(k+0.5)*h
	    s=s+f(x)
	  enddo
	  t=(t1+s*h)/2.0
	  p=abs(t-t1)
	  t1=t
	  n=2*n
	  h=h/2.0
	enddo 
	end subroutine