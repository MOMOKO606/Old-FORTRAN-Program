*********************************************************************
*  自适应梯形积分子程序 （Adaptive Trapezoidal Integration） 
*  b,a——分别为积分上下限。
*  eps——为精度要求。
*  f——双精度实型函数子程序名，输入参数，用于计算被积函数f(x)。
*  t——返回积分值。
*  以下两个子程序前者“从左向右”分解，后者“从右向左”分解。
*********************************************************************
	subroutine ATI(a,b,eps,f,t)
	dimension s(30,7)
	double precision a,b,f,t,s,h,x,f0,t1,t2
	t=0.0
	k=1
	s(k,1)=a
	s(k,2)=b
	s(k,3)=f(a)
	s(k,4)=f(b)
	s(k,5)=b-a
	s(k,6)=(b-a)*(f(a)+f(b))/2.0
	s(k,7)=eps
*10	if(k.ne.0)then
	do while(k.ne.0)
	  h=s(k,5)
	  x=s(k,1)+h/2.0
	  f0=f(x)
	  t1=h*(s(k,3)+f0)/4.0
	  t2=h*(f0+s(k,4))/4.0
	  if((abs(s(k,6)-t1-t2).lt.s(k,7)).or.k.ge.29)then
	    t=t+t1+t2
		k=k-1
	  else
	    s(k+1,1)=s(k,1)
		s(k+1,2)=x
		s(k,1)=x
		s(k+1,3)=s(k,3)
		s(k+1,4)=f0
		s(k,3)=f0
		s(k,5)=s(k,5)/2.0
		s(k+1,5)=s(k,5)
		s(k,6)=t2
		s(k+1,6)=t1
		s(k,7)=s(k,7)/1.4
		s(k+1,7)=s(k,7)
		k=k+1
	  endif
	enddo
*	  goto 10
*	endif	
	end subroutine  


***************************************************************************
	subroutine ATI(a,b,eps,f,t)
	dimension s(30,7)
	double precision a,b,f,t,s,h,x,f0,f1,f3,p,t1,t2
	t=0.0
	f0=f(a)
	f1=f(b)
	p=h*(f0+f1)/2.0
	k=1
	s(k,1)=a
	s(k,2)=b
	s(k,3)=b-a
	s(k,4)=f0
	s(k,5)=f1
	s(k,6)=p
	s(k,7)=eps
10	if(k.ne.0)then
	  h=s(k,3)
	  x=s(k,1)+h/2.0
	  f3=f(x)
	  t1=h*(s(k,4)+f3)/4.0
	  t2=h*(f3+s(k,5))/4.0
	  if((abs(s(k,6)-t1-t2).lt.s(k,7)).or.(k.ge.29))then
	    t=t+t1+t2
		k=k-1
	  else
	    s(k+1,1)=x
		s(k+1,2)=s(k,2)
		s(k,2)=x
		s(k,3)=h/2.0
	    s(k+1,3)=s(k,3)
	    s(k+1,4)=f3
		s(k+1,5)=s(k,5)
		s(k,5)=f3
		s(k,7)=s(k,7)/1.4
		s(k+1,7)=s(k,7)
		s(k,6)=t1
		s(k+1,6)=t2
		k=k+1
	  endif
	  goto 10
	endif
	end subroutine
