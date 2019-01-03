********************************************************************
*自适应梯形求积分子程序测试程序	
********************************************************************
	external f
	double precision f,a,b,t
	a=-1.0
	b=1.0
	eps=0.000001
	call ATI(a,b,eps,f,t)
	write(*,*) t
	end


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
10	if(k.ne.0)then
	  h=s(k,5)
	  x=s(k,1)+h/2.0
	  f0=f(x)
	  t1=h*(s(k,3)+f0)/4.0
	  t2=h*(f0+s(k,4))/4.0
	  if(abs(s(k,6)-t1-t2).lt.s(k,7).or.k.ge.29)then
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
	  goto 10
	endif	
	end subroutine  

	double precision function f(x)
	double precision x
	f=1.0/(1.0+25*x*x)
	end function
		

