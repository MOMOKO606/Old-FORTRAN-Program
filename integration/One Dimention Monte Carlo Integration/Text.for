*************************************************************
*蒙特卡洛法一维积分子程序测试程序
*************************************************************	
	external f,nrnd1
	real nrnd1
	double precision a,b,f,s
	a=2.5
	b=8.4
	call MCI_1D(a,b,f,s)
	write(*,*) s
	end
	
	subroutine MCI_1D(a,b,f,s)
	double precision  a,b,f,r,x,s,k
	real nrnd1
	r=1.0
	s=0.0
	do k=1,10000
	  x=a+(b-a)*nrnd1(r)
	  s=s+f(x)
	enddo
	s=(b-a)*s/10000.0
	end subroutine

	function f(x)
	double precision x
	f=x*x+sin(x)
	end function

	real function nrnd1(r)
	double precision s,u,v,r
	s=65536.0
	u=2053.0
	v=13849.0
	r=mod(r,s)
!   取余的另一种写法
!	m=r/s
!	r=r-m*s
	r=u*r+v
	r=mod(r,s)
!   取余的另一种写法
!	m=r/s
!	r=r-m*s
	nrnd1=r/s
	end function