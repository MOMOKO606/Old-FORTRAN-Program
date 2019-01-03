***************************************************************
*  蒙特卡洛多重积分测试程序
***************************************************************
	external f
	dimension x(3),a(3),b(3)
	double precision x,a,b,s
	data a,b/3*1.0,3*2.0/
	n=3
	call MCI_nD(a,b,n,x,f,s)
	write(*,*) s
	end

	
	subroutine MCI_nD(a,b,n,x,f,s)
	dimension a(n),b(n),x(n)
	double precision a,b,x,r,s,t,f
	real nrnd1
	m=10000
	t=1.0
	r=1.0
	do i=1,m  
	  do j=1,n
	    x(j)=a(j)+(b(j)-a(j))*nrnd1(r)
 	  enddo   
	  s=s+f(n,x)
	enddo
	do j=1,n
	  t=t*(b(j)-a(j))  
	enddo
	s=s*t/m
	end subroutine

	real function nrnd1(r)
	double precision s,u,v,r
	s=65536.0
	u=2053.0
	v=13849.0
	m=r/s
	r=r-m*s
	r=u*r+v
	m=r/s
	r=r-m*s
	nrnd1=r/s
	end function

	double precision function f(n,x)
	dimension x(n)
	double precision x
	f=0.0
	do i=1,n
	  f=f+x(i)*x(i)
	enddo
	end function