*********************************************************************
*  蒙特卡洛多重积分（Monte Carlo Multiple Integration）
*  b(j),a(j)——分别为第j层积分的上下限。
*  n——为积分重数。
*  f——双精度实型函数子程序名，输入参数，用于计算被积函数f(n,x)。
*  s——返回积分值。
*********************************************************************	
      subroutine MCI_nD(a,b,n,x,f,s)
	dimension a(n),b(n),x(n)
	double precision a,b,x,r,s,f,q
	real nrnd1,m
	m=10000.0
	q=10000.0
	r=1.0
10	if (m+1.0.ne.1.0)then
	  do j=1,n
	    x(j)=a(j)+(b(j)-a(j))*nrnd1(r)
 	  enddo   
	  m=m-1.0
	  s=s+f(n,x)/q
	  goto 10
	endif
	do j=1,n
	  s=s*(b(j)-a(j))  
	enddo
	end subroutine


**********************************************************************
*无if写法：
**********************************************************************
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