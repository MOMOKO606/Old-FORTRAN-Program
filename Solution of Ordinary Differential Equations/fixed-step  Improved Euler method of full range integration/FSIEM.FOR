**************************************************************************************************
*  定步长改进欧拉方法子程序FSIEM（fixed-step Improved Euler method of full range integration）
*  m——方程组中方程的个数。
*  n——步数。
*  h——步长。
*  t——微分方程组的初值点。
*  y——双精度实型一维数组，在初值点t处，m个未知函数的初值。
*  f——子程序名，输入参数，用于计算方程组中各方程的右端函数值。
*  z——返回值，双精度实型二维数组，存放m个函数，n个点上的值。
**************************************************************************************************	
	subroutine FSIEM(m,n,h,t,y,f,z)
	dimension y(m),z(m,n),p(m),q(m),d(m)
	double precision y,z,t,h,p,d,tt,q
	do i=1,m
	  z(i,1)=y(i)
	enddo
	t=t-h
	do j=2,n
	  t=t+h
	  call f(t,y,m,d)
	  do i=1,m
	    p(i)=z(i,j-1)+h*d(i)
	  enddo
	  tt=t+h
	  call f(tt,p,m,d)
	  do i=1,m
	    q(i)=z(i,j-1)+h*d(i)
	    z(i,j)=(p(i)+q(i))/2.0
          y(i)=z(i,j)
	  enddo
	enddo   
	end subroutine
********************************************************************
*  另一种写法：
********************************************************************
	subroutine FSIEM(m,n,h,t,y,f,z)
	dimension y(m),z(m,n),d(m)
	double precision y,z,t,h,d
	do i=1,m
	  z(i,1)=y(i)
	enddo
	do j=2,n
	  t=t+(j-2)*h
	  call f(t,y,m,d)
	  do i=1,m
	    y(i)=z(i,j-1)+h*d(i)
	  enddo
	  t=t+h
	  call f(t,y,m,d)
	  do i=1,m
	    d(i)=z(i,j-1)+h*d(i)
	  enddo
	  do i=1,m
	    y(i)=(y(i)+d(i))/2.0
          z(i,j)=y(i)
	  enddo
	enddo   
	end subroutine