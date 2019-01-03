********************************************************************
*  定步长纬梯方法子程序FSWN（fixed-step witty method）
*  m——方程组中方程的个数。
*  n——步数。
*  h——步长。
*  t——微分方程组的初值点。
*  y——双精度实型一维数组，在初值点t处，m个未知函数的初值。
*  f——子程序名，输入参数，用于计算方程组中各方程的右端函数值。
*  z——返回值，双精度实型二维数组，存放m个函数，n个点上的值。
********************************************************************
	subroutine FSWN(t,y,h,m,n,f,z)
	dimension y(m),z(m,n),d(m),a(m)
	double precision y,z,d,a,t,h,x
	do i=1,m
	  z(i,1)=y(i)
	enddo
	call f(t,y,m,d)
	do j=1,n-1
	  do i=1,m
	    y(i)=z(i,j)+0.5*h*d(i)
	    a(i)=d(i)
	  enddo
	  x=t+(j-0.5)*h
	  call f(x,y,m,d)
	  do i=1,m
	    z(i,j+1)=z(i,j)+h*d(i)
	    d(i)=2*d(i)-a(i)
	  enddo
	enddo
	end subroutine

