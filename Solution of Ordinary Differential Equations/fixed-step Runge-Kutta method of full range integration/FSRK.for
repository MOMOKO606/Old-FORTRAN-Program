**************************************************************************************************
*  全区间积分的定步长龙格-库塔法子程序FSRK（fixed-step Runge-Kutta method of full range integration）
*  m——方程组中方程的个数。
*  n——步数。
*  h——步长。
*  t——微分方程组的初值点。
*  y——双精度实型一维数组，在初值点t处，m个未知函数的初值。
*  f——子程序名，输入参数，用于计算方程组中各方程的右端函数值。
*  z——返回值，双精度实型二维数组，存放m个函数，n个点上的值。
**************************************************************************************************	
	subroutine FSRK(t,y,m,n,h,f,z)
	dimension z(m,n),y(m),d(m),k1(m),k2(m),k3(m)
	double precision z,y,d,k1,k2,k3,t,h,x
	do i=1,m
	  z(i,1)=y(i)
	enddo
	t=t-h
	do j=2,n
	  t=t+h
	  call f(t,y,m,d)
	  x=t+h*0.5
	  do l=1,2
	    do i=1,m
	      y(i)=z(i,j-1)+0.5*h*d(i)
	      if(l.eq.1)then
		    k1(i)=d(i)
	      else
	        k2(i)=d(i)
	      endif
	    enddo
	    call f(x,y,m,d)
	  enddo
	  x=t+h
	  do i=1,m
	    k3(i)=d(i)
	    y(i)=z(i,j-1)+h*d(i)
	  enddo
	  call f(x,y,m,d)
	  do i=1,m
	    z(i,j)=z(i,j-1)+(k1(i)+2*(k2(i)+k3(i))+d(i))*h/6.0
	    y(i)=z(i,j)
	  enddo
	enddo
	end subroutine
*************************************************************************************
*  写法2：
*************************************************************************************
	subroutine FSRK(t,y,m,n,h,f,z)
	dimension z(m,n),y(m),d(m),b(m),a(4)
	double precision z,y,d,b,a,h,t,x,tt
	a(1)=h/2.0
	a(2)=a(1)
	a(3)=h
	a(4)=h
	x=t
	do i=1,m
	  z(i,1)=y(i)
	enddo
	do j=2,n
	  call f(t,y,m,d)
	  do i=1,3
	    b(i)=y(i)
	  enddo
	  do k=1,3
	    do i=1,m
		  y(i)=z(i,j-1)+a(k)*d(i)
	      b(i)=b(i)+a(k+1)*d(i)/3.0
	    enddo
	    tt=t+a(k)
	    call f(tt,y,m,d)
        enddo
	  do i=1,m
	    y(i)=b(i)+h*d(i)/6.0
	    z(i,j)=y(i)
	  enddo
	  t=t+h
	enddo
	t=x
	end subroutine
