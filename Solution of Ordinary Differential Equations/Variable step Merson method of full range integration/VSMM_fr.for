**************************************************************************************************
*  全区间积分的变步长默森方法子程序VSMM_fr（Variable step Merson method of full range integration）
*  m——方程组中方程的个数。
*  n——步数。
*  h——步长。
*  t——微分方程组的初值点。
*  y——双精度实型一维数组，在初值点t处，m个未知函数的初值。
*  f——子程序名，输入参数，用于计算方程组中各方程的右端函数值。
*  eps——精度要求。
*  z——返回值，双精度实型二维数组，存放m个函数，n个点上的值。
**************************************************************************************************
	subroutine VSMM_fr(t,y,m,n,h,f,eps,z)
	dimension z(m,n),y(m),r(m),a(m),b(m),c(m),d(m)
	double precision z,y,r,a,b,c,d,t,t1,tt,x,hh,h,p,temp,s
	do i=1,m
	  z(i,1)=y(i)
	enddo
	t1=t
	do j=2,n
	  t1=t+(j-2)*h
	  nn=1
	  hh=h
	  p=1.0+eps
	  do while(p.ge.eps)
	    do i=1,m
	      r(i)=y(i)
	      y(i)=z(i,j-1)
	    enddo
	    tt=t1-hh
	    do l=1,nn
	      tt=tt+hh
	      call f(tt,y,m,d)
	      do i=1,m
	        a(i)=d(i)
	        y(i)=y(i)+hh*d(i)/3.0
	      enddo
	      x=tt+hh/3.0
	      call f(x,y,m,d)
	      do i=1,m
	        b(i)=d(i)
	        y(i)=y(i)+hh*(d(i)-a(i))/6.0
	      enddo
		  call f(x,y,m,d)
	      do i=1,m
!			b(i)=d(i)
		    y(i)=y(i)+3.0*hh*(d(i)-4.0*(b(i)+a(i)/4.0)/9.0)/8.0
110			b(i)=d(i)
		  enddo
		  x=tt+hh/2.0
		  call f(x,y,m,d)
		  do i=1,m
		    temp=d(i)-15.0*(b(i)-a(i)/5.0)/16.0  
	        y(i)=y(i)+2.0*hh*temp
		    c(i)=d(i)
		  enddo
	      x=tt+hh
	      call f(x,y,m,d)
	      do i=1,m
	        temp=d(i)-8.0*(c(i)-9.0*(b(i)-2.0*a(i)/9.0)/8.0)
	        y(i)=y(i)+hh*temp/6.0
	      enddo
	    enddo
	    p=0.0
	    do i=1,m
	      s=abs(y(i)-r(i))
	      if(s.gt.p) p=s
	    enddo
	    hh=hh/2.0
	    nn=nn+nn
	  enddo
	  do i=1,m
	    z(i,j)=y(i)
	  enddo
	enddo 
	end subroutine 

