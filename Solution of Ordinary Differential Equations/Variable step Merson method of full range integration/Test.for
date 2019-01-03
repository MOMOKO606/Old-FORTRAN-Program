**************************************************************************
*  全区间积分的变步长默森方法子程序测试程序：
**************************************************************************
	external f
	dimension y(2),z(2,11)
	double precision t,y,z,h
	data y/0.0,1.0/
	m=2
	n=11
	h=0.1
	t=0.0
	eps=0.00001
	call VSMM_fr(t,y,m,n,h,f,eps,z)
	do j=1,n
	  t=(j-1)*h
	  write(*,10) t
	  write(*,20) (z(i,j),i=1,m)
	enddo
10    format(1x,'t=',f7.3)
20    format(1x,'y1=',d13.6,5x,'y2=',d13.6)
	end

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

	subroutine f(t,y,m,d)
	dimension y(m),d(m)
	double precision t,y,d
	d(1)=60.0*y(2)*(0.06+t*(t-0.6))
	d(2)=-60.0*y(1)*(0.06+t*(t-0.6))
	end subroutine