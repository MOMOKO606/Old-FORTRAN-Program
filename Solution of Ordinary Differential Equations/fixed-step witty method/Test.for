**************************************************
*  定步长纬梯方法子程序测试程序：	
**************************************************	
	external f
	dimension y(3),z(3,11),d(3)
	double precision y,z,d,t,h
	data y/-1.0,0.0,1.0/
	t=0.0
	h=0.1
	m=3
	n=11
	call FSWN(t,y,h,m,n,f,z)
	t=t-h
	do j=1,n
	    t=t+h
	    write(*,*) 't=',t
	    write(*,20) (z(i,j),i=1,m)
	enddo	
20    format(1x,'y(1)=',d13.6,3x,'y(2)=',d13.6,3x,'y(3)=',d13.6) 
	end

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

	subroutine f(t,y,m,d)
	dimension y(m),d(m)
	double precision t,y,d
	d(1)=y(2)
	d(2)=-y(1)
	d(3)=-y(3)
	end subroutine
