**************************************************
*  定步长改进欧拉方法子程序测试程序：	
**************************************************
	external f
	dimension y(3),z(3,11),d(3)
	double precision y,z,t,h,d
	m=3
	n=11
	t=0.0
	h=0.01
	y(1)=-1.0
	y(2)=0.0
	y(3)=1.0
	call fsiem(m,n,h,t,y,f,z)
      do j=1,n
	  x=(j-1)*h
	  write(*,20) x
 	  write(*,30) (z(i,j),i=1,m)
	enddo
20    format(1x,'t=',f)
30    format(1x,'y(1)=',d13.6,3x,'y(2)=',d13.6,3x,'y(3)=',d13.6)
	end 
	
	subroutine fsiem(m,n,h,t,y,f,z)
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

	subroutine f(t,y,m,d)
	dimension y(m),d(m)
	double precision t,y,d
	d(1)=y(2)
	d(2)=-y(1)
	d(3)=-y(3)
	end subroutine