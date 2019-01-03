**************************************************
*  全区间积分的定步长龙格-库塔法子程序测试程序：	
**************************************************
	external f
	dimension z(3,11),y(3),d(3)
	double precision t,y,h,f,z
	data y/-1.0,0.0,1.0/
	t=0.0
	h=0.01
	m=3
	n=11
	call FSRK(t,y,m,n,h,f,z)
	do j=1,n
	  t=(j-1)*h
	  write(*,10) t
	  write(*,20) (z(i,j),i=1,3)
	enddo
10    format(1x,'t=',f7.3)
20    format(1x,'y(1)=',d13.6,3x,'y(2)=',d13.6,3x,'y(3)=',d13.6)
	end
	
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

	subroutine f(t,y,m,d)
	dimension y(m),d(m)
	double precision t,y,d
	d(1)=y(2)
	d(2)=-y(1)
	d(3)=-y(3)
	end subroutine
