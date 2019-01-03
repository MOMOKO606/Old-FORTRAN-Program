**************************************************
*  积分一步的变步长欧拉方法子程序测试程序：	
**************************************************	
	external f
	dimension y(3),d(3),yy(3),r(3),p(3)
	double precision t,y,h
	t=0.0
	data y/-1.0,0.0,1.0/
	m=3
	h=0.01
	eps=0.00001
	write(*,30) t,y(1),y(2),y(3)
	do i=1,10
	  call VSEM_1S(t,y,m,h,eps,f)
	  t=t+h
	  write(*,30) t,y(1),y(2),y(3)
	enddo
30	format(1x,'t=',f4.2,3x,'y(1)=',d13.6,3x,'y(2)=',d13.6,3x,
     &       'y(3)=',d13.6)
	end

	subroutine VSEM_1S(t,y,m,h,eps,f)
	dimension y(m),yy(m),p(m),d(m),r(m)
	double precision yy,y,h,t,tt,p,d,r,hh
	hh=h
	n=1
	sum=1.0+eps
	do i=1,m
	  yy(i)=y(i)
	  r(i)=y(i)
	enddo
10	do j=1,n
        tt=t+(j-1)*hh
	  call f(tt,y,m,d)
	  do i=1,m
	    p(i)=y(i)+hh*d(i)
	  enddo
	  tt=t+j*hh
	  call f(tt,p,m,d)
	  do i=1,m
	    d(i)=y(i)+hh*d(i)
	    y(i)=(p(i)+d(i))/2.0
	  enddo
	enddo  
	sum=0.0
	do i=1,m
	  x=abs(y(i)-r(i))
	  if(x.gt.sum) sum=x
	enddo
	if(sum.ge.eps)then
	  hh=hh/2.0
	  n=n+n
	  do i=1,m
	    r(i)=y(i)
	    y(i)=yy(i)
	  enddo
	  goto 10	
	endif
	end subroutine 

	subroutine f(t,y,m,d)
	dimension y(m),d(m)
	double precision t,y,d
	d(1)=y(2)
	d(2)=-y(1)
	d(3)=-y(3)
	end subroutine