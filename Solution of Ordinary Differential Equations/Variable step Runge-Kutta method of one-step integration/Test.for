***********************************************************************
*  积分一步的变步长龙格-库塔法子程序测试程序
***********************************************************************
	external f
	dimension y(2)
	double precision t,y,h
	data y/0.0,1.0/
	m=2
	t=0.0
	h=0.1
	eps=0.00001
	do j=1,10
	  call VSRK_1s(t,y,m,h,f,eps)
	  t=t+h
	  write(*,20) t
	  write(*,30) y(1),y(2)
	enddo
20    format(1x,'t=',f4.2)
30    format(3x,'y(1)=',d13.6,3x,'y(2)=',d13.6)
	end

	subroutine VSRK_1s(t,y,m,h,f,eps)
	dimension b(m),y(m),yy(m),s(m),a(4),d(m),r(m)
	double precision b,d,y,yy,r,s,a,hh,h,t,tt,x,p,q
	hh=h
	tt=t-hh
	n=1
	do i=1,m
	  b(i)=y(i)
	  s(i)=y(i)
	enddo
10	a(1)=hh/2.0
	a(2)=a(1)
	a(3)=hh
	a(4)=a(3)
	do j=1,n
	  do i=1,m
	    r(i)=y(i)
	  enddo
	  tt=tt+hh
	  call f(tt,y,m,d)
	  do k=1,3
	    do i=1,m
	      yy(i)=r(i)+a(k)*d(i)
	      y(i)=y(i)+a(k+1)*d(i)/3.0
          enddo
	    x=tt+a(k)
	    call f(x,yy,m,d)
	  enddo
	  do i=1,m
	    y(i)=y(i)+hh*d(i)/6.0
        enddo
      enddo
	p=0.0
      do i=1,m
	  q=abs(y(i)-b(i))
	  if(q.gt.p) p=q
	enddo
	if(p.ge.eps)then
	  n=n+n
	  hh=hh/2.0
	  do i=1,m
	    b(i)=y(i)
	    y(i)=s(i)
	  enddo
	  goto 10
	endif  	  
	end subroutine

	subroutine f(t,y,m,d)
	dimension y(m),d(m)
	double precision t,y,d
	d(1)=y(2)
	d(2)=-y(1)
	end subroutine
