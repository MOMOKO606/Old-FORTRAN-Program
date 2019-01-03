**************************************************************************
*  全区间积分的变步长基尔方法子程序测试程序：
**************************************************************************
	external f
	dimension y(3),z(3,11)
	double precision t,y,z,h
	data y/0.0,1.0,1.0/
	m=3
	n=11
	t=0.0
	h=0.1
	eps=0.000001
	call VSGM_fr(t,y,m,n,h,eps,f,z)
	do j=1,n
	  t=(j-1)*h
	  write(*,10) t
	  write(*,20) (z(i,j),i=1,3)
	enddo
10    format(1x,'t=',f7.3)
20    format(1x,'y1=',d13.6,5x,'y2=',d13.6,5x,'y3=',d13.6)
	end	 
	  
	subroutine VSGM_fr(t,y,m,n,h,eps,f,z)
	dimension z(m,n),a(4),b(4),c(4),qq(m),q(m),y(m),r(m),d(m)
	double precision z,a,b,c,qq,q,y,r,d,t,t1,tt,h,hh,x,temp,s,p
	data a/0.0,0.5,0.5,1.0/
	data b/0.5,0.29289321881,1.7071067812,0.166666667/
	data c/2.0,1.0,1.0,2.0/
	do i=1,m
	  z(i,1)=y(i)
	  qq(i)=0.0
	enddo
	do j=1,n-1
	  t1=t+(j-1)*h
	  n2=1
	  hh=h
	  s=eps+1.0
	  do while(s.ge.eps)
	    do i=1,m
	      r(i)=y(i)
	      y(i)=z(i,j)
	      q(i)=qq(i)
	    enddo
	    tt=t1-hh
	    do l=1,n2
	      tt=tt+hh
	      do k=1,4
	        x=tt+a(k)*hh
	        call f(x,y,m,d)
	        do i=1,m
	          d(i)=hh*d(i)
		      temp=b(k)*(d(i)-c(k)*q(i))
	          y(i)=y(i)+temp
	          if(k.le.3)then
	            q(i)=q(i)+3*temp-b(k)*d(i)
	          else
	            q(i)=q(i)+3*temp-0.5*d(i)
	          endif
	        enddo
	      enddo
	    enddo
	    s=0.0
	    do i=1,m
	      p=abs(y(i)-r(i))
	      if(p.gt.s) s=p
	    enddo
	    n2=n2+n2
	    hh=hh/2.0
	  enddo
	  do i=1,m
	    z(i,j+1)=y(i)
	    qq(i)=q(i)
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
