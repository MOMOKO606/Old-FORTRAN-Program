**************************************************************************************************
*  积分一步的变步长龙格-库塔法子程序VSRK_1S（Variable step Runge-Kutta method of one-step integration）
*  m——方程组中方程的个数。
*  h——步长。
*  t——微分方程组的初值点。
*  y——双精度实型一维数组，在初值点t处，m个未知函数的初值,同时也用来存储返回值。
*  f——子程序名，输入参数，用于计算方程组中各方程的右端函数值。
*  eps——精度要求。
**************************************************************************************************
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
**************************************************************************************************
*  do while写法：
**************************************************************************************************
	subroutine VSRK_1s(t,y,m,h,f,eps)
	dimension y(m),yy(m),a(4),r(m),c(m),b(m),d(m)
	double precision y,yy,a,r,c,b,d,h,hh,p,q,x,t,tt
	hh=h
	n=1
	p=1.0+eps
	do i=1,m
	  yy(i)=y(i)
	enddo
	do while(p.ge.eps)
	  a(1)=hh/2.0
	  a(2)=a(1)
	  a(3)=hh
	  a(4)=a(3)
	  tt=t-hh
	  do i=1,m
	    r(i)=y(i)	    
		y(i)=yy(i)
	  enddo
	  do j=1,n
	    tt=tt+hh
	    call f(tt,y,m,d)
	    do i=1,m
	      c(i)=y(i)
	      b(i)=y(i)
	    enddo
	    do k=1,3
	      do i=1,m
	        y(i)=c(i)+a(k)*d(i)
	        b(i)=b(i)+a(k+1)*d(i)/3.0
	      enddo
	      x=tt+a(k)
	      call f(x,y,m,d)
	    enddo
	    do i=1,m
	      y(i)=b(i)+hh*d(i)/6.0
	    enddo		
	  enddo  
	  p=0.0
	  do i=1,m
	    q=abs(y(i)-r(i))
	    if(q.gt.p) p=q
        enddo
	  n=n+n
	  hh=hh/2.0
	enddo
	end subroutine

