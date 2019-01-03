**************************************************************************************************
*  积分一步的变步长基尔法子程序VSGM_1S（Variable step Gill method of one-step integration）
*  m——方程组中方程的个数。
*  h——步长。
*  t——微分方程组的初值点。
*  y——双精度实型一维数组，在初值点t处，m个未知函数的初值,同时也用来存储返回值。
*  f——子程序名，输入参数，用于计算方程组中各方程的右端函数值。
*  eps——精度要求。
**************************************************************************************************
	subroutine VSGM_1s(t,y,m,h,f,q,eps)
	dimension a(4),b(4),c(4),y(m),yy(m),r(m),q(m),d(m),qq(m)
	double precision a,b,c,y,yy,r,q,d,sum,h,hh,t,tt,x,temp,qq
	data b/0.5,0.29289322,1.70710678,0.16666667/
	data c/2.0,1.0,1.0,2.0/
	do i=1,m
	  yy(i)=y(i)
	  qq(i)=q(i)
	enddo
	sum=1.0+eps
	hh=h
	n=1
10	if(sum.ge.eps)then
	  do i=1,m
	    r(i)=y(i)
	    y(i)=yy(i)
	    q(i)=qq(i)
	  enddo
	  a(1)=0.0
	  a(2)=hh/2.0
	  a(3)=a(2)
	  a(4)=hh
	  tt=t-hh
	  do j=1,n
	    tt=tt+hh
	    do k=1,4
	      x=tt+a(k)
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
	  sum=0.0
	  do i=1,m
	    p=abs(y(i)-r(i))
	    if(p.gt.sum) sum=p
	  enddo
	  n=n+n
	  hh=hh/2.0
	  goto 10
	endif
	end subroutine

	  
