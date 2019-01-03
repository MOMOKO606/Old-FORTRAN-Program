**************************************************************************************************
*  全区间积分的变步长基尔方法子程序VSGM_fr（Variable step Gill method of full range integration）
*  m——方程组中方程的个数。
*  n——步数。
*  h——步长。
*  t——微分方程组的初值点。
*  y——双精度实型一维数组，在初值点t处，m个未知函数的初值。
*  f——子程序名，输入参数，用于计算方程组中各方程的右端函数值。
*  eps——精度要求。
*  z——返回值，双精度实型二维数组，存放m个函数，n个点上的值。
**************************************************************************************************	
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
	      
	
