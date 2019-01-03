************************************************************************
*  高斯法多重积分子程序MIGM（Multiple Integral Gauss Method）
*  n——n重积分。
*  js——第i重积分区间分为js(i)个子区间。
*  x——双精度实型一维数组，输入参数，n个积分变量。x(i),i=1,2,…,n。
*  fs——双精度实型函数子程序名，计算第i重积分上下限。
*  f——双精度实型函数子程序名，输入参数，用于计算被积函数值。
*  s——返回积分值。
************************************************************************
	subroutine MIGM(n,js,x,fs,f,s)
	dimension js(n),x(n)
	dimension t(5),c(5),d(2,11),is(2,11),cc(11)
	double precision x,f,s,d,cc,is,t,c,p,dn,up
	data t/-0.9061798459,-0.5384693101,0.0,0.5384693101,0.9061798459/
	data c/0.2369268851,0.4786286705,0.5688888889,0.4786286705,
     &       0.2369268851/
	d(1,n+1)=1.0
	d(2,n+1)=1.0
	m=1
10	do j=m,n
	  call fs(j,n,x,dn,up)
	  d(1,j)=0.5*(up-dn)/js(j)
	  cc(j)=d(1,j)+dn
	  x(j)=d(1,j)*t(1)+cc(j)
	  d(2,j)=0.0
	  is(1,j)=1.0
	  is(2,j)=1.0
	enddo
	j=n
30	k=is(1,j)
	if(j.eq.n)then
	  p=f(n,x)
	else
	  p=1.0
	endif
	d(2,j)=d(2,j)+d(1,j+1)*d(2,j+1)*p*c(k)
	is(1,j)=is(1,j)+1
	if(is(1,j).gt.5)then
	  if(is(2,j).ge.js(j))then
	    j=j-1
		if(j.eq.0)then
	      s=d(2,1)*d(1,1)
		endif
	    goto 30
	  endif
	  is(2,j)=is(2,j)+1
	  cc(j)=cc(j)+2.0*d(1,j)
	  is(1,j)=1.0
	endif
	k=is(1,j)
	x(j)=d(1,j)*t(k)+cc(j)
	if(j.eq.n)goto 30
	m=j+1
	goto 10
	end subroutine
**********************************************************************
* do while写法：
**********************************************************************
	subroutine MIGM(n,js,x,fs,f,s)
	dimension x(n),t(6),c(6),d(2,11),is(2,11),cc(11)
	double precision x,t,c,d,is,dn,up,cc,p,f,s
	data t/-0.9061798459,-0.5384693101,0.0,0.5384693101,0.9061798459,
     &       0.0/
	data c/0.2369268851,0.4786286705,0.5688888889,0.4786286705,
     &       0.2369268851,0.0/
	d(1,n+1)=1.0
	d(2,n+1)=1.0
	do j=1,n
	  call fs(j,n,x,dn,up)
	  d(1,j)=0.5*(up-dn)/js
	  cc(j)=d(1,j)+dn
	  x(j)=d(1,j)*t(1)+cc(j)
	  is(1,j)=1.0
	  is(2,j)=1.0
	  d(2,j)=0.0
	enddo
	j=n
	do while(j.gt.0)
	  do while(is(2,j).le.js)
	    do while(is(1,j).le.5)
	      if((q.eq.1.0).and.(j.ne.n))then
		    q=0.0 
	      else
	        if(j.eq.n)then
	          p=f(n,x)
	        else
	          p=1.0
	        endif
	        k=is(1,j)
 		    d(2,j)=d(2,j)+d(1,j+1)*d(2,j+1)*p*c(k) 
		    is(1,j)=is(1,j)+1
		    k=is(1,j)
		    x(j)=d(1,j)*t(k)+cc(j)
	      endif
		  if(j.ne.n)then
	        m=j+1
	        do j=m,n
	          call fs(j,n,x,dn,up)
	          d(1,j)=0.5*(up-dn)/js
	          cc(j)=d(1,j)+dn
	          x(j)=d(1,j)*t(1)+cc(j)
	          is(1,j)=1
	          is(2,j)=1
	  	      d(2,j)=0.0
		    enddo
	        j=n
	      endif
	    enddo
	    is(2,j)=is(2,j)+1
	    cc(j)=cc(j)+2.0*d(1,j)
	    x(j)=d(1,j)*t(1)+cc(j) 
	    is(1,j)=1
	    q=1.0
	  enddo
	  j=j-1
	  q=0.0
	enddo
	s=d(2,1)*d(1,1)
	end subroutine	   



