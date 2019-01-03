*****************************************************************************
*  高斯法多重积分子程序测试程序
*****************************************************************************	
	external fs,f
	dimension js(3),x(3)
	double precision x,f,s
	data js/4,4,4/
	n=3
	call MIGM(n,js,x,fs,f,s)
	write(*,*) s
	end
	
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
	      return
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

	subroutine fs(j,n,x,dn,up)
	dimension x(n)
	double precision x,t,dn,up
	if (j.eq.1)then
	  dn=0.0
	  up=1.0
	elseif(j.eq.2)then
	  dn=0.0
	  up=sqrt(1.0-x(1)*x(1))
	elseif(j.eq.3)then
	  t=x(1)*x(1)+x(2)*x(2)
	  dn=sqrt(t)
	  up=sqrt(2.0-t)
	endif
	end subroutine

	double precision function f(n,x)
	dimension x(n)
	double precision x
	f=x(3)*x(3)
	end function


	
