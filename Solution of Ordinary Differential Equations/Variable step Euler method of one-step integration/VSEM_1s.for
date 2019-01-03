**************************************************************************************************
*  ����һ���ı䲽���Ľ�ŷ�������ӳ���VSEM_1S��Variable step Euler method of one-step integration��
*  m�����������з��̵ĸ�����
*  h����������
*  t����΢�ַ�����ĳ�ֵ�㡣
*  y����˫����ʵ��һά���飬�ڳ�ֵ��t����m��δ֪�����ĳ�ֵ,ͬʱҲ�����洢����ֵ��
*  f�����ӳ�������������������ڼ��㷽�����и����̵��Ҷ˺���ֵ��
*  eps��������Ҫ��
**************************************************************************************************
*  д��1��
**************************************************************************************************	
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
**************************************************************************************************
*  д��2(do whileд��)��
**************************************************************************************************
	subroutine VSEM_1S(t,y,m,h,eps,f)
	dimension y(m),yy(m),p(m),d(m),r(m)
	double precision yy,y,h,t,tt,p,d,r,hh
	hh=h
	n=1
	sum=1.0+eps
	do i=1,m
	  yy(i)=y(i)
	enddo
	do while(sum.ge.eps)
	  do i=1,m
	    r(i)=y(i)
	    y(i)=yy(i)
	  enddo
	  do j=1,n
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
	  hh=hh/2.0
	  n=n+n
	enddo
	end subroutine 
**************************************************************************************************
*  д��3(do whileд��)��
**************************************************************************************************
	subroutine VSEM_1S(t,y,m,h,eps,f,a,b,c,d)
	dimension y(m),a(m),c(m),d(m),b(m)
	double precision a,y,h,t,x,c,d,b,hh
	hh=h
	n=1
	p=1.0+eps
	do i=1,m
	  a(i)=y(i)
	enddo
10	if(p.ge.eps)then
	  do i=1,m
	    b(i)=y(i)
	    y(i)=a(i)
	  enddo
	  do j=1,n
	    do i=1,m
	      c(i)=y(i)
	    enddo
          x=t+(j-1)*hh
	    call f(x,y,m,d)
	    do i=1,m
	      y(i)=c(i)+hh*d(i)
	    enddo
	    x=t+j*hh
	    call f(x,y,m,d)
	    do i=1,m
	      d(i)=c(i)+hh*d(i)
	      y(i)=(y(i)+d(i))/2.0
	    enddo
	  enddo  
	  p=0.0
	  do i=1,m
	    q=abs(y(i)-b(i))
	    if(q.gt.p) p=q
	  enddo
	  hh=hh/2.0
	  n=n+n
	  goto 10
	endif
	end subroutine 