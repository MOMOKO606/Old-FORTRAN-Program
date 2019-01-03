**************************************************************************************************
*  ȫ������ֵı䲽��Ĭɭ�����ӳ���VSMM_fr��Variable step Merson method of full range integration��
*  m�����������з��̵ĸ�����
*  n����������
*  h����������
*  t����΢�ַ�����ĳ�ֵ�㡣
*  y����˫����ʵ��һά���飬�ڳ�ֵ��t����m��δ֪�����ĳ�ֵ��
*  f�����ӳ�������������������ڼ��㷽�����и����̵��Ҷ˺���ֵ��
*  eps��������Ҫ��
*  z��������ֵ��˫����ʵ�Ͷ�ά���飬���m��������n�����ϵ�ֵ��
**************************************************************************************************
	subroutine VSMM_fr(t,y,m,n,h,f,eps,z)
	dimension z(m,n),y(m),r(m),a(m),b(m),c(m),d(m)
	double precision z,y,r,a,b,c,d,t,t1,tt,x,hh,h,p,temp,s
	do i=1,m
	  z(i,1)=y(i)
	enddo
	t1=t
	do j=2,n
	  t1=t+(j-2)*h
	  nn=1
	  hh=h
	  p=1.0+eps
	  do while(p.ge.eps)
	    do i=1,m
	      r(i)=y(i)
	      y(i)=z(i,j-1)
	    enddo
	    tt=t1-hh
	    do l=1,nn
	      tt=tt+hh
	      call f(tt,y,m,d)
	      do i=1,m
	        a(i)=d(i)
	        y(i)=y(i)+hh*d(i)/3.0
	      enddo
	      x=tt+hh/3.0
	      call f(x,y,m,d)
	      do i=1,m
	        b(i)=d(i)
	        y(i)=y(i)+hh*(d(i)-a(i))/6.0
	      enddo
		  call f(x,y,m,d)
	      do i=1,m
!			b(i)=d(i)
		    y(i)=y(i)+3.0*hh*(d(i)-4.0*(b(i)+a(i)/4.0)/9.0)/8.0
110			b(i)=d(i)
		  enddo
		  x=tt+hh/2.0
		  call f(x,y,m,d)
		  do i=1,m
		    temp=d(i)-15.0*(b(i)-a(i)/5.0)/16.0  
	        y(i)=y(i)+2.0*hh*temp
		    c(i)=d(i)
		  enddo
	      x=tt+hh
	      call f(x,y,m,d)
	      do i=1,m
	        temp=d(i)-8.0*(c(i)-9.0*(b(i)-2.0*a(i)/9.0)/8.0)
	        y(i)=y(i)+hh*temp/6.0
	      enddo
	    enddo
	    p=0.0
	    do i=1,m
	      s=abs(y(i)-r(i))
	      if(s.gt.p) p=s
	    enddo
	    hh=hh/2.0
	    nn=nn+nn
	  enddo
	  do i=1,m
	    z(i,j)=y(i)
	  enddo
	enddo 
	end subroutine 

