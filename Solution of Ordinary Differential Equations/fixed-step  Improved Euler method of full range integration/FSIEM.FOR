**************************************************************************************************
*  �������Ľ�ŷ�������ӳ���FSIEM��fixed-step Improved Euler method of full range integration��
*  m�����������з��̵ĸ�����
*  n����������
*  h����������
*  t����΢�ַ�����ĳ�ֵ�㡣
*  y����˫����ʵ��һά���飬�ڳ�ֵ��t����m��δ֪�����ĳ�ֵ��
*  f�����ӳ�������������������ڼ��㷽�����и����̵��Ҷ˺���ֵ��
*  z��������ֵ��˫����ʵ�Ͷ�ά���飬���m��������n�����ϵ�ֵ��
**************************************************************************************************	
	subroutine FSIEM(m,n,h,t,y,f,z)
	dimension y(m),z(m,n),p(m),q(m),d(m)
	double precision y,z,t,h,p,d,tt,q
	do i=1,m
	  z(i,1)=y(i)
	enddo
	t=t-h
	do j=2,n
	  t=t+h
	  call f(t,y,m,d)
	  do i=1,m
	    p(i)=z(i,j-1)+h*d(i)
	  enddo
	  tt=t+h
	  call f(tt,p,m,d)
	  do i=1,m
	    q(i)=z(i,j-1)+h*d(i)
	    z(i,j)=(p(i)+q(i))/2.0
          y(i)=z(i,j)
	  enddo
	enddo   
	end subroutine
********************************************************************
*  ��һ��д����
********************************************************************
	subroutine FSIEM(m,n,h,t,y,f,z)
	dimension y(m),z(m,n),d(m)
	double precision y,z,t,h,d
	do i=1,m
	  z(i,1)=y(i)
	enddo
	do j=2,n
	  t=t+(j-2)*h
	  call f(t,y,m,d)
	  do i=1,m
	    y(i)=z(i,j-1)+h*d(i)
	  enddo
	  t=t+h
	  call f(t,y,m,d)
	  do i=1,m
	    d(i)=z(i,j-1)+h*d(i)
	  enddo
	  do i=1,m
	    y(i)=(y(i)+d(i))/2.0
          z(i,j)=y(i)
	  enddo
	enddo   
	end subroutine