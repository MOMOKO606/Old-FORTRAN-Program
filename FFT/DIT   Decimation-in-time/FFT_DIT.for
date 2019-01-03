*****************************************************************************************
*     ʱ���ѡ����DIT����Decimation In Time
*     I:pr,pi �������ݵ�ʵ�����鲿
*     P: 1����Ȼ��תΪ��λ��
*        2������WN
*        3����һ����������
*        4����������
*     O:fr,fi ������ݵ�ʵ�����鲿
*	
*     n�����ݵ��ܵ���
*     k: n=2**k	
*     l:l����0ʱ��Fourier���任��l��������ʱ��Fourier���任
******************************************************************************************	
	subroutine FFT_DIT(pr,pi,n,k,fr,fi,l)
	dimension pr(n),pi(n),fr(n),fi(n)
	double precision pr,pi,fr,fi,p,q,s,vr,vi,poddr,poddi
*     real(8) pr(n),pi(n),fr(n),fi(n)
*	real(8) p,q,s,vr,vi,poddr,poddi
	do it=0,n-1
	  m=it
	  is=0
	  do i=0,k-1
	    j=m/2
	    is=is*2+(m-j*2)	
*         is=is+(m-j*2)*2**(k-1-i)	
		m=j
	  enddo
	  fr(it+1)=pr(is+1)
	  fi(it+1)=pi(is+1)
	enddo
	pr(1)=1.0
	pi(1)=0.0
	pr(2)=cos(6.283185306/n)
	pi(2)=-sin(6.283185306/n)
	if(l.ne.0) pi(2)=-pi(2)
	do i=3,n
	  p=pr(i-1)*pr(2)
	  q=pi(i-1)*pi(2)
	  s=(pr(i-1)+pi(i-1))*(pr(2)+pi(2))
	  pr(i)=p-q
	  pi(i)=s-p-q
	enddo
	do it=0,n-2,2
	  vr=fr(it+1)
	  vi=fi(it+1)
	  fr(it+1)=vr+fr(it+2)
	  fi(it+1)=vi+fi(it+2)
	  fr(it+2)=vr-fr(it+2)
	  fi(it+2)=vi-fi(it+2)
	enddo
	m=n/2
	nv=2
	do lo=0,k-2,1
	  m=m/2
	  nv=2*nv
        do it=0,(m-1)*nv,nv
	    do j=0,(nv/2)-1
	      p=pr(m*j+1)*fr(it+j+1+nv/2)
		  q=pi(m*j+1)*fi(it+j+1+nv/2)
		  s=pr(m*j+1)+pi(m*j+1)
		  s=s*(fr(it+j+1+nv/2)+fi(it+j+1+nv/2))
		  poddr=p-q
		  poddi=s-p-q
		  fr(it+j+1+nv/2)=fr(it+j+1)-poddr
		  fi(it+j+1+nv/2)=fi(it+j+1)-poddi
		  fr(it+j+1)=fr(it+j+1)+poddr
	      fi(it+j+1)=fi(it+j+1)+poddi
		enddo
	  enddo
	enddo
	if(l.ne.0)then
	  do i=1,n
	    fr(i)=fr(i)/n
	    fi(i)=fi(i)/n
	  enddo
	endif
	end subroutine  

