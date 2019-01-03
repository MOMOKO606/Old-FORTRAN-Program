*****************************************************************************************
*     时间抽选法（DIT）：Decimation In Time
*     I:pr,pi 输入数据的实部、虚部
*     P: 1、自然序转为倒位序
*        2、计算WN
*        3、第一级蝶形运算
*        4、蝶形运算
*     O:fr,fi 输出数据的实部、虚部
*	
*     n：数据的总点数
*     k: n=2**k	
*     l:l等于0时做Fourier正变换，l不等于零时做Fourier反变换
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

