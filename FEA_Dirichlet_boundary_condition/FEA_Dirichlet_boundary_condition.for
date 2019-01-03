	program FEA_Dirichlet_boundary_condition
	character *40 cmdfile,xyfile,i3file,nb1_u1file,ufile
	integer nd,ne,nd1,iw
	real,allocatable::xy(:,:),i3(:,:),sk(:,:),nb1(:),u1(:),u(:)
	call read_cmd(cmdfile,nd,ne,nd1,xyfile,i3file,nb1_u1file,ufile)
	allocate(xy(1:2,1:nd),i3(1:3,1:ne),nb1(1:nd1),u1(1:nd1),u(1:nd))
	call input_data(nd,ne,nd1,xyfile,i3file,nb1_u1file,
     &xy,i3,nb1,u1)
	call mbw(ne,i3,iw)
	allocate(sk(1:nd,1:iw))
	call uk1(nd,ne,iw,i3,xy,sk)
	call ub1(nd1,nb1,u1,nd,iw,sk,u)
	call ldlt(sk,nd,iw,u,ie)
	call output_u(nd,u,ufile)
	deallocate(xy,i3,sk,nb1,u1,u)
	end program

	subroutine read_cmd(cmdfile,nd,ne,nd1,xyfile,i3file,nb1_u1file,
     &ufile)
	character *(*) cmdfile,xyfile,i3file,nb1_u1file,ufile
	integer nd,ne,nd1
	write(*,*) 'please input cmdfile'
	read(*,*) cmdfile
	open(11,file=cmdfile,status='old')
	read(11,*) ne,nd,nd1
	read(11,*) xyfile
	read(11,*) i3file
	read(11,*) nb1_u1file
	read(11,*) ufile
	close(11)
	end subroutine

	subroutine input_data(nd,ne,nd1,xyfile,i3file,nb1_u1file,
     &xy,i3,nb1,u1)
	character*(*) xyfile,i3file,nb1_u1file
	integer nd,ne,nd1
	real xy(2,nd),i3(3,ne),nb1(nd1),u1(nd1)
	open(12,file=xyfile,status='old')
	do j=1,nd
	  read(12,*) (xy(i,j),i=1,2)
	enddo
	close(12)
	open(13,file=i3file,status='old')
	do j=1,ne
	  read(13,*) (i3(i,j),i=1,3)
	enddo
	close(13)
	open(14,file=nb1_u1file,status='old')
	do i=1,nd1
	  read(14,*) nb1(i),u1(i)
	enddo
	close(14)
	end subroutine

	subroutine mbw(ne,i3,iw)
	real i3(3,ne)
	iw=0
	do i=1,ne
	  m=max(abs(i3(1,i)-i3(2,i)),abs(i3(2,i)-i3(3,i)),
     &abs(i3(3,i)-i3(1,i)))
	if(m+1.gt.iw) iw=m+1
	enddo
	end subroutine

	subroutine uke1(x,y,ke)
	dimension x(3),y(3),a(3),b(3)
	real ke(3,3)
	a(1)=y(2)-y(3)
	a(2)=y(3)-y(1)
	a(3)=y(1)-y(2)
	b(1)=x(3)-x(2)
	b(2)=x(1)-x(3)
	b(3)=x(2)-x(1)
	s=2.0*(a(1)*b(2)-a(2)*b(1))
	do i=1,3
	  do j=1,i
	    ke(i,j)=(a(i)*a(j)+b(i)*b(j))/s
	  enddo
	enddo
	end subroutine

	subroutine uk1(nd,ne,iw,i3,xy,sk)
	real i3(3,ne),xy(2,nd),sk(nd,iw)
	dimension x(3),y(3)
	real ke(3,3)
	do i=1,nd
	  do j=1,iw
	    sk(i,j)=0
	  enddo
	enddo
	do l=1,ne
	  do j=1,3
	    i=i3(j,l)
		x(j)=xy(1,i)
		y(j)=xy(2,i)
	  enddo
	  call uke1(x,y,ke)
	  do j=1,3
	    nj=i3(j,l)
		do k=1,j
		  nk=i3(k,l)
		  if(nj.lt.nk) then
			nj=nj-nk+iw
			sk(nk,nj)=sk(nk,nj)+ke(j,k)
			nj=nj+nk-iw
		  else
			nk=nk-nj+iw
			sk(nj,nk)=sk(nj,nk)+ke(j,k)
		  endif
		enddo
	  enddo
	enddo
	end subroutine

	subroutine ub1(nd1,nb1,u1,nd,iw,sk,u)
	real nb1(nd1),u1(nd1),sk(nd,iw),u(nd)
	do i=1,nd
	  u(i)=0
	enddo
	do i=1,nd1
	  j=nb1(i)
	  sk(j,iw)=sk(j,iw)*1.E10
	  u(j)=sk(j,iw)*u1(i)
	enddo
	end subroutine

	SUBROUTINE LDLT(A,N,IW,P,IE)
	DIMENSION A(N,IW),P(N)
	DO 15 I=1,N
	IF(I.LE.IW) GOTO 20
	IT=I-IW+1
	GOTO 30
20	IT=1
30	K=I-1
	IF(I.EQ.1) GOTO 40
	DO 25 L=IT,K
	IL=L+IW-I
	B=A(I,IL)
	A(I,IL)=B/A(L,IW)
	P(I)=P(I)-A(I,IL)*P(L)
	MI=L+1
	DO 25 J=MI,I
	IJ=J+IW-I
	JL=L+IW-J
25	A(I,IJ)=A(I,IJ)-A(J,JL)*B
40	IF(A(I,IW).EQ.0.) GOTO 100
15	CONTINUE
	DO 45 J=1,N
	IF(J.LE.IW) GOTO 60
	IT=N-J+IW
	GOTO 70
60	IT=N
70	I=N-J+1
	P(I)=P(I)/A(I,IW)
	IF(J.EQ.1) GOTO 45
	K=I+1
	DO 65 MJ=K,IT
	IJ=I-MJ+IW
65	P(I)=P(I)-P(MJ)*A(MJ,IJ)
45	CONTINUE
	IE=0
	GOTO 110
100	IE=1
110	RETURN
	END SUBROUTINE

	subroutine output_u(nd,u,ufile)
	character*(*) ufile
	dimension u(nd)
	open (15,file=ufile,status='unknown')
	do i=1,nd
	  write(15,*) u(i)
	enddo
	close(15)
	end subroutine


