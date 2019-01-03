**********************************************************************	
*函数功能:计算一维数据中的数据个数n。	
**********************************************************************	
	integer function get_n(filename,k)
	real temp(1:k)
	character*(*) filename
	get_n=0
	open(11,file=filename,status='old')
	do while(.not.eof(11))
	  read(11,*) (temp(i),i=1,k)
	  get_n=get_n+1 
	end do
	close(11)
	end function

***********************************************************************
*子函数形式。
***********************************************************************
	subroutine get_n(filename,k,n)
	real temp(1:k)
	character*(*) filename
	n=0
	open(11,file=filename,status='old')
	do while(.not.eof(11))
	  read(11,*) (temp(i),i=1,k)
	  n=n+1 
	end do
	close(11)
	end subroutine
