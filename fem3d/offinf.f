
	program offinf

c shows content of offline data file

	implicit none

	integer it,nkn,nel,nrec,iu
	character*60 name

	name = 'off.dat'

	nrec = 0
	iu = 1

	open(iu,file='off.dat',status='old',form='unformatted',err=99)

    1	continue

          read(iu,err=98,end=2) it,nkn,nel
	  nrec = nrec + 1
	  write(6,*) nrec,it,nkn,nel

	  read(iu)
	  read(iu)
	  read(iu)
	  read(iu)
	  read(iu)

	  goto 1
    2	continue

	close(iu)

	stop
   98	continue
	stop 'error stop offinf: reading file'
   99	continue
	stop 'error stop offinf: opening file'
	end

