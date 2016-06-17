
default:

prep: colormap
colormap:
	cd matplotlib; cp matplotlib.dat ..
	cd oceanlib; make prep; cp oceanlib.dat ..
	cd oceanlib_txt; make prep; cp oceanlib_txt.dat ..
	cat matplotlib.dat oceanlib.dat oceanlib_txt.dat > colormap.dat

clean:
	-rm -f matplotlib.dat oceanlib.dat oceanlib_txt.dat
	cd matplotlib; make clean
	cd oceanlib; make clean
	cd oceanlib_txt; make clean

