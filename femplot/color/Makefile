
default:

colormap:
	cd oceanlib; ./prep.sh; cp oceanlib.dat ..
	cd matplotlib; cp matplotlib.dat ..
	cat matplotlib.dat oceanlib.dat > colormap.dat

clean:
	-rm -f oceanlib.dat matplotlib.dat

