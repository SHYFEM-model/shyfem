#
# $Id: Makefile,v 1.00 2020-11-06 :gmica ASC CMCC Exp $
#
#--------------------	FILES	---------------------

DIRS=external/gotm src tools/partitioner tools/pre tools/2nc

#-------------------------------------------------------------------------------------------------

all:
	mkdir -p libs/mod bin
	for DIR in $(DIRS) ;do make -C $$DIR; done	  

clean:
	for DIR in $(DIRS) ;do make -C $$DIR clean; done	  

#----------------------------------------------------------------------------------------------------------------
