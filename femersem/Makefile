##-----------------------------------------------------------------------
## Makefile for making a new distribution release
##-----------------------------------------------------------------------


RELEASE=bfm-2.5e-system-2.1

fem:
	cd src; make fem
cleanall: clean
clean:
	cd src; make clean

all: distclean sources runenv

distclean:
	make -C src/ $@
	@$(RM)  ./bin/* 

sources:
	@sed -e "s/_BFMRELEASE_/$(RELEASE)/" ../bfm_env_proto.sh > ./bfm_env.sh
	@sed -e "s/_BFMRELEASE_/$(RELEASE)/" ../bfm_env_proto.csh > ./bfm_env.csh
	@chmod +x ./bfm_env.csh ./bfm_env.sh
	@mv ../bfm-dev ../$(RELEASE)
	@tar -zcvf ../$(RELEASE).tgz --exclude='*~' ../$(RELEASE)
	@mv ../$(RELEASE) ../bfm-dev

runenv:
	@tar -zcvf ../run-$(RELEASE).tgz --exclude='*.nc' ../bfm-run

list:
depend:

##-----------------------------------------------------------------------
##
##-----------------------------------------------------------------------
