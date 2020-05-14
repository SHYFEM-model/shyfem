
#------------------------------------------------------------------------
#
#    Copyright (C) 1985-2020  Georg Umgiesser
#
#    This file is part of SHYFEM.
#
#------------------------------------------------------------------------

FROM debian:stretch

RUN apt-get update && apt-get install -y \
	g++ \
	gfortran \
	libxt-dev \
	make

WORKDIR /root/shyfem-src

COPY . .

# install (soft-links current directory to /root/shyfem)
RUN make install

# build model
RUN make clean && make fem

ENV PATH "/root/shyfem/bin:${PATH}"

# run tests
CMD ./examples/mar_menor/run.sh

