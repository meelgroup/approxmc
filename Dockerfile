FROM ubuntu:16.04 as builder

LABEL maintainer="Mate Soos"
LABEL version="1.0"
LABEL Description="Approxmc"

# get curl, etc
RUN apt-get update && apt-get install --no-install-recommends -y software-properties-common
# RUN add-apt-repository -y ppa:ubuntu-toolchain-r/test
# RUN apt-get update
RUN apt-get install --no-install-recommends -y libboost-program-options-dev gcc g++ make cmake zlib1g-dev wget make

# get M4RI
WORKDIR /
RUN wget https://bitbucket.org/malb/m4ri/downloads/m4ri-20140914.tar.gz
RUN tar -xvf m4ri-20140914.tar.gz
WORKDIR m4ri-20140914
RUN ./configure
RUN make \
    && make install

# build CMS
WORKDIR /
RUN wget https://github.com/msoos/cryptominisat/archive/5.6.7.tar.gz
RUN tar -xvf 5.6.7.tar.gz
WORKDIR /cryptominisat-5.6.7
RUN mkdir build
WORKDIR /cryptominisat-5.6.7/build
RUN cmake -DSTATICCOMPILE=ON -DUSE_GAUSS=ON ..
RUN make -j6 \
    && make install

# build approxmc
USER root
COPY . /home/solver/approxmc
WORKDIR /home/solver/approxmc
RUN mkdir build
WORKDIR /home/solver/approxmc/build
RUN cmake -DSTATICCOMPILE=ON ..
RUN make -j6 \
    && make install

# set up for running
FROM alpine:latest
COPY --from=builder /usr/local/bin/approxmc /usr/local/bin/
ENTRYPOINT ["/usr/local/bin/approxmc"]

# --------------------
# HOW TO USE
# --------------------
# on file through STDIN:
#    zcat mizh-md5-47-3.cnf.gz | docker run --rm -i -a stdin -a stdout msoos/approxmc

# on a file:
#    docker run --rm -v `pwd`/myfile.cnf.gz:/in msoos/approxmc in

# echo through STDIN:
#    echo "1 2 0" | docker run --rm -i -a stdin -a stdout msoos/approxmc

# hand-written CNF:
#    docker run --rm -ti -a stdin -a stdout msoos/approxmc

