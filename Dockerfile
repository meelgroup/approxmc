FROM ubuntu:16.04 as builder

LABEL maintainer="Mate Soos"
LABEL version="1.0"
LABEL Description="Approxmc"

# get curl, etc
RUN apt-get update && apt-get install --no-install-recommends -y software-properties-common \
    && rm -rf /var/lib/apt/lists/*

RUN add-apt-repository -y ppa:ubuntu-toolchain-r/test \
    && rm -rf /var/lib/apt/lists/*

RUN apt-get update \
    && apt-get install --no-install-recommends -y libboost-program-options-dev gcc g++ make cmake zlib1g-dev wget autoconf automake make libtool \
    && rm -rf /var/lib/apt/lists/*

# get M4RI
WORKDIR /
RUN wget https://bitbucket.org/malb/m4ri/downloads/m4ri-20140914.tar.gz \
    && tar -xvf m4ri-20140914.tar.gz
WORKDIR m4ri-20140914
RUN ./configure \
    && make \
    && make install \
    && make clean

# build CMS
WORKDIR /
RUN wget https://github.com/msoos/cryptominisat/archive/5.6.6.tar.gz \
    && tar -xvf 5.6.6.tar.gz
WORKDIR /cryptominisat-5.6.6
RUN mkdir build
WORKDIR /cryptominisat-5.6.6/build
RUN cmake -DSTATICCOMPILE=ON -DUSE_GAUSS=ON .. \
    && make -j6 \
    && make install \
    && rm -rf *

# build approxmc
USER root
COPY . /home/solver/approxmc
WORKDIR /home/solver/approxmc
RUN mkdir build
WORKDIR /home/solver/approxmc/build
RUN cmake -DSTATICCOMPILE=ON .. \
    && make -j6 \
    && make install \
    && rm -rf *

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

