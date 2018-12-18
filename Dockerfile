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

# set up build env
# RUN groupadd -r solver -g 433
# RUN useradd -u 431 -r -g solver -d /home/solver -s /sbin/nologin -c "Docker image user" solver
# RUN mkdir -p /home/solver/approxmc
# RUN chown -R solver:solver /home/solver

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
RUN git clone --depth 1 https://github.com/msoos/cryptominisat.git
WORKDIR /cryptominisat
RUN mkdir build
WORKDIR /cryptominisat/build
RUN cmake .. \
    && make -j2 \
    && make install \
    && rm -rf *

# build approxmc
USER root
COPY . /home/solver/approxmc
WORKDIR /home/solver/approxmc
RUN mkdir build
WORKDIR /home/solver/approxmc/build
RUN cmake .. \
    && make -j2 \
    && make install \
    && rm -rf *

# set up for running
FROM alpine:latest
COPY --from=builder /usr/local/bin/* /usr/local/bin/
COPY --from=builder /usr/local/bin/lib/* /usr/local/lib/
ENTRYPOINT ["/usr/local/bin/approxmc"]

# --------------------
# HOW TO USE
# --------------------
# on file through STDIN:
#    zcat mizh-md5-47-3.cnf.gz | docker run --rm -i -a stdin -a stdout approxmc

# on a file:
#    docker run --rm -v `pwd`/myfile.cnf.gz:/in approxmc in

# echo through STDIN:
#    echo "1 2 0" | docker run --rm -i -a stdin -a stdout approxmc

# hand-written CNF:
#    docker run --rm -ti -a stdin -a stdout approxmc

