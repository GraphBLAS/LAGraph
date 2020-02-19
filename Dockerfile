ARG SS_RELEASE=3.2.0draft28
ARG BASE_CONTAINER=graphblas/suitesparse-graphblas:${SS_RELEASE}
FROM ${BASE_CONTAINER}

RUN apt update && apt install -y \
    build-essential \
    git \
    make \
    cmake \
    curl \
    sudo \
    zlib1g-dev \
    m4 \
    && rm -rf /var/lib/apt/lists/*

COPY . /LAGraph
WORKDIR /LAGraph
RUN make clean
RUN make library
# RUN make tests
RUN make install
