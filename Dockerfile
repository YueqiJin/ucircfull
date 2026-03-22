FROM ubuntu:24.04 AS builder

ARG DEBIAN_FRONTEND=noninteractive
ARG APP_PREFIX=/opt/ucircfull

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        build-essential \
        cmake \
        pkg-config \
        ca-certificates \
        git \
        cargo \
        rustc \
        libboost-all-dev \
        libseqan3-dev \
        libseqan2-dev \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /src
COPY . .

RUN cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
    && cmake --build build -j"$(nproc)" \
    && cmake --install build --prefix "${APP_PREFIX}"

FROM ubuntu:24.04 AS runtime

ARG DEBIAN_FRONTEND=noninteractive
ARG APP_PREFIX=/opt/ucircfull

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        ca-certificates \
        build-essential \
        python3 \
        python3-pip \
        python3-venv \
        samtools \
        seqkit \
        libboost-all-dev \
    && python3 -m venv /opt/porechop-venv \
    && /opt/porechop-venv/bin/pip install --no-cache-dir --upgrade pip \
    && /opt/porechop-venv/bin/pip install --no-cache-dir https://github.com/rrwick/Porechop/archive/refs/tags/v0.2.4.tar.gz \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder ${APP_PREFIX} ${APP_PREFIX}

ENV PATH="/opt/porechop-venv/bin:${APP_PREFIX}/bin:${PATH}"
WORKDIR /data

ENTRYPOINT ["ucircfull"]
CMD ["--help"]