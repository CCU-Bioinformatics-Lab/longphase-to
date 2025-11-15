# ========= Build stage =========
FROM ubuntu:20.04 AS build
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
      git g++ gcc autoconf automake libtool make \
      zlib1g-dev libbz2-dev liblzma-dev ca-certificates  && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /src
RUN git clone --depth 1 \
    https://github.com/CCU-Bioinformatics-Lab/longphase-to.git .

# make longphase-to binary
RUN autoreconf -i && \
    ./configure --prefix=/usr/local && \
    make -j"$(nproc)" && \
    mv ./longphase-to /usr/local/bin/longphase-to && \
    strip --strip-unneeded /usr/local/bin/longphase-to || true

# ========= Runtime stage =========
FROM ubuntu:20.04 AS runtime
ARG DEBIAN_FRONTEND=noninteractive
ARG VERSION=v1.0.0
# install runtime dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
      zlib1g libbz2-1.0 liblzma5 libstdc++6 libgomp1 && \
    rm -rf /var/lib/apt/lists/*

# copy installed files
COPY --from=build /usr/local /usr/local

# OCI labels
LABEL org.opencontainers.image.title="longphase-to"
LABEL org.opencontainers.image.source="https://github.com/CCU-Bioinformatics-Lab/longphase-to"
LABEL org.opencontainers.image.version="${VERSION}"
LABEL org.opencontainers.image.authors="CCU Bioinformatics Lab"

ENV PATH="/usr/local/bin:${PATH}"
WORKDIR /work

ENTRYPOINT ["longphase-to"]
CMD ["--help"]
