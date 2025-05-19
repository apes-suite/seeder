FROM ghcr.io/apes-suite/apes-pyenv:2025.4
ARG SEEDER_VERSION
ENV SEEDER_VERSION=${SEEDER_VERSION}

RUN bin/waf configure install --notests --prefix=$VIRTUAL_ENV
