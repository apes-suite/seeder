FROM ghcr.io/apes-suite/apes-pyenv:2025.4
ARG SEEDER_VERSION
ENV SEEDER_VERSION=${SEEDER_VERSION}

COPY . /tmp/seeder
RUN (cd /tmp/seeder && bin/waf configure install --notests --prefix=$VIRTUAL_ENV) && rm -rf /tmp/seeder
