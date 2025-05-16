ARG PYENV_VERSION=main
FROM ghcr.io/apes-suite/apes-pyenv:${PYENV_VERSION}

SHELL ["/bin/bash", "-c"]

WORKDIR /home/apes
COPY . ./seeder
RUN cd seeder && bin/waf configure install --prefix=$VIRTUAL_ENV
