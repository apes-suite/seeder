ARG PYENV_VERSION=main
FROM ghcr.io/apes-suite/apes-pyenv:${PYENV_VERSION}

SHELL ["/bin/bash", "-c"]

COPY --chown=apes . ./seeder
RUN cd seeder && FC=mpif90 bin/waf configure install --mpicmd 'mpiexec --oversubscribe' --prefix=$VIRTUAL_ENV
