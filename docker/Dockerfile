FROM python:3.8-buster

MAINTAINER Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>
LABEL org.opencontainers.image.source https://github.com/bihealth/scelvis

ARG app_git_url=https://github.com/bihealth/scelvis.git
ARG app_git_tag
ARG app_git_depth=1

ENV DEBIAN_FRONTEND noninteractive
ENV CUSTOM_STATIC_DIR /usr/src/app/local-static

# Copy source code into Docker image.
RUN mkdir -p /usr/src
RUN git clone --depth $app_git_depth --branch $app_git_tag $app_git_url /usr/src/app

# Install Python dependencies.
RUN cd /usr/src/app && \
    pip install -e .

# Expose port and set default command.
EXPOSE 3050
CMD ["scelvis", "run", "--fake-data", "--port", "3050"]
