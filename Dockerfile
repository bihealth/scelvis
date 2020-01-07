# Dockerfile for building scelvis Docker image.
#
# The official Docker image is built automatically through Bioconda and is available on quay.io.  However, in
# some situations, it's useful to manually build a Docker image, e.g., when putting intermediate versions into
# private Docker registries.

FROM python:3.7.4

COPY . /app
WORKDIR /app

RUN pip install -e .

EXPOSE 3050
CMD ["scelvis", "run", "--fake-data", "--port", "3050"]
