FROM python:3.6-slim

USER root

WORKDIR /app

ADD . /app
ARG dash_datadir=datasets/hnc
ADD ${dash_datadir} /app


RUN ls /app

RUN pip install --trusted-host pypi.python.org -r requirements.txt

EXPOSE 80

ENV NAME World
ENV DASH_DATADIR=${dash_datadir}

CMD ["python", "app.py"]
