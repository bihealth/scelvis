FROM python:3.6-slim

USER root

WORKDIR /app

ADD . /app

RUN ls /app

RUN pip install --trusted-host pypi.python.org -r requirements.txt

EXPOSE 80

ENV NAME World

#CMD ["python", "app.py"]
CMD ["python", "hnc.py"]
