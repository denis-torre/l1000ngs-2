FROM ubuntu:17.10

RUN apt-get update && apt-get install -y python3 \
 python3-pip \
 python3-dev \
 python3-setuptools \
 nginx \
 uwsgi-core \
 libpcre3 \
 libpcre3-dev 

RUN mkdir website
COPY requirements.txt /website
WORKDIR /website

RUN pip3 install -r requirements.txt

COPY . .
RUN chmod +x boot.sh; chmod -R 777 /website/app/static;

ENTRYPOINT ./boot.sh