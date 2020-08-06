FROM python:3.8.2-buster

WORKDIR /app/

ADD requirements.txt /app/requirements.txt

RUN pip3 install -r requirements.txt

ADD bot.py /app/bot.py
