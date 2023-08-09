# docker build -t draw-base-image:0.1 .

FROM python:3.10-bullseye

WORKDIR /app

COPY . .

RUN pip install --no-cache-dir -r requirements.txt

RUN pip install .