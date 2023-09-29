
# Draw

A pipeline for the analysis of DREEM and SEISMIC data.

## Installation

... TODO ...

## Docker Build

Local Build:
docker build -t draw-base-image:0.8 .

Push to Google Cloud:
make sure cloudbuild.yaml has correct name and docker tag
gcloud builds submit --config=cloudbuild.yaml
<!-- docker tag draw-base-image:0.8 gcr.io/draw-385021/draw-base-image:0.8
docker push gcr.io/draw-385021/draw-base-image:0.8 -->

## Documentation

The documentation is available on [Github Pages](https://rouskinlab.github.io/draw).

## Contributors

Yves Martin, Casper L'espérance-Kerckhoff
