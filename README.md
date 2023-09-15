
# Draw

A pipeline for the analysis of DREEM and SEISMIC data.

## Installation

... TODO ...

## Docker Build

Local Build:
docker build -t draw-base-image:0.7 .

Push to Google Cloud:
gcloud builds submit --config=cloudbuild.yaml
docker tag draw-base-image:0.7 gcr.io/draw-385021/draw-base-image:0.7
docker push gcr.io/draw-385021/draw-base-image:0.7

## Documentation

The documentation is available on [Github Pages](https://rouskinlab.github.io/draw).

## Contributors

Yves Martin, Casper L'espérance-Kerckhoff
