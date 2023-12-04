
# Draw

A pipeline for the analysis of DREEM and SEISMIC data.

## Installation

... TODO ...

## Docker Build

Local Build:
docker build -t draw-base-image:0.14 .

Push to Google Cloud:
make sure cloudbuild.yaml has correct name and docker tag
gcloud builds submit --config=cloudbuild.yaml

## Documentation

The documentation is available on [Github Pages](https://rouskinlab.github.io/draw).

## Contributors

Yves Martin, Casper L'espérance-Kerckhoff
