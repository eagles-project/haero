#!/usr/bin/env bash

docker login
./build-ext-docker-image.sh Debug double 1
./build-ext-docker-image.sh Debug double 4
./build-ext-docker-image.sh Release double 1
./build-ext-docker-image.sh Release double 4
./build-ext-docker-image.sh Debug single 1
./build-ext-docker-image.sh Debug single 4
./build-ext-docker-image.sh Release single 1
./build-ext-docker-image.sh Release single 4

docker image push coherellc/haero-tpl:Debug-double-pack-size-1
docker image push coherellc/haero-tpl:Debug-double-pack-size-4
docker image push coherellc/haero-tpl:Release-double-pack-size-1
docker image push coherellc/haero-tpl:Release-double-pack-size-4
docker image push coherellc/haero-tpl:Debug-single-pack-size-1
docker image push coherellc/haero-tpl:Debug-single-pack-size-4
docker image push coherellc/haero-tpl:Release-single-pack-size-1
docker image push coherellc/haero-tpl:Release-single-pack-size-4
