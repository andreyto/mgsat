#!/bin/bash

## This will build all MGSAT and derived Docker images
## The Dockerfile is a multi-stage file, and therefore
## will not work on RHEL 7 olden Docker. It needs
## the newer one from docker.com.
## In dire need, you can split the Dockerfile into
## several single-stage files and build separately.
## Set envar MGSAT_REGISTRY if you are going to push
## the resulting images to the registry.

this_dir=$(cd $(dirname $0); pwd)

MGSAT_TAG=2.8
if [ -n "$MGSAT_REGISTRY" ]; then
    registry_pref=${MGSAT_REGISTRY}/
else
    registry_pref=""
fi
MGSAT_DOCKERFILE="${this_dir}/Dockerfile.mgsat_rocker_ver"

docker_build() {
    target=$1
    docker build -f "$MGSAT_DOCKERFILE" \
        --tag $registry_pref$target:latest \
        --tag $registry_pref$target:$MGSAT_TAG \
        --target $target .
}

rm -rf mgsat_build_context
mkdir -p mgsat_build_context

pushd mgsat_build_context

docker_build mgsat-deps
docker_build mgsat
docker_build cgfease
docker_build multiomig

popd
