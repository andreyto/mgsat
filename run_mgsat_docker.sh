#!/bin/bash
## If you set MGSAT_DEV to any non-empty value, the current
## dir on the host will be bind-mounted to /home/rstudio/work.
## This is considered a development mode because mgsat directory
## in the image will be masked by whatever is in the current
## dir on the host (presumably the current dir on the host
## will contain checked-out mgsat subdirectory). Same for other
## code dirs under /home/rstudio/work in any derived containers.
##
## If MGSAT_DEV is not set, then the current host dir will be 
## bind-mounted into /home/work.
##
## WARNING: Rstudio startup script will chown everything under
## /home/rstudio to USERID, which is our script is set to the
## UID on the host. If you run with a set MGSAT_DEV, your
## entire current dir will be chowned. This might take a while
## and be not what you need. To avoid this, run with MGSAT_DEV
## unset.
set -e
echo -n "Set desired password for RStudio Web interface and press Enter: "
read -s MGSAT_DOCKER_PASSWD
if [ -n "$MGSAT_DEV" ]; then
    bind_mounts="-v $(pwd):/home/rstudio/work"
else
    bind_mounts="-v $(pwd):/home/work"
fi
echo
echo "Bind mounts are: $bind_mounts"
port=${MGSAT_DOCKER_PORT:-8787}
image=${MGSAT_DOCKER_IMAGE:-andreyto/mgsat:latest}
echo "Using Docker image: ${image}"
docker run --rm -p ${port}:8787 \
	-d -P \
	$bind_mounts \
	-e USERID=$UID -e PASSWORD="$MGSAT_DOCKER_PASSWD" \
	$image
echo "Point your Web browser to: $(hostname):${port}"
echo "Loging with user 'rstudio' and the password that you just entered"
