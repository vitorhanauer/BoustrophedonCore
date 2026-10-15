#!/bin/sh
set -eu

image=${BP_BENCH_IMAGE:-boustrophedon-bench}
repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)

printf 'Building %s from %s\n' "$image" "$repo_root" >&2
docker build --tag "$image" "$repo_root" >&2

printf 'Host: ' >&2
uname -a >&2
printf 'Host CPU model: ' >&2
awk -F: '/model name/ { sub(/^ /, "", $2); print $2; exit }' /proc/cpuinfo >&2
printf 'Host online CPUs: ' >&2
getconf _NPROCESSORS_ONLN >&2
printf 'Image: ' >&2
docker image inspect --format '{{.Id}}' "$image" >&2
printf 'Container UID and limits (memory.max, memory.swap.max, cpu.max):\n' >&2
docker run --rm \
    --memory=512m \
    --memory-swap=512m \
    --cpus=1 \
    --entrypoint /bin/sh \
    "$image" -c 'id -u; cat /sys/fs/cgroup/memory.max /sys/fs/cgroup/memory.swap.max /sys/fs/cgroup/cpu.max' >&2

if [ "$#" -eq 0 ]; then
    set -- --all
fi

exec docker run --rm \
    --memory=512m \
    --memory-swap=512m \
    --cpus=1 \
    "$image" "$@"
