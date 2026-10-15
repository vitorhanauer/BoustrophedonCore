FROM debian:bookworm-slim AS build

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update \
    && apt-get install --no-install-recommends -y \
        cmake \
        gcc \
        libc-dev \
        make \
        time \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /src
COPY . .

RUN cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON \
    && cmake --build build --parallel 2 \
    && ctest --test-dir build --output-on-failure \
    && cmake --install build --prefix /usr/local \
    && test -x build/bp_runner \
    && mkdir -p /artifacts \
    && install -m 0755 build/bp_runner /artifacts/bp_runner \
    && if [ -x build/bp_example ]; then install -m 0755 build/bp_example /artifacts/bp_example; fi

FROM debian:bookworm-slim AS runtime

COPY --from=build /usr/local/lib/ /usr/local/lib/
COPY --from=build /usr/local/include/boustrophedon.h /usr/local/include/boustrophedon.h
COPY --from=build /artifacts/ /usr/local/bin/

RUN ldconfig

USER 10001:10001
ENTRYPOINT ["/usr/local/bin/bp_runner"]
CMD ["--all"]
