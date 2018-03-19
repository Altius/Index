# Build bedGraphToBigWig
FROM alpine:3.7 as kentutils-build
RUN apk add --no-cache \
      g++ \
      gcc \
      git \
      libpng-dev \
      make \
      mysql-dev \
      zlib-dev
# Get an archive of kentUtils, remove a file that doesn't build, and compile
RUN wget https://github.com/ENCODE-DCC/kentUtils/archive/v302.0.0.tar.gz \
      && tar xf v302.0.0.tar.gz \
      && cd kentUtils-302.0.0 \
      && sed -i 's/fof.o //' src/lib/makefile \
      && make

# Build masterlist container
FROM alpine:3.7
RUN apk add --no-cache \
  bash \
  g++ \
  R \
  R-dev
RUN Rscript -e 'install.packages("caTools", repo="http://cran.rstudio.com")'

# Get Bedops
RUN wget -O - https://github.com/bedops/bedops/releases/download/v2.4.31/bedops_linux_x86_64-v2.4.31.tar.bz2 \
  | tar -C /usr/local -xjf -

# Copy kentutils in
COPY --from=kentutils-build /kentUtils-302.0.0/bin/bedToBigBed /usr/local/bin/

ENV PATH=/masterlist:$PATH
ADD . /masterlist/
WORKDIR /data
