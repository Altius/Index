FROM hotspot2
RUN apk add --no-cache \
  g++ \
  perl \
  R \
  R-dev
RUN Rscript -e 'install.packages("caTools", repo="http://cran.rstudio.com")'

WORKDIR masterlist

ADD . .
