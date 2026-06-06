FROM rocker/r-ver:4.3.0

# Install underlying Linux compilation tools and dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    libxml2-dev libssl-dev libcurl4-openssl-dev libgsl-dev libglpk-dev \
    libpng-dev zlib1g-dev make g++ \
    && rm -rf /var/lib/apt/lists/*

RUN R -e "options(repos = c(CRAN = 'https://cloud.r-project.org')); install.packages('BiocManager')"
RUN R -e "BiocManager::install(c('Seurat', 'ggplot2', 'dplyr', 'Matrix', 'patchwork', 'limma', 'gprofiler2', 'RColorBrewer', 'scales'))"

WORKDIR /workspace

CMD ["Rscript", "main.R"]