FROM r-base

RUN apt-get update -qq && \
    apt-get -y --no-install-recommends install \
        libfontconfig1-dev \
        libcurl4-openssl-dev \
        gcc \
        g++ \
        libfreetype6-dev \
        libglib2.0-dev \
        libcairo2-dev \
        libharfbuzz-dev \
        libfribidi-dev \
        libtiff5-dev \
        libxml2-dev \
        libudunits2-dev \
        libgdal-dev \
        libproj-dev \
        libgeos-dev \
        libtool && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN R -e "install.packages(c(
    'httpuv','shiny','shinyDashboardThemeDIY','dashboardthemes',
    'shinydashboard','shinyvalidate','DT','ggiraph',
    'tidyverse','markdown','ggnewscale','khroma',
    'units','sf','spdep','sfsmisc'
), repos='https://cloud.r-project.org')"

WORKDIR /app
COPY ./app /app

EXPOSE 3838

CMD ["R", "-e", "shiny::runApp('.', host='0.0.0.0', port=3838)"]
