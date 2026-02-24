# Base image https://hub.docker.com/
FROM r-base

# System libraries of general use
## install debian packages
RUN apt-get update -qq && apt-get -y --no-install-recommends install \	
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
	libtool 

# Create the Shiny app directory in the container
RUN mkdir ./app

# Set permissions for the app directory
RUN chown -R 1001:1001 /app

# Install any needed packages
RUN R -e "install.packages(c('httpuv', 'shiny', 'shinyDashboardThemeDIY', 'dashboardthemes', 'shinydashboard', \
	'shinyvalidate', 'DT', 'ggiraph', 'tidyverse', 'markdown', 'ggnewscale', 'khroma', 'units', \
	'sf', 'spdep', 'sfsmisc'))"

# Copy the app files to the container
cp ./app ./app

WORKDIR /app

# Expose the port that the app will be running on
EXPOSE 3838

# Run the app on container start
CMD ["R", "-e", "shiny::runApp('.', host = '0.0.0.0', port = 3838)"]

