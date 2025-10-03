##
# Globals
#

# load packages
require('shinydashboard')
require('shinyvalidate')

# Load helper functions
files = list.files(path="./R", pattern="\\.R$", full.names=TRUE, recursive=TRUE)
sapply(files, source)
