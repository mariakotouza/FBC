FROM rocker/r-ver:3.6.0

RUN mkdir /home/analysis
RUN mkdir /home/analysis/output

RUN apt-get update
RUN apt-get install libxml2-dev -y
RUN apt-get install r-cran-xml -y

RUN R -e "install.packages(c('parallel','plyr','pryr','dplyr','data.table','stringr','tidyr','entropy','ggplot2','ggseqlogo','gridExtra','cluster','seqinr','collapsibleTree','data.tree','DiagrammeR','stringdist','igraph','networkD3','plsgenomics','xtable','tictoc','optparse'), repos='http://cran.rstudio.com/')"

COPY /app/. /home/analysis/.