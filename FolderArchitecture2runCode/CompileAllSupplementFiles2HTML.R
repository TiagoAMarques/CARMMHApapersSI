library(knitr)
library(rmarkdown)

for(i in 1:8){
  #create all the html files
render(input=paste0("Supplement_S",i,".Rmd"),output_format=c("html_document"))
  #create all the pdf files
  # does not work - why - says '"pdflatex"' not found
  # render(input=paste0("Supplement_S",i,".Rmd"),output_format=c("html_document"))
  
  }