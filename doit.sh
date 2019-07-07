#!/bin/sh

#pdflatex -interaction nonstopmode -halt-on-error -file-line-error thesis
#bibtex thesis 
#pdflatex -interaction nonstopmode -halt-on-error -file-line-error thesis
#pdflatex -interaction nonstopmode -halt-on-error -file-line-error thesis

#open -a Preview thesis.pdf

pdflatex -interaction scrollmode thesis 
bibtex thesis 
pdflatex -interaction scrollmode thesis 
pdflatex -interaction scrollmode thesis 
open -a Preview thesis.pdf

#latex thesis && \
#biber thesis && \
#latex thesis

# with biblatex, it is only necessary to run latex once after bibtex, not twice.
# with older tex installations, you may need to change 'biber' to 'bibtex'
