#!/bin/bash

# run: chmod +x compile_latex.sh
# run: ./compile_latex.sh Proposal.tex
# Check if the correct number of arguments is provided
if [ "$#" -ne 1 ]; then
  echo "Usage: $0 <filename.tex>"
  exit 1
fi

# Get the filename without the extension
filename=$(basename "$1" .tex)

# Run pdflatex, bibtex, and pdflatex twice
pdflatex "$filename.tex"
bibtex "$filename"
pdflatex "$filename.tex"
pdflatex "$filename.tex"

echo "Compilation finished. Check the output for any errors."