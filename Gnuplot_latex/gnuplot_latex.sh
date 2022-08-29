function gnuplot_latex() {
  gnuplot -e ${3:-'x=0'}"; texoutput='"${2:-$1}"'" $1.plt
  sed -i 's/calc/calc,physics,siunitx,mathabx,color/g' ${2:-$1}.tex
  pdflatex -interaction=nonstopmode ${2:-$1}.tex
  rm ${2:-$1}-inc.eps ${2:-$1}-inc-eps-converted-to.pdf ${2:-$1}.aux ${2:-$1}.log ${2:-$1}.log ${2:-$1}.tex
}