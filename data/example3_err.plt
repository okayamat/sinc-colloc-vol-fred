set encoding iso_8859_1
set term postscript eps enhanced "Times-Roman" 20
set output "example3_err.eps"
set size 0.9
set logscale x
set logscale y
set key spacing 3
set xrange [10:1000]
set yrange [1e-16:1]
set format y "10^{%L}"
set xlabel "{/Times-Italic=24 n}{/Times-Roman=24 =2}{/Times-Italic=24 N}{/Times-Roman=24 +1}"
set ylabel "{/Times-Roman=24 maximum error}"
plot "SE_orig_ex3.dat" using 1:2 w lp title "Original Sinc-collocation (SE)" pt 2 ps 1, "SE_new_ex3.dat" using 1:2 w lp title "New Sinc-collocation (SE)" pt 8 ps 1, "DE_orig_ex3.dat" using 1:2 w lp title "Original Sinc-collocation (DE)" pt 4 ps 1, "DE_new_ex3.dat" using 1:2 w lp title "New Sinc-collocation (DE)" pt 6 ps 1, "tabella_nystrom.dat" using 1:2 w lp title "product Gauss rules" pt 1 ps 1
