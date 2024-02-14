set term postscript color enhanced eps 25

set out "Fig_R_Fukushima.eps"
set size 1,1.2

set title "e^+ e^- production rate, Leading Landau Level approx\n\n (Fukushima) M^2=(200)^2 [MeV^2], B=10 M_{pi}^2"
set xlabel "P_T [MeV]"
set ylabel "Const R/dQ^4" offset 2,0
set key bottom
set mxtics 10
set mytics 10

plot [:][0.5:1.5]\
"B10MPI_FUKUSHIMA/epem_perp.dat"     using 2:3 title "p_T in x-y plane (perp to z = B)" w l lt 1 lc 1 lw 3,\
"B10MPI_FUKUSHIMA/epem_parallel.dat" using 2:3 title "P_T in z-y plane (para to z = B)" w l lt 1 lc 2 lw 3

set out "Fig_R_N1.eps"
set size 1,1.2

set title "e^+ e^- production rate, N_1 only approx\n\n M^2=(200)^2 [MeV^2], B=10 M_{pi}^2"
set xlabel "P_T [MeV]"
set ylabel "Const R/dQ^4" offset 2,0
set key bottom
set mxtics 10
set mytics 10

plot [:][0.5:1.5]\
"B10MPI_N1/epem_perp.dat"     using 2:3 title "p_T in x-y plane (perp to z = B)" w l lt 1 lc 1 lw 3,\
"B10MPI_N1/epem_parallel.dat" using 2:3 title "P_T in z-y plane (para to z = B)" w l lt 1 lc 2 lw 3

set out "Fig_R_N1N0.eps"
set size 1,1.2

set title "e^+ e^- production rate, N_1+N_0 only approx\n\n M^2=(200)^2 [MeV^2], B=10 M_{pi}^2"
set xlabel "P_T [MeV]"
set ylabel "Const R/dQ^4" offset 2,0
set key bottom
set mxtics 10
set mytics 10

plot [:][0.5:1.5]\
"B10MPIN1N0/epem_perp.dat"     using 2:3 title "p_T in x-y plane (perp to z = B)" w l lt 1 lc 1 lw 3,\
"B10MPIN1N0/epem_parallel.dat" using 2:3 title "P_T in z-y plane (para to z = B)" w l lt 1 lc 2 lw 3

set out "Fig_R_FULL.eps"
set size 1,1.2

set title "e^+ e^- production rate, Full contribution\n\n M^2=(200)^2 [MeV^2], B=10 M_{pi}^2"
set xlabel "P_T [MeV]"
set ylabel "Const R/dQ^4" offset 2,0
set key bottom
set mxtics 10
set mytics 10

plot [:][0.5:1.5]\
"B10MPI_FULL/epem_perp.dat"     using 2:3 title "p_T in x-y plane (perp to z = B)" w l lt 1 lc 1 lw 3,\
"B10MPI_FULL/epem_parallel.dat" using 2:3 title "P_T in z-y plane (para to z = B)" w l lt 1 lc 2 lw 3

#B10MPIN1N0
#B10MPI_FUKUSHIMA
#B10MPI_FULL
#B10MPI_N1
#Fig_R_Fukushima.eps
#epem_pair_prod_asym_test
#makefig.plt
#makefig.plt~
