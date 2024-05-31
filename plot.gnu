reset
set colorsequence classic
set terminal epslatex standalone color size 9cm,6cm
set out 'Nutrients.tex'
set xr[0:10]
set xtics 2
set yr[-0.01:1.01]
set ytics 0.2
set xlabel 'time $(\times 1000)$' offset 0,0.5
set ylabel '$N/N_0$' offset 0.5,0
set pointsize 2
set key top right samplen 2
plot "1_Disk-like/data.txt" u ($1*1000):($4/0.25) w l lc 0 dt 1 lw 2 title "Disk", "2_DBM-like/data.txt" u ($1*1000):($4/0.071) w l lc 1 dt 6 lw 2 title "DMB", "3_DLA-like/data.txt" u ($1*1000):($4/0.087) w l lc 2 dt 3 lw 2 title "DLA", "4_ConcentricRing/data.txt" u ($1*1000):($4/0.1) w l lc 3 dt 4 lw 2 title "Rings"

reset
set colorsequence classic
set terminal epslatex standalone color size 9cm,6cm
set out 'All.tex'
set xr[0:10]
set xtics 2
set yr[-2:202]
set ytics 50
unset key
set xlabel 'time $(\times 1000)$' offset 0,0.5
set ylabel '$(B+I)/B_0$' offset 0.5,0
set pointsize 2
plot "1_Disk-like/data.txt" u ($1*1000):($3/0.001254504) w l lc 0 dt 1 lw 2, "2_DBM-like/data.txt" u ($1*1000):($3/0.001254504) w l lc 1 dt 6 lw 2, "3_DLA-like/data.txt" u ($1*1000):($3/0.001254504) w l lc 2 dt 3 lw 2, "4_ConcentricRing/data.txt" u ($1*1000):($3/0.001254504) w l lc 3 dt 4 lw 2

# reset
# set colorsequence classic
# set terminal epslatex standalone color size 9cm,6cm
# set out 'Active.tex'
# set xr[0:51]
# set xtics 10
# set yr[0.01:61]
# set ytics 10
# set key top right samplen 2
# set xlabel 'time $(\times 1000)$' offset 0,0.5
# set ylabel '$B/B_0$' offset 0.5,0
# set pointsize 2
# plot "1_Disk-like/data.txt" u ($1*1000):($2/0.001254504) w l lc 0 dt 1 lw 2 title "Disk", "2_DBM-like/data.txt" u ($1*1000):($2/0.001254504) w l lc 1 dt 6 lw 2 title "DMB", "3_DLA-like/data.txt" u ($1*1000):($2/0.001254504) w l lc 2 dt 3 lw 2 title "DLA", "4_ConcentricRing/data.txt" u ($1*1000):($2/0.001254504) w l lc 3 dt 4 lw 2 title "Rings"
