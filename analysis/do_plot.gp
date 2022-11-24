set term postscript enhanced solid color
set output 'TRAJs.ps'
set encoding iso_8859_1

E0=-0.1228747317E+004
unset key
ntrajs=40

set xlabel "Time (fs)"
set title "Distance C-S"
set ylabel "Distance ({\305})"
plot "TRAJ0/distance.time" u ($0*.5):1 w l lt 1 lw 3,\
for [i=1:ntrajs] "TRAJ".i."/distance.time" u ($0*.5):1 w l lt 3
set title "Singlet population"
set ylabel "Singlet probability"
plot [][-.01:1.01] "TRAJ0/pop_d.out" u 1:(sum [i=2:5] column(i)) w l lw 3,\
for [i=1:ntrajs] "TRAJ".i."/pop_d.out" u 1:(sum [i=2:5] column(i)) w l lt 3

set ylabel "Energy (eV)"
do for [i=0:ntrajs] {
set title "TRAJ".i.""
p [][-.1:16] "TRAJ".i."/energy.out" u 1:($2+$3-E0)*27.211 w l lw 3,"" u 1:($3-E0)*27.211 w p ps 2,\
"" u 1:($4-E0)*27.211 w l lt 3,\
"" u 1:($5-E0)*27.211 w l lt 3,\
"" u 1:($6-E0)*27.211 w l lt 3,\
"" u 1:($7-E0)*27.211 w l lt 3,\
"" u 1:($8-E0)*27.211 w l lt 3,\
"" u 1:($9-E0)*27.211 w l lt 3,\
"" u 1:($10-E0)*27.211 w l lt 3,\
"" u 1:($11-E0)*27.211 w l lt 3,\
"" u 1:($12-E0)*27.211 w l lt 3,\
"" u 1:($13-E0)*27.211 w l lt 3,\
"" u 1:($14-E0)*27.211 w l lt 3,\
"" u 1:($15-E0)*27.211 w l lt 3,\
"" u 1:($16-E0)*27.211 w l lt 3,\
"" u 1:($17-E0)*27.211 w l lt 3,\
"" u 1:($18-E0)*27.211 w l lt 3,\
"" u 1:($19-E0)*27.211 w l lt 3
}

set ylabel "Diabatic Populations"
set key
do for [i=0:ntrajs] {
set title "TRAJ".i.""
p [][-.01:1.01] \
"TRAJ".i."/pop_d.out" u 1:2 title "S0" w l,\
"TRAJ".i."/pop_d.out" u 1:3 title "S1" w l,\
"TRAJ".i."/pop_d.out" u 1:4 title "S2" w l,\
"TRAJ".i."/pop_d.out" u 1:5 title "S3" w l,\
"TRAJ".i."/pop_d.out" u 1:($6+$7+$8) title "T1" w l,\
"TRAJ".i."/pop_d.out" u 1:($9+$10+$11) title "T2" w l,\
"TRAJ".i."/pop_d.out" u 1:($12+$13+$14) title "T3" w l,\
"TRAJ".i."/pop_d.out" u 1:($15+$16+$17) title "T4" w l
}
