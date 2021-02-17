set terminal gif animate delay 5
set output 'Correlator.gif'

FILE = 'Correlator20.dat'
stats FILE nooutput
set xrange [-2: 2]
set yrange [-2:20]
set ylabel "Abs"
set xlabel "x"
percentile="P5 P10 P20 P25 P50 P75"
cir=7;


do for [i=1:int(STATS_blocks)] {
    stats FILE index (i) 
    #print(STATS_max_y);
    #a = STATS_max_y;
    plot FILE index (i) u ($1):($3*$3 + $2*$2) w lp pt cir ps 1.5 lt rgb "red" title columnheader, \
    FILE index (i) u ($1):($2) w lp pt cir ps 1.5 lt rgb "blue" title "Re", FILE index (i) u ($1):($3) w lp pt cir ps 1.5 lt rgb "green" title "Re"
}




