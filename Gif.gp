set terminal gif animate delay 10
set output 'Correlator.gif'

FILE = 'Correlator.dat'
stats FILE nooutput
set xrange [*: *]
set yrange [0.5:1]
set ylabel "Sz"
set xlabel "x"
percentile="P5 P10 P20 P25 P50 P75"
cir=7;

f(x) = x;
max_y = 0.5;
do for [i=1:int(STATS_blocks)] {
    stats FILE index (i) 
    #print(STATS_max_y);
    a = STATS_max_y;
    plot FILE index (i) u ($1):($3*$3 + $2*$2) w lp pt cir ps 1.5 lt rgb "red" title columnheader , a title sprintf("max %f" , STATS_max_y)
}


set terminal gif animate delay 1
set output 'Sz_profile_avergage.gif'
stats 'Sz_average_profile.dat' nooutput
set xrange [-1:1]
set yrange [-0.6:0.6]
set ylabel "Sz avrg"
set xlabel "x/t"
percentile="P5 P10 P20 P25 P50 P75"
cir=7;

f(x) = x;
do for [i=1:int(STATS_blocks)] {
    plot "Sz_average_profile.dat" index (i) u ($1/(f(i))):($2) w lp pt cir ps 1.5 lt rgb "red" title columnheader
}

set terminal gif animate delay 1
set output 'Sz_profile_avergage_sqrt.gif'
stats 'Sz_average_profile.dat' nooutput
set xrange [-10:10]
set yrange [-0.6:0.6]
set ylabel "Sz avrg"
set xlabel "x/sqrt t"
percentile="P5 P10 P20 P25 P50 P75"
cir=7;

f(x) = x**(0.5);
do for [i=1:int(STATS_blocks)] {
    plot "Sz_average_profile.dat" index (i) u ($1/(f(i))):($2) w lp pt cir ps 1.5 lt rgb "red" title columnheader
}


set terminal gif animate delay 1
set output 'Sz_profile_avergage_not_rescaled.gif'
stats 'Sz_average_profile.dat' nooutput
set xrange [-200:200]
set yrange [-0.6:0.6]
set ylabel "Sz avrg"
set xlabel "x"
percentile="P5 P10 P20 P25 P50 P75"
cir=7;

f(x) = x;
do for [i=1:int(STATS_blocks)] {
    plot "Sz_average_profile.dat" index (i) u ($1):($2) w lp pt cir ps 1.5 lt rgb "red" title columnheader
}



##, "" index (i) u ((abs($1)<=f(i))? $1 : NaN ):(0.0) w l lw 10 ti "x=t"

