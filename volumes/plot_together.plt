#!/bin/bash


# ------- MERGE ALL TXT FILES INTO A SINGLE ONE ------
OUT="together.txt"

: > "$OUT"                              # start with an empty output
for f in $(ls -1t -- *.txt); do                      # loop through all .txt in the cwd
    [[ "$f" == "$OUT" ]] && continue    # skip the output itself if it lives here
    cat "$f" >> "$OUT"                  # add exactly two blank lines after each file
    echo >> "$OUT"
    echo >> "$OUT"
done


gnuplot <<- EOF
        # ---------------------- VOLUME DISTRIBUTION ------------------------ 
        set terminal pngcairo enhanced font 'Arial,18' size 1080,1080 enhanced

        # --- STATS ---
        unset xrange
        unset yrange
        stats "$OUT" u 5 name "V" nooutput

        # --- BIN FUN ---
        # binwidth = 0.2;
        binwidth = (ceil(V_max) - floor(V_min))/50.0;
        bin(x,width)=width*floor(x/width) + width/2.0


        set output "together_0.png"

        set title "Distribution of cell volumes"
        set xlabel 'Volume (cm^3)'
        set ylabel 'Frequency'

        # set xtics 1
        set grid
        set boxwidth binwidth
        set style fill empty border -1
        set xrange[0:ceil(V_max)]
        plot '$OUT' using (bin(\$5,binwidth)):(1.0/V_records) i 0 smooth freq with boxes lw 3 lc rgb "#0072BD" title '3 year old' , \
             '$OUT' using (bin(\$5,binwidth)):(1.0/V_records) i 1 smooth freq with boxes lw 3 lc rgb "#D95319" title '8 year old' , \
             '$OUT' using (bin(\$5,binwidth)):(1.0/V_records) i 2 smooth freq with boxes lw 3 lc rgb "#EDB120" title '13 year old', \
             '$OUT' using (bin(\$5,binwidth)):(1.0/V_records) i 3 smooth freq with boxes lw 3 lc rgb "#7E2F8E" title 'Adult'



        # ---------------------- VOLUME VS Z ------------------------ 

         # --- STATS ---
        unset xrange
        unset yrange
        stats '$OUT' u 4 i 0 nooutput name "Z0"
        stats '$OUT' u 4 i 1 nooutput name "Z1"
        stats '$OUT' u 4 i 2 nooutput name "Z2"
        stats '$OUT' u 4 i 3 nooutput name "Z3"

        # --- FIT ---
        f0(x) = m0*x + b0; m0 = 1; b0 = 1
        f1(x) = m1*x + b1; m1 = 1; b1 = 1
        f2(x) = m2*x + b2; m2 = 1; b2 = 1
        f3(x) = m3*x + b3; m3 = 1; b3 = 1

        set fit quiet
        set fit logfile '/dev/null'

        fit f0(x) "$1" using 4:5 i 0 via m0,b0 
        fit f1(x) "$1" using 4:5 i 1 via m1,b1 
        fit f2(x) "$1" using 4:5 i 2 via m2,b2 
        fit f3(x) "$1" using 4:5 i 3 via m3,b3 

        set output "together_1.png"

        set xtics 2
        # set xrange[-26:-2]
        # set yrange[0:8]
        set title "Volume and height"
        set xlabel 'z coordinate of seed (cm)'
        set ylabel 'Volume of Voronoi cell (cm^3)'

        c0 = "#0072BD"
        c1 = "#D95319"
        c2 = "#EDB120"
        c3 = "#7E2F8E"
        
        plot \
            "$OUT" i 0 using 4:5 with points lw 0.2 lt 6 lc rgb c0 notitle, \
            [Z0_min:Z0_max] f0(x) with lines lw 4 lc rgb c0 title "3 yar old", \
            "$OUT" i 1 using 4:5 with points lw 0.1 lt 6 lc rgb c1 notitle, \
            [Z1_min:Z1_max] f1(x) with lines lw 4 lc rgb c1 title "8 year old", \
            "$OUT" i 2 using 4:5 with points lw 0.1 lt 6 lc rgb c2 notitle, \
            [Z2_min:Z2_max] f2(x) with lines lw 4 lc rgb c2 title "13 year old", \
            "$OUT" i 3 using 4:5 with points lw 0.1 lt 6 lc rgb c3 notitle, \
            [Z3_min:Z3_max] f3(x) with lines lw 4 lc rgb c3 title "Adult"

EOF
