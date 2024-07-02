set term png # setting terminal type

# Plot 1
set out "z1.png" # setting output file name for plot 1
set xl "t" # x-axis label for plot 1
set yl "x(t)" # y-axis label for plot 1
set title "Wychylenie x(t) - Przypadek 1" # plot title for plot 1
set grid # enabling grid for plot 1
plot "output1.txt" u 1:2 w p lt 3 pt 6 t ""

# Plot 2
set out "z2.png" # setting output file name for plot 2
set title "Wychylenie x(t) - Przypadek 2" # plot title for plot 2
plot "output2.txt" u 1:2 w p lt 4 pt 6 t ""

# Plot 3
set out "z3.png" # setting output file name for plot 3
set title "Wychylenie x(t) - Przypadek 3" # plot title for plot 3
plot "output3.txt" u 1:2 w p lt 5 pt 6 t ""
