set sty dat l
set xrange [0:1029]
set yrange [-2:3]
set xlabel 'parcel index'
set ylabel 'scalar value'
plot 'Data_00000.dat' us 1 ti "Sc=1/16", 'Data_00000.dat' us 2 ti "Sc=1", 'Data_00000.dat' us 3 ti "Sc=16"; pause -1
plot 'Data_00001.dat' us 1 ti "Sc=1/16", 'Data_00001.dat' us 2 ti "Sc=1", 'Data_00001.dat' us 3 ti "Sc=16"; pause -1
plot 'Data_00002.dat' us 1 ti "Sc=1/16", 'Data_00002.dat' us 2 ti "Sc=1", 'Data_00002.dat' us 3 ti "Sc=16"; pause -1
plot 'Data_00003.dat' us 1 ti "Sc=1/16", 'Data_00003.dat' us 2 ti "Sc=1", 'Data_00003.dat' us 3 ti "Sc=16"; pause -1
plot 'Data_00004.dat' us 1 ti "Sc=1/16", 'Data_00004.dat' us 2 ti "Sc=1", 'Data_00004.dat' us 3 ti "Sc=16"; pause -1
