### set xrange [0:50]
### plot 'balance.txt' using ($1*179.5):2 with lines,\
###      'balance.txt' using ($1*179.5):3 with lines,\
###      'balance.txt' using ($1*179.5):4 with lines,\
###      'balance.txt' using ($1*179.5):5 with lines,\
###      'mkm' using 2:3 with points,\
###      'mkm' using 2:4 with points,\
###      'mkm' using 2:7 with points,\
###      'mkm' using 2:8 with points

### pause -1

plot 'rey_regions.txt' using ($1*179.5):(-$2*4200/179.5*4200/179.5) with lines,\
     'rey_regions.txt' using ($1*179.5):(-$3*4200/179.5*4200/179.5) with lines,\
     'rey_regions.txt' using ($1*179.5):(-$4*4200/179.5*4200/179.5) with lines,\
     'rey.txt' using ($1*179.5):(-$2*4200/179.5*4200/179.5) with points

pause -1

### set autoscale
### plot 'velbar.txt' using 1:2 with lines,\
###      'velbar.txt' using 1:3 with lines,\
###      'velbar.txt' using 1:4 with lines
### pause -1

### plot 'vortbar.txt' using 1:2 with lines,\
###      'vortbar.txt' using 1:3 with lines,\
###      'vortbar.txt' using 1:4 with lines
### pause -1

### plot 'velrms.txt' using 1:2 with lines,\
###      'velrms.txt' using 1:3 with lines,\
###      'velrms.txt' using 1:4 with lines
### pause -1

### plot 'vortrms.txt' using 1:2 with lines,\
###      'vortrms.txt' using 1:3 with lines,\
###      'vortrms.txt' using 1:4 with lines
### pause -1

### plot 'vorttrms.txt' using 1:2 with lines,\
###      'vorttrms.txt' using 1:3 with lines,\
###      'vorttrms.txt' using 1:4 with lines,\
###      'vortrms.txt' using 1:2 with points,\
###      'vortrms.txt' using 1:3 with points,\
###      'vortrms.txt' using 1:4 with points
### pause -1

### plot 'drg.txt' using 1:4 with lines
### pause -1

### set xrange [0:179.5]
### plot 'veldrms.txt' using ($1*179.5):($2*4200/179.5) with lines,\
###      'veldrms.txt' using ($1*179.5):($3*4200/179.5) with lines,\
###      'veldrms.txt' using ($1*179.5):($4*4200/179.5) with lines
### pause -1

### plot 'veltrms.txt' using ($1*179.5):($2*4200/179.5) with lines,\
###      'veltrms.txt' using ($1*179.5):($3*4200/179.5) with lines,\
###      'veltrms.txt' using ($1*179.5):($4*4200/179.5) with lines
### pause -1

### plot 'rey.txt' using ($1*179.5):(-$2*4200/179.5*4200/179.5) with lines,\
###      'rey_dt.txt' using ($1*179.5):(-$2*4200/179.5*4200/179.5) with lines,\
###      'rey_dt.txt' using ($1*179.5):(-$3*4200/179.5*4200/179.5) with lines
### pause -1


