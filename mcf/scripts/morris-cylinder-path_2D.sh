awk '$1<=0.052&&$1>=0.048{print} !NF{print}' particles.dat > path1.dat
awk '$1>=0.096{print} !NF{print}' particles.dat > path2.dat
awk '(($1<0.03||$1>0.07)&&($2>=0.05&&$2<=0.054)) ||( $2>=0.05 && ($1>=0.03&&$1<=0.07) && (($1-0.05)**2+($2-0.05)**2)>=0.02**2 && (($1-0.05)**2+($2-0.05)**2)<=0.024**2 ){print} !NF {print}' particles.dat > path3.dat
awk '$2>=0.096{print} !NF{print}' particles.dat > path4.dat