for F in *.f90; do
    cp $F ${F/.f90/.F90}
done