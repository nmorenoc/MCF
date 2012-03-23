#! /bin/bash

if [ -z "$1" ]; then
    dirname=.
else
    dirname=$1
fi

if [ -z "$2" ]; then
    time=50
else
    time=$2
fi

echo "dir :" $dirname;
echo "time:" $time;

perl $HOME/MCF/mcf/scripts/boundary_adapt_t_show_step.perl ${dirname} ${time}
