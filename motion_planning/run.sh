#!/bin/bash

ns=15
nt=250

# solving initial triangle (run only initially)
python3 ./init/bubble_init.py $ns $ns
python3 ./init/bubble_moc.py $ns $ns

# solving characteristic steps
python3 bubble_init.py $ns $nt

for i in $ns # 2 4 6 8 10 12
do echo $i
   python3 bubble_moc.py 
done

cp -r ./data ./post/
# # # python3 bubble_plot.py
# # # python3 bubble_plot_error.py 
