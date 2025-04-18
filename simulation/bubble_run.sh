#!/bin/bash

ns=15

for i in $ns # 2 4 6 8 10 12
do echo $i
   python3 bubble_moc.py $i
done

# python3 bubble_plot_solution.py $ns
# python3 bubble_plot_error.py 
