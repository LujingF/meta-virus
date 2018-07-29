lefse_path/format_input.py for_lefse.matrix lefse_input.txt.in -c 1 -o 1000000
lefse_path/lefse/run_lefse.py lefse_input.txt.in lefse_input.txt.in.res
lefse_path/plot_cladogram.py lefse_input.txt.in.res lefse_input.txt.in.res.caladogram.png --format png --dpi 400
awk '{if(NF==5) print}' lefse_input.txt.in.res | grep 'S__' | awk '{print $1}' | awk -F '.' '{print $NF}' | sort | uniq > lefse_diff_species.result
