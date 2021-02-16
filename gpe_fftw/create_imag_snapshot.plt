set terminal png size 600,600
set output "imag_snapshot.png"
set pm3d map
set palette defined ( 0 '#000090',1 '#000fff',2 '#0090ff',3 '#0fffee',4 '#90ff70',5 '#ffee00',6 '#ff7000',7 '#ee0000',8 '#7f0000')
NX = 100
set title "Static GPE Solution"
stats "./data/latest/flux_imag.bin" binary record=(NX/5,-1) format="%*int%6double%*int" using 4:5 nooutput
scale = sqrt(STATS_max_x**2 + STATS_max_y**2)

splot "./data/latest/wf_imag_fin.bin" binary record=(NX,-1) format="%*int%6double%*int" using 1:2:4 with pm3d t "", \
"./data/latest/flux_imag.bin" binary record=(NX/5,-1) format="%*int%6double%*int" using 1:2:(0):($4/scale):($5/scale):(0) with vectors t "" lc rgb "black"

print "Saved as 'imag_snapshot.png'"
