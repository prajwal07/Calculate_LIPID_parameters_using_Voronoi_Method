set terminal jpeg size 2050,1000
set output 'Analysis_APL_vernoi_plot_frame_1.jpg'
set multiplot layout 1,2 rowsfirst

set title "TOP LAYER" font "Helvetica,20"
plot 'vtm-analysis/frame_1_layer1_pro.crd', 'vtm-analysis/frame_1_layer1_lip.crd', 'vtm-analysis/frame_1_layer1_rnd.crd', 'vtm-analysis/frame_1_layer1_b.vor' with lines, 'vtm-analysis/frame_1_layer1_nb.vor' with lines

set title "BOTTOM LAYER" font "Helvetica,20"
plot 'vtm-analysis/frame_1_layer2_pro.crd', 'vtm-analysis/frame_1_layer2_lip.crd', 'vtm-analysis/frame_1_layer2_rnd.crd', 'vtm-analysis/frame_1_layer2_b.vor' with lines, 'vtm-analysis/frame_1_layer2_nb.vor' with lines
