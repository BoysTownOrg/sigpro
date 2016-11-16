; tst_tb.plt
sizfac=5 : yhor=y : ticdir=in : clip=y
xllc=1 : yllc=1 : xlen=9 : ylen=6
include tst_tb.txt
stat 
xmin=0 : xmax=ceil($x_max/100)*100
%repeat (rate>0)&(size>0)
xmax=1000*size/rate
xlab=(msec)
data
include tst_tb.txt
%%
xint=4 : xper=80
ymin=-1 : ymax=1 : yfmt=f.2
yint=8 : yper=90 : yanskp=1
pltyp=line
plot
%msg 0.5 (ylen+0.1) "tst_tb: rate=%.0f size=%.0f" rate size
