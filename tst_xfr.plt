; tst_xfr.plt - plot transfer-function results
;
head=n : ticdir=in : yhor=y : clip=y
pltyp=lines : sizfac=6
xmin=0.03 : xmax=30 : xcyc=2 : xper=100
ymin=-40  : ymax=0 : yint=4 : yper=66.7
zdat=0
ylen=3.0 : yllc=4.5
xanskp=-1 : ydata=$2
;
select=$2>-90
xlabel=
ylabel=magnitude (dB)
tlabel=transfer function
pltcol=1
include transfer.txt
plot
0.5 1.2 "
    fourth-order, low-pass
|_0,pltcol=1| Butterworth
"
;
newframe
ymin=-1 : ymax=0 : yint=2.2 : yper=100*2/3
ylen=3.0 : yllc=1.5 : yfmt=f.1
xanskp=0 : ydata=$3
xlabel=frequency (kHz)
ylabel=phase (cyc)
tlabel=
pltcol=1
include transfer.txt
plot
