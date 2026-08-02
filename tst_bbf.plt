; tst_bbf.plt - plot Butterworth band-pass filter results
;
head=n : ticdir=in : yhor=y : clip=y
pltyp=lines : sizfac=6
xmin=0.01 : xmax=1 : xcyc=2 : xper=90
ymin=-40  : ymax=0 : yint=4 : yper=66.7
zdat=0
ylen=3.0 : yllc=4.5
xanskp=-1 : ydata=$2
;
select=$2>-90
xlabel=
ylabel=magnitude (dB)
pltcol=1
include tst_bbf.txt
plot
0.5 1.2 "
    fourth-order
|_0,pltcol=1| Butterworth
"
;
newframe
ymin=0 : ymax=50 : yint=5 : yper=100*15/17
ylen=3.0 : yllc=1.5 : yfmt=i
xanskp=0 : ydata=$4
xlabel=frequency / half_sample_rate
ylabel=delay (samples)
pltcol=1
include tst_bbf.txt
plot
