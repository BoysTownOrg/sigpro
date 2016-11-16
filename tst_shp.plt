; tst_shp.plt

head=n : yhor=y : ticdir=in : clip=y
sizfac=5 : pltyp=line 
xmin=0 : xmax=10 : xdat=$1/1000
ymin=-60 : ymax=20 : yint=8 : yper=80
zdat=0 : yanskp=1
xllc=2 : xlen=5 : ylen=3
;
newframe
yllc=4.5 : xanskp=-1
xlab=
ylab=spectral magnitude (dB)
lintyp=6
include spec1.txt
stat
n1=($n-1)*2
plot
lintyp=0
include spec2.txt
stat
n2=($n-1)*2
plot
%msg 0.5 0.52 "|_6| filter (n=%.0f)" n1
%msg 0.5 0.28 "|_0| waveform  (n=%.0f)" n2
;
newframe
yllc=1.5 : xanskp=0
xlab=frequency (kHz)
ylab=spectral magnitude (dB)
include spec3.txt
plot
%msg 0.5 0.4 "|_0| filtered waveform" n1
