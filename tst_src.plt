; tst_src.plt

head=n : yhor=y : ticdir=in : clip=y
sizfac=5 : 
xmin=0 : xmax=1 : xint=10 : xper=90
ymin=-1 : ymax=1 : yint=8 : yper=90
zdat=0
xllc=1 : xlen=8.5 : ylen=6.5
yllc=1 : xanskp=-1 : yanskp=-1 : axlwt=0
xlab=
ylab=
tlab=tst_src: sample-rate conversion
;
pltyp=line 
include s0.txt
plot
pltyp=symb : symbol=0 : SYMsiz=1
include s1.txt
plot
pltyp=symb : symbol=11 : symsiz=0.5
include s2.txt
plot
mhkey=2.3
4  -0.0 "|_0| original signal"
mvkey=0.7
4  -0.2 " |00|  down-sampled signal: 3 samples/cycle"
mvkey=0.5
4  -0.4 " |11|  up-sampled signal: 18.2 samples/cycle"
