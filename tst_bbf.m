function tst_bbf
gn=1e-3;
wn=[0.0295,0.0589];
fo=2;
fprintf('Butterworth band-pass: fo=%d wn=[%.4f %.4f]\n',fo,wn);
[b,a]=butter(fo,wn);
fprintf('b =');
fprintf(' %8.4f',b/gn);
fprintf(' * %.0e\n',gn);
fprintf('a =');
fprintf(' %8.4f',a);
fprintf('\n');
fprintf('Butterworth filter analog prototype: fo=%d\n',fo);
[z,p,k]=buttap(fo);
fprintf('Butterworth transfer function\n');
[num,den]=zp2tf(z,p,k);
return
