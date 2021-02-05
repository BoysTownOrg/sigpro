function tst_bbf
[b,a]=butter(2,[0.029462,0.058925]);
fprintf('b =');
fprintf(' %8.4f',b*1e6);
fprintf('\n');
fprintf('a =');
fprintf(' %8.4f',a);
fprintf('\n');
return
