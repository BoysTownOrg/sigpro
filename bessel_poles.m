function bessel_poles
fp = fopen('bessel_poles.txt','w');
for n=1:16
   [z,p,k]=besselap(n);
   p=sort(p);
   fprintf(fp,'static float p%02d[] = {',n);
   k = fix(n/2);
   for m=1:k
      pp = p(n - 2 * (m - 1));
      fprintf(fp,'%12.9f,%12.9f',real(pp),imag(pp));
      if (2*m < n)
         fprintf(fp,',');
      end
   end
   if rem(n,2)
      fprintf(fp,'%12.9f',real(p(1)));
   end
   fprintf(fp,'};\n');
end
fclose(fp);
