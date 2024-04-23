%build cauchy matrices with formula: 
ridx = (1:m).'; cidx = 1:n;

Cout = ww.^(cidx).*(expm1(1i*lam(ridx)*n))./expm1(1i*lam(ridx)-2*pi*1i*(cidx)/n);
Cout = Cout; 


isgood = 1-(nodes(ridx)).*(ww.^(2*(cidx)));
isgood = abs(isgood) < 1e-9; 
[cidxm, ridxm] = meshgrid(cidx, ridx);

Cout(isgood) = ww.^(cidxm(isgood)).*(1 - ( nodes(ridxm(isgood)).*ww.^(2*cidxm(isgood))).^n)./...
            (1-nodes(ridxm(isgood)).*(ww.^(2*cidxm(isgood))));

Cout = Cout/sqrt(n);