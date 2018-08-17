function yout = pdf_fit(xcp,ycp,Nbins,xstr,ystr)


x_var = linspace(min(xcp),max(xcp),101);
y_var = linspace(min(ycp),max(ycp),101);


sigma_x = sqrt(var(xcp));
sigma_y = sqrt(var(ycp));
mean_x = mean(xcp);
mean_y = mean(ycp);

norm_x = normpdf(x_var,mean_x,sigma_x);
norm_y = normpdf(y_var,mean_y,sigma_y);


figure,
subplot(211)
hx = histogram(xcp,Nbins);
fct_x = trapz(hx.BinEdges(1:end-1),hx.Values);

bar(hx.BinEdges(1:end-1),hx.Values./fct_x);
hold on;
plot(x_var,norm_x,'r-');
xlabel(xstr); ylabel('PDF of x');


subplot(212)
hy = histogram(ycp,Nbins);
fct_y = trapz(hy.BinEdges(1:end-1),hy.Values);

bar(hy.BinEdges(1:end-1),hy.Values./fct_y);
hold on;
plot(y_var,norm_y,'r-');
xlabel(ystr); ylabel('PDF of y');
