function [xc, yc] = drift_adj(x,y)

% Mohammad Asif Zaman 
% Aug. 15, 2018

% Input arguments: tracked x any y position of the practicle with drift
% Output argument: drift corrected x and y position of the particle


%  Frame indices

nf = 1:length(x);

% Linear fit coefficients
px = polyfit(nf,x,1);  
py = polyfit(nf,y,1);

xc = x - polyval(px,nf);
yc = y - polyval(py,nf);

