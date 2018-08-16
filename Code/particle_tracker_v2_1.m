clear all;
close all;
clc;
clf;

% March 8, 2017: Mohammad Asif Zaman
% First draft of particle traker program. It has been tested for single
% particle tracking.


% Aug. 15, 2018
% The current version works with pixels rather than actual dimensions. The
% results can be converted to lenght units by multiplying with dx and dy
% and the end. It is better to keep the intemediate variables in pixle
% units.


% Control variables
Ctrl_drift_adjust = 'y';
Ctrl_centering = 'y';
Ctrl_video_playback = 'n';
Nbins = 12; % number of bins in the histogram


% Indicate the path of the video file

source_file = 'Video/test1.avi';
% source_file = 'Video\Z#7C70N2c3_D500_P196b.avi';



% Read video file
videoObj = VideoReader(source_file);
vidFrames = read(videoObj);

% Identifying frame size
[frameSizeH, frameSizeW, numColors, numFrames] = size(vidFrames);

% Initializing output frames
iFrames = zeros(frameSizeH, frameSizeW, numFrames);


% Dimension parameters
dx = (10e-6/196)*1e9;  % Size of one pixel
dia_guess = 25; % Guess of the diameter in pixels.


nf = 1:numFrames;  % Frame indices

tic
for m = 1:numFrames
    iFrames(:,:,m) = rgb2gray(vidFrames(:,:,:,m));  % convert RGB to intensity map (gray scale image)
    
%     a = iFrames(:,:,m);
%     figure(1); imagesc(a); colormap(gray);
    
    bFrame = bpass(iFrames(:,:,m),1,dia_guess);  % Filtering the raw frame
    %     figure(2);imagesc(b); colormap(gray);
    
    mx = max(max(bFrame));
    threshold = 0.8*mx;
    
    % Find the pixel of center of object
    pk = pkfnd(bFrame,threshold,dia_guess+1);  % The second term should be slightly larger than dia_guess for noisy data
    
    % Find the center of the object at subpixel resolution
    cnt = cntrd(bFrame,pk,dia_guess +2);
    xpixel(m) = cnt(1);  
    ypixel(m) = cnt(2);
    
    
    % Displaing the tracking 
    % =================================
%     imagesc(bFrame)
    subplot(121),imagesc(iFrames(:,:,m));
    title('Raw video (iFrames)');
    colormap gray;
    rectangle('Position',[xpixel(m)-dia_guess/2 ypixel(m)-dia_guess/2 dia_guess dia_guess],'edgecolor','r')
    xlim([1 frameSizeW]);
    ylim([1 frameSizeH]);
    xlabel('x (pixels)'); ylabel('y (pixels)');
    axis square;
    
    subplot(122), imagesc(bFrame);
    title('Filtered video (bFrames)');
    rectangle('Position',[xpixel(m)-dia_guess/2 ypixel(m)-dia_guess/2 dia_guess dia_guess],'edgecolor','r')
    xlim([1 frameSizeW]);
    ylim([1 frameSizeH]);
    xlabel('x (pixels)'); ylabel('y (pixels)');
    axis square;
    
    iFrames2(:,:,m) = getframe(gcf); 
    clf;   
    
    % This clf significantly reduces processing time. Othewise, the frames
    % were getting superimposed on top of each other and getframe command
    % was running slow.
    
    % ====================================
%     
%     RGB = insertShape(vidFrames(:,:,:,m),'circle',[xpixel(m) ypixel(m) dia_guess/2]);
%     iFrames2(:,:,m) = rgb2gray(RGB);
end
toc




pos_lst = [xpixel' ypixel' (1:length(xpixel))'];
tr = track(pos_lst, 3);
x_tr_cnt = 0.5*(max(tr(:,1)) + min(tr(:,1)));
y_tr_cnt = 0.5*(max(tr(:,2)) + min(tr(:,2)));


% Drift correction
if Ctrl_drift_adjust == 'y'
    [xcp, ycp] = drift_adj(xpixel,ypixel);
else
    xcp = xpixel; ycp = ypixel;
end

% Centering
% The mean position of the particle is set at (0,0) by subtracting the mean
% of the raw data.
if Ctrl_centering == 'y'
    xcp = xcp - mean(xcp);
    ycp = ycp - mean(ycp);
end
   

% Plotting raw tracked data
figure,
subplot(211),plot(nf,xpixel);
title('Raw position data');
xlabel('Frame no.'); ylabel('x (pixels)');
subplot(212),plot(nf,ypixel);
xlabel('Frame no.'); ylabel('y (pixels)');


 % Plotting drift corrected data 
if Ctrl_drift_adjust == 'y'
    figure,
    subplot(211),plot(nf,xcp,'r');
    title('Drift corrected data');
    xlabel('Frame no.'); ylabel('x (pixels)');
    subplot(212),plot(nf,ycp,'r');
    xlabel('Frame no.'); ylabel('y (pixels)');
end

% Plotting x, y position data
figure,
subplot(121),plot(xpixel,ypixel,'bo');
title('Raw position data');
xlabel('x (pixels)'); ylabel('y (pixels)');

if Ctrl_drift_adjust == 'y'
    subplot(122),plot(xcp,ycp,'ro');
    title('Drift corrected data');
    xlabel('x (pixels)'); ylabel('y (pixels)');
end


% Plotting histograms
figure,
subplot(121), hist(xcp,Nbins);
title('PDF of x position');
xlabel('x (pixels)'); ylabel('Count');
subplot(122), hist(ycp,Nbins);
title('PDF of y position');
xlabel('y (pixels)'); ylabel('Count');


% Video playback
if Ctrl_video_playback == 'y'
    implay(iFrames2)
end






% x_tr_adj = (tr(:,1) - x_tr_cnt)*dx;
% y_tr_adj = (tr(:,2) - y_tr_cnt)*dx;
% 
% 
% x_cnt = 0.5*(max(xpixel) + min(xpixel));
% y_cnt = 0.5*(max(ypixel) + min(ypixel));
% 
% x_adj = (xpixel - x_cnt)*dx;
% y_adj = (ypixel - y_cnt)*dx;
% 
% subplot(211), plot(x_adj,'b');
% hold on; plot(x_tr_adj,'r.');
% ylabel('x (nm)');
% xlabel('Frame number');
% 
% subplot(212), plot(y_adj,'b');
% hold on; plot(y_tr_adj,'r.');
% ylabel('y (nm)');
% xlabel('Frame number');
% 
% figure, plot(x_adj,y_adj,'bo-');
% % hold on, plot(x_tr_adj,y_tr_adj,'ro');
% 
% figure, subplot(121),imagesc(iFrames(:,:,end)); colormap(gray); axis square;
% 
% subplot(122),imagesc(bFrame); colormap(gray); axis square;

%%% 

% figure,
% imagesc(iFrames(:,:,1))
% colormap gray;
% hold on;
% rectangle('Position',[xpixel(1)-dia_guess/2 ypixel(1)-dia_guess/2 dia_guess dia_guess],'edgecolor','r')
% 





