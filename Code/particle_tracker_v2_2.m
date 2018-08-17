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
% =========================================================================
Ctrl_drift_adjust = 'y';
Ctrl_centering = 'y';
Ctrl_video_playback = 'n';
Ctrl_disp_tracking = 'y';

Ctrl_use_dim = 'y';
Ctrl_vid_write = 'y';
Ctrl_data_write = 'y';
Ctrl_pdf_fit = 'y';
Ctrl_frame_fraction = 1;

Nbins = 12; % number of bins in the histogram

% =========================================================================

% Indicate the path of the video file

source_file = 'Video/test1.avi';
% source_file = 'Video/Images.avi';

% source_file = 'Video\Z#7C70N2c3_D500_P196b.avi';



% Read video file
videoObj = VideoReader(source_file);
vidFrames = read(videoObj);

% Identifying frame size
[frameSizeH, frameSizeW, numColors, numFrames] = size(vidFrames);
frame_rate = videoObj.FrameRate;
% Initializing output frames
iFrames = zeros(frameSizeH, frameSizeW, numFrames);


% Dimension parameters
dx = 42/1000;  % Size of one pixel
dt = 1/frame_rate;

dia_guess = 25; % Guess of the diameter in pixels.


nf = 1:floor(numFrames*Ctrl_frame_fraction);  % Frame indices

tic
for m = 1:length(nf)
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
    if Ctrl_disp_tracking == 'y';
        %     imagesc(bFrame)
        subplot(121),imagesc(iFrames(:,:,m));
        title(['Raw video (iFrames = ', num2str(m),')']);
        colormap gray;
        rectangle('Position',[xpixel(m)-dia_guess/2 ypixel(m)-dia_guess/2 dia_guess dia_guess],'edgecolor','r','linewidth',1.1)
        xlim([1 frameSizeW]);
        ylim([1 frameSizeH]);
        xlabel('x (pixels)'); ylabel('y (pixels)');
        %     axis equal;
        pbaspect([1 1 1])
        
        subplot(122), imagesc(bFrame);
        title(['Filtered video (bFrames = ', num2str(m),')']);
        rectangle('Position',[xpixel(m)-dia_guess/2 ypixel(m)-dia_guess/2 dia_guess dia_guess],'edgecolor','r','linewidth',1.1)
        xlim([1 frameSizeW]);
        ylim([1 frameSizeH]);
        xlabel('x (pixels)'); ylabel('y (pixels)');
        %     axis equal;
        pbaspect([1 1 1])
        
        iFrames2(:,:,m) = getframe(gcf);
        clf;
        
        % This clf significantly reduces processing time. Othewise, the frames
        % were getting superimposed on top of each other and getframe command
        % was running slow.
        
        % ====================================
    end
%     
%     RGB = insertShape(vidFrames(:,:,:,m),'circle',[xpixel(m) ypixel(m) dia_guess/2]);
%     iFrames2(:,:,m) = rgb2gray(RGB);
end
toc




% pos_lst = [xpixel' ypixel' (1:length(xpixel))'];
% tr = track(pos_lst, 3);
% x_tr_cnt = 0.5*(max(tr(:,1)) + min(tr(:,1)));
% y_tr_cnt = 0.5*(max(tr(:,2)) + min(tr(:,2)));


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
   

if Ctrl_use_dim == 'y'
    scl_1 = dx;
    scl_2 = dt;
    xstr = 'x (\mum)';
    ystr = 'y (\mum)';
    tstr = 't (sec)';
else
    scl_1 = 1;
    scl_2 = 1;
    xstr = 'x (pixels)';
    ystr = 'y (pixels)';
    tstr = 'n (Frame no.)';
end

xpixel = xpixel * scl_1;
ypixel = ypixel * scl_1;
xcp = xcp * scl_1;
ycp = ycp * scl_1;
nf = nf *scl_2;


% Calculating standard deviation
fprintf('Standard deviation of x position = %1.4f units\n',sqrt(var(xcp)));
fprintf('Standard deviation of y position = %1.4f units\n',sqrt(var(ycp)));

if Ctrl_use_dim == 'y'
    fprintf('Units in microns.\n');
else
    fprintf('Units in pixels.\n');
end


% Plotting raw tracked data
figure,
subplot(211),plot(nf,xpixel);
title('Raw position data');
xlabel(tstr); ylabel(xstr);
subplot(212),plot(nf,ypixel);
xlabel(tstr); ylabel(ystr);


 % Plotting drift corrected data 
if Ctrl_drift_adjust == 'y'
    figure,
    subplot(211),plot(nf,xcp,'r');
    title('Drift corrected data');
    xlabel(tstr); ylabel(xstr);
    subplot(212),plot(nf,ycp,'r');
    xlabel(tstr); ylabel(ystr);
end

% Plotting x, y position data
figure,
subplot(121),plot(xpixel,ypixel,'bo');
title('Raw position data');
xlabel(xstr); ylabel(ystr);

if Ctrl_drift_adjust == 'y'
    subplot(122),plot(xcp,ycp,'ro');
    title('Drift corrected data');
    xlabel(xstr); ylabel(ystr);
end


% Plotting histograms
figure,
subplot(211), hist(xcp,Nbins);
title('PDF of x position');
xlabel(xstr); ylabel('Count');
subplot(212), hist(ycp,Nbins);
title('PDF of y position');
xlabel(ystr); ylabel('Count');


% Video playback
if Ctrl_video_playback == 'y'
    implay(iFrames2)
end

% Write video file
if Ctrl_vid_write == 'y'
    v = VideoWriter('tracked_video.avi');
    v.FrameRate = frame_rate;
    v.Quality = 25;
    open(v);
    writeVideo(v,iFrames2);
    close(v);
end

% PDF fitting
if Ctrl_pdf_fit == 'y'
    pdf_fit(xcp,ycp,Nbins,xstr,ystr);
end

if Ctrl_data_write == 'y'
    Mdata = [xpixel' ypixel' xcp' ycp'];
    csvwrite('data_file.csv',Mdata);
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





