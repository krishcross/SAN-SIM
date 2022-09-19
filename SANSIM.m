% Program for reconstruction using 'SAN SIM'
% Abbreviation of 'SAN SIM' = Saturable Absorption Nonlinear Structured Illumination Microscopy
% Copyright (c) Department of Physics, New Delhi, India
% Indian Institute of Technology Delhi (IITD)
% Year - 2022
%
% The algorithm was proposed in the article https://doi.org/10.1364/OL.460502
% Webpage - https://github.com/krishcross/SAN-SIM/
%
% Please cite the following paper if this program is used:
%
% K. Samanta, A. Tiwari, S. Joseph, and J. Joseph, "Saturable absorption assisted 
% nonlinear structured illumination microscopy", Optics Letters 47(11), 2702-2705 (2012). 
% --------------------------------------------------------------------
% Default parameters are set in the program. 
% You may like to change them as needed or keep them as they are.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forward Model %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;clc;
% close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Mesh
sps=40e-9;
[X,Y] = meshgrid(sps*(-256:255),sps*(-256:255));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PSF
NA=0.65;
lambda=560e-9;k=2*pi/lambda;
scale = NA*k; 
R=sqrt(X.^2+Y.^2);
PSF=abs(2*besselj(1,scale*R+eps,1)./(scale*R+eps)).^2; 
% figure,imagesc(PSF);axis image;colormap hot;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OTF
OTF2d=fft2(PSF);
OTF2dmax = max(max(abs(OTF2d)));
OTF2d = OTF2d./OTF2dmax;
OTF = abs(fftshift(OTF2d));
% figure,imagesc(OTF);axis image;colormap gray;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Synthetic Sample
[X,Y] = meshgrid(sps*(-256:255),sps*(-256:255));
[theta,r]=cart2pol(X,Y);
sample=double(0.5*(1+cos(120*theta))); % Siemens Star
%figure;imagesc(sample);axis image;c=0.2*morgenstemning(4000)+0.8*hot(4000);colormap(c);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Illumination Patterns
[xx,yy]=meshgrid(sps*(-256:255));
k=0.8*k;N=12;
amp = 30e8;
gR1 = amp*0.5*(1+sin(k*(cosd(0)*xx+sind(0)*yy)+0*2*pi/N));
gR2 = amp*0.5*(1+sin(k*(cosd(60)*xx+sind(60)*yy)+0*2*pi/N));
gR3 = amp*0.5*(1+sin(k*(cosd(120)*xx+sind(120)*yy)+0*2*pi/N));
for i=1:(N-1)
    gR1(:,:,i+1)= amp*0.5*(1+sin(k*(cosd(0)*xx+sind(0)*yy)+i*2*pi/N));
    gR2(:,:,i+1)= amp*0.5*(1+sin(k*(cosd(60)*xx+sind(60)*yy)+i*2*pi/N));
    gR3(:,:,i+1)= amp*0.5*(1+sin(k*(cosd(120)*xx+sind(120)*yy)+i*2*pi/N));
end


alpha = 3.7e4; 
beta = 4;
Is = 20e8;
factor2 = 0.002;
for i=1:N
factor1 = alpha./(1+gR1(:,:,i)./Is);
gRm1(:,:,i) = exp(-factor1.*factor2).*gR1(:,:,i);
factor1 = alpha./(1+gR2(:,:,i)./Is);
gRm2(:,:,i) = exp(-factor1.*factor2).*gR2(:,:,i);
factor1 = alpha./(1+gR3(:,:,i)./Is);
gRm3(:,:,i) = exp(-factor1.*factor2).*gR3(:,:,i);
end
gRm1 = gRm1./max(gRm1(:));gRm2 = gRm2./max(gRm2(:));gRm3 = gRm3./max(gRm3(:));
% figure,imagesc(gRm1(:,:,1));axis image;colormap gray;
% FT =  abs(fftshift(fft2(ifftshift(gRm1(:,:,1)))));
% FT = FT./max(FT(:));
% figure,plot(FT(257,:));axis tight;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Synthetic Raw Data 
for i = 1:N
    s1(:,:,i) = sample.*gRm1(:,:,i);
    s2(:,:,i) = sample.*gRm2(:,:,i);
    s3(:,:,i) = sample.*gRm3(:,:,i);
end

for i = 1:N
    Im1(:,:,i) = abs(ifftshift(ifft2(fftshift(fftshift(fft2(ifftshift(s1(:,:,i)))).*OTF))));
    Im2(:,:,i) = abs(ifftshift(ifft2(fftshift(fftshift(fft2(ifftshift(s2(:,:,i)))).*OTF))));
    Im3(:,:,i) = abs(ifftshift(ifft2(fftshift(fftshift(fft2(ifftshift(s3(:,:,i)))).*OTF))));
end
%figure();imagesc(Im1(:,:,1));axis square;axis off;c=0.7*morgenstemning(4000)+0.3*hot(4000);colormap(c);
%figure();imagesc(Im2(:,:,1));axis square;axis off;c=0.7*morgenstemning(4000)+0.3*hot(4000);colormap(c);
%figure();imagesc(Im3(:,:,1));axis square;axis off;c=0.7*morgenstemning(4000)+0.3*hot(4000);colormap(c);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inverse Model %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PSFd = abs(fftshift( ifft2(fftshift(OTF.^3)) ));
PSFd = PSFd/max(max(PSFd));
PSFd = PSFd/sum(sum(PSFd));
h = 10;
PSFe = PSFd(256-h+1:256+h,256-h+1:256+h);
for i = 1:N
    I1(:,:,i) = edgetaper(Im1(:,:,i),PSFe);
    I2(:,:,i) = edgetaper(Im2(:,:,i),PSFe);
    I3(:,:,i) = edgetaper(Im3(:,:,i),PSFe);
end

I(:,:,:,1)=I1;I(:,:,:,2)=I2;I(:,:,:,3)=I3;

DI0T = abs(ifftshift(ifft2(fftshift( fftshift(fft2(ifftshift(sample))).*OTF))) ); 
DI0T = edgetaper(DI0T,PSFe);
figure,imagesc(DI0T);axis image;c=0.7*morgenstemning(4000)+0.3*hot(4000);colormap(c);

MSI = zeros(512, 512, N, 3);
SPIO = zeros(512, 512, 3); % shift phase image per orientation
WFIO= mean(I, 3); % wide field image per orientation
wideFieldImage = mean(mean(I,3),4); % wide field image 

for iOrientation = 1:3
    for iPhase = 1:N
        MSI(:,:,iPhase,iOrientation) = (I(:,:,iPhase,iOrientation) -  WFIO(:,:,1,iOrientation)).^beta;
    end
    SPIO(:,:,iOrientation) = (mean(MSI(:,:,:,iOrientation), 3)).^(1/beta);
end
  
NLSIM = mean(SPIO,3);
NLSIM = NLSIM./max(NLSIM(:));

figure();imagesc(NLSIM);axis square;axis off;c=0.7*morgenstemning(4000)+0.3*hot(4000);colormap(c);


% FTSAMPLE = abs(fftshift(fft2(ifftshift(sample))));
% FTWF = abs(fftshift(fft2(ifftshift(sample))).*OTF);
% FTNLSIM = abs(fftshift(fft2(ifftshift(NLSIM))));
% figure;imagesc((FTNLSIM).^0.3);axis square;axis off;
% figure();imagesc((FTWF).^0.1);axis square;axis off;
% figure();imagesc((FTNLSIM).^0.1);axis square;axis off;

% imwrite(NLSIM,'NLSIM.tiff','WriteMode','append');
% imwrite(DI0T,'WF.tiff','WriteMode','append');