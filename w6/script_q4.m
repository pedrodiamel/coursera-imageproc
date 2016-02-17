%% Image and video processing: 
%  From Mars to Hollywood with a stop at the hospital
%
%  Instructions
%  ------------
%
% Implement the basic equation of the active contours and test it on some 
% images with simple objects. Initialize the evolving curve in different 
% forms for more testing.
%
%@autor: Pedro Marrero
%@date: 17/02/2016
%

%% Initialization
clear ; close all; clc

%% Load image

I = imread('https://upload.wikimedia.org/wikipedia/commons/thumb/0/08/Flag_of_Switzerland_%28Pantone%29.svg/320px-Flag_of_Switzerland_%28Pantone%29.svg.png');
I = I./max(I(:));
I = 255 - I.*255;

if size(I,3) == 3 
I = rgb2gray(I);
end
I = double(I);
I = imresize(I, 1/3);

figure(1); imshow(I,[]); title('Image original');


%% Initialization phi0

[n,m] = size(I);

% Iinitialize LSF as binary step function
c0=2; sqsize = 5;
initLSF = c0*ones(n,m);
vs{1} = (1:sqsize) + fix(n/2); vs{2} = (1:sqsize) + fix(m/2); 
initLSF(vs{:})=-c0; 
phi0 = initLSF;

% Show phi0
figure(2);
mesh(-phi0);   % for a better view, the LSF is displayed upside down
hold on;  contour(phi0, [0,0], 'r','LineWidth',2);
title('Initial level set function');
view([-80 35]);
hold off;

figure(3);
imagesc(I,[0, 255]); axis off; axis equal; colormap(gray); hold on;  contour(phi0, [0,0], 'r');
title('Initial zero level contour');


%% Image force calculate


% A driving expansion (or deflation) force which is synthesized from the
% image gradient.
sigma = 3.25;
Is = imgaussfilt(I,sigma, 'FilterSize',3);
[Ix,Iy] = gradient(Is);
gI = 1./(1 + Ix.^2 + Iy.^2);
gI = (gI./max(gI(:))).*255;


% A force which attracts the surface towards the boundary, which has a
% stabilizing effect, especially when there is a large variation in the 
% image gradient value.
[gIx,gIy] = gradient(gI);

% show g_I
figure(4); imshow(gI,[]);

%% Iteration methods 

sig = 1; dt = 0.0005; N = 500;

phi = double(phi0);
[n,m] = size(phi);

for t=1:N
      
    % boundary condition
    pphi = padarray(phi,[1 1],'symmetric'); 
    
    
    % Updated according to the partial differential equation (PDE)    

    %Central different aproximationds
    i = 2:(n+1); j = 2:(m+1); ds = 1.0;
    gx  = (pphi(i,j+1) - pphi(i,j-1))./(2*ds);
    gy  = (pphi(i+1,j) - pphi(i-1,j))./(2*ds);  
    gxx = (pphi(i,j+1) + pphi(i,j-1) - 2*pphi(i,j))./(ds^2);
    gyy = (pphi(i+1,j) + pphi(i-1,j) - 2*pphi(i,j))./(ds^2);
    gxy = (pphi(i+1,j+1) + pphi(i-1,j-1) - pphi(i-1,j+1) - pphi(i+1,j-1) )./( 4*ds^2 );
   
    % \div ( \grad \phi / |\grad \phi| )
    kappa = (gxx.*(gy.^2) + gyy.*(gx.^2) - 2*gxy.*gx.*gy)./((gx.^2 + gy.^2 + eps).^(1.5));
            
    % |\grad phi|
    gNorm = (gx.^2 + gy.^2).^(0.5);      
    
    % Flow evolution ecuation
    % \phi_t + g_I(1-\epsilon \kappa)|\grad \phi| + \beta \grad gI x \grad \phi
    
    epsilon = 0.25;
    beta = 0.9;
    F = -(gI.*(1-epsilon.*kappa).*gNorm) + beta.*(gIx.*gx + gIy.*gy); 
        
    % update
    phi = phi + (sig).*F.*dt;    
    
    figure(5); clf;
    imagesc(I,[0, 255]); axis off; axis equal; colormap(gray); hold on;
    contour(phi, [0 0], 'b','LineWidth',1);      
    drawnow;
    
    fprintf('Update t=%d\n', t*dt);
    
    
end







