
%% Image and video processing: 
%  From Mars to Hollywood with a stop at the hospital
%
%  Instructions
%  ------------
%
% Use the level-sets method to implement constant motion. Consider simply 
% an image as the embedding function and deform it according to the 
% corresponding equation, I_t=|\divI|. Implement also I_t =?|\div I|. 
% Observe the result of both cases for different evolution intervals.
%
%@autor: Pedro Marrero
%@date: 17/02/2016
%
 

%% Iinitialize
clear; close all; clc;

%% Load image
I = imread('https://upload.wikimedia.org/wikipedia/en/2/24/Lenna.png');
I = max(double(I),0);
I = (I./max(I(:)));

if size(I,3) == 3 
I = rgb2gray(I);
end
I = double(I);
% I = imresize(I, 1/2);

In = imnoise(I,'salt & pepper',0.2);
In = In.*255;

figure(1); imshow(I,[]); title('Original image');
figure(2); imshow(In,[]); title('Noise image');


%%
[n,m] = size(In);
phi = double(In);
sig = +1; dt = 0.05; N = 200;

for t=1:N
          
    
    % boundary condition
    pphi = padarray(phi,[1 1],'symmetric');
    
    
    % updated according to the partial differential equation
    %   d phi        |          |
    %  ------- = +/- | grad(phi)|       (1)
    %    dt          |          |

    %Central different aproximationds
    % grad(phi) = dphi/dx *ex + dphi/dy * ey
    % grad(phi) = <gx,gy>
        
    i = 2:(n+1); j = 2:(m+1); ds = 1.0;
    gx  = (pphi(i,j+1) - pphi(i,j-1))./(2*ds);
    gy  = (pphi(i+1,j) - pphi(i-1,j))./(2*ds);  
    gxx = (pphi(i,j+1) + pphi(i,j-1) - 2*pphi(i,j))./(ds^2);
    gyy = (pphi(i+1,j) + pphi(i-1,j) - 2*pphi(i,j))./(ds^2);
       
    F = gxx + gyy;
	%F = del2(phi);
    
    
    % update
    phi = phi + (sig).*F.*dt;    
    
    figure(3); imshow(pphi,[]); title(['Smoothing t:' num2str(t*dt)]); drawnow;
    fprintf('Update t=%d\n', t*dt);
    
end

