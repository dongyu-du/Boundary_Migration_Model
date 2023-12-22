% Reconstruction procedures for "A boundary migration model for imaging within volumetric scattering media"
% by Dongyu Du, Xin Jin, and Rujia Deng et.al.
% Questions can be addressed to Dongyu Du(dudy19@mails.tsinghua.edu.cn)
%
% Copyright (c) 2023 Dongyu Du, Tsinghua University
% 
%    Licensed under the Apache License, Version 2.0 (the "License");
%    you may not use this file except in compliance with the License.
%    You may obtain a copy of the License at
% 
%        http://www.apache.org/licenses/LICENSE-2.0
% 
%    Unless required by applicable law or agreed to in writing, software
%    distributed under the License is distributed on an "AS IS" BASIS,
%    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%    See the License for the specific language governing permissions and
%    limitations under the License.

clear;
clc;

scene=8;
    
    switch scene             
        case {1}
             dataname = '8cmFoam_LambertianT_7cmFoam.mat';  
             load(['./data/' dataname]);
             width = 0.11;
             pad_length=148;rect1=14;rect2=79;   gamma =2.6;
        case {2}
             dataname = '7cmFoam_LambertianT_7cmFoam.mat';  
             load(['./data/' dataname]);
             width = 0.1;
             pad_length=126;rect1=15;rect2=67;   gamma =1.7;
        case {3}
             dataname = '6cmFoam_LambertianT_7cmFoam.mat';  
             load(['./data/' dataname]);
             width = 0.1;
             pad_length=120;rect1=18;rect2=68;   gamma =1.7;
        case {4}
             dataname = '5cmFoam_LambertianT_7cmFoam.mat';  
             load(['./data/' dataname]);
             width = 0.05;
             pad_length=119;rect1=18;rect2=64;   gamma =2.2;
        case {5}
             dataname = '4cmFoam_LambertianT_7cmFoam.mat';  
             load(['./data/' dataname]);
             width = 0.1;
             pad_length=112;rect1=22;rect2=65;   gamma =2.5;  
        case {6}
             dataname = '3cmFoam_LambertianT_7cmFoam.mat';  
             load(['./data/' dataname]);
             width = 0.1;
             pad_length=117;rect1=16;rect2=65;   gamma =1;              
        case {7}
             dataname = '2cmFoam_LambertianT_7cmFoam.mat';  
             load(['./data/' dataname]);
             width = 0.1;
             pad_length=65;rect1=18;rect2=42;   gamma =1.5;    
        case {8}
             dataname = 'Fog_Bunny.mat';  
             load(['./data/' dataname]);
             width = 0.1;
             pad_length=140;rect1=11;rect2=68;   gamma =1.7;  
                                            
    end        


    % Data Rectification
    raw_data= permute(raw_data(end:-1:1,end:-1:1,:),[2,1,3]);           
    rect_data = zeros(size(raw_data,1),size(raw_data,2),pad_length);
    rect_data(:,:,1:rect2-rect1+1)=raw_data(:,:,rect1:rect2);
    rect_data = rect_data(:,:,:);   
    

    % Constants Configuration
    u_a=0;                    % Ignore the absorption of the scattering media
    g=0;                      % Phase function parameter
    u_s=313;                  % Scattering coefficient
    n=1.23;                   % Relative refractive index
    D=1/(3*(u_a+(1-g)*u_s));  % Diffusion coefficient
    c = 3e8;                  % Speed of light in free-space
    ce=c/2/n;                 % Speed of light in the scattering media
    temporal_resolution =55e-12;   % Temporal resolution of the system        
    N = size(rect_data,1);            % Spatial resolution of data
    M = size(rect_data,3);            % Temporal resolution of data
    range = M.*c.*temporal_resolution;     % Maximum range

  
    % Data pad
    data = permute(rect_data,[3 2 1]);
    grid_z = repmat(linspace(0,1,M)',[1 N N]);
    data = data .* grid_z.^2;
    data = sqrt(data);
    
 
    % Phi(x,y,z=0,t) = measurement
    Phi_t=zeros(2*M,N,N);
    Phi_t(1:end/2,:,:)=data;
    % Construct the matrix A in Eq.8
    t=(0:2*M-1)';
    f_1=((0:2*M-1));
    A = exp(-t*(4*pi*pi*D.*(f_1.^2)+ce^2*u_a)./ce)./(2*pi); 
    % 1D numerical transform over t
    Phi_f = pagemtimes(pinv(A),Phi_t);
    % 2D Fourier transform over (x,y) to get Phi_bar(kx,ky,f)
    Phi_bar = zeros(2.*M,2.*N,2.*N);
    Phi_bar(:,1:end/2,1:end/2) = Phi_f;
    Phi_bar = fft(Phi_bar,'',2);
    Phi_bar = fftshift(fft(Phi_bar,'',3));      
    % Frequency interpolation: Phi_bar(kx,ky,f) -> Phi_prime(kx,ky,kz)
    [f,ky,kx] = ndgrid(-M:M-1,-N:N-1,-N:N-1);
    f = f./M; ky = ky./N; kx = kx./N;
    Phi_prime = interpn(f,ky,kx,Phi_bar,sqrt(abs((((N.*range)./(M.*width.*4)).^2).*(kx.^2+ky.^2)+f.^2)),ky,kx,'linear',0);
    Phi_prime = Phi_prime.*(f > 0);
    Phi_prime = Phi_prime.*abs(f)./max(sqrt(abs((((N.*range)./(M.*width.*4)).^2).*(kx.^2+ky.^2)+f.^2)),1e-6);
    % Phi_prime(kx,ky,kz) -> Phi(x,y,z,t=0) = object
    Object = ifftn(ifftshift(Phi_prime));


    Object = abs(Object).^2;    
    result = abs(Object(1:end./2,1:end./2,1:end./2));  
    result = result(:,:,end:-1:1);
    x = linspace(width,-width,size(result,3));
    y = linspace(width,-width,size(result,2));


    % Show the results 
    figure
    imagesc(x,y,squeeze(max(result(:,:,:),[],1)).^gamma);
    title('Front view of result');
    set(gca,'XTick',linspace(min(x),max(x),3));
    set(gca,'YTick',linspace(min(y),max(y),3));
    xlabel('x (m)');
    ylabel('y (m)');
    colormap('hot');
    axis square;