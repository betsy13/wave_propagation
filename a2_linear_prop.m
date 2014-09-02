%%%%%% a2_linear_prop
%%%%%% R. Matt Leone, Fall 2007
%%%%%% This program demonstrates linear propagation
%%%%%% of light along z via the paraxial equation.
close all; clear all; format short; clc;
%% Inputs, size in microns.
waist = 10;
lambda = 1;
r_obstacle = 3;
z_max = pi*waist^2/lambda / 4;
extend_factor = 4;
angular_momentum = 1; %% L for Laguerre term.
%% Functions to create input field.
f_gauss = inline('(exp(-(X.^2+Y.^2)/waist.^2))','X','Y','waist');
f_laguerre = inline('(X+sqrt(-1)*Y).^L','X','Y','L');
f_mask = inline('1-exp(-((X.^2+Y.^2)/r_obstacle.^2).^16)','X','Y','r_obstacle');
%% Tricky scaling section.
n_points = 200;
x_max = waist / 2 * extend_factor;
dx = 2 * x_max / n_points;
k_max = pi / dx;
dk = k_max / n_points * 2;
dz = z_max;
Vx = dx*[-n_points:n_points];
[X,Y] = meshgrid(Vx,Vx);
Vk = ifftshift(dk*[-n_points:n_points]); %% DFT shift here.
[Kx,Ky] = meshgrid(Vk,Vk);
k_wave_number = 2 * pi / lambda;
Kz = -(Kx.^2 + Ky.^2) / (2 * k_wave_number);
Freq_Shift = exp(sqrt(-1) * dz * Kz);
%% Get fields.
Field = zeros(size(X));
Field(:,:) = f_gauss(X,Y,waist);
%Field(:,:) = f_gauss(X,Y,waist).*f_laguerre(X,Y,angular_momentum).*f_mask(X,Y,r_obstacle);
Field = Field ./ max(max(Field));
%% Propagate
figure(1);
imagesc(abs(Field.^2),[0 1]);
set(gca,'FontSize',15);
axis off;
title('Initial Intensity (r space)');
Field = fft2(ifft2(Field).*Freq_Shift);
Field = Field ./ max(max(Field));
figure(2);
imagesc(abs(Field.^2),[0 1]);
set(gca,'FontSize',15);
axis off;
title('z=waist/5 (intensity scaled to show dispersion)');