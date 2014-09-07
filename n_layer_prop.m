close all; clear all; format short; clc;
I = sqrt(-1);
%% Parameters of incident ray
theta = pi/6; % incident ray angle
lambda = 1.0; % incident ray wavelength
c = 1.0; % speed of light
H0 = 1.0; % amplitude of incident ray
%% Parameters of medium
N_layers = 5; % number of layers
dl = 1.0; % width of layer
d = dl*[0:N_layers]; % distances to the begining of layers along x axis
eps = [1.0, 3.0, 7.0, 5.0, 10.0, 4.0, 1.0]; % permittivities of layers
mu = [1.0, 3.0, 2.0, 2.5, 4.0, 3.0, 1.0]; % permeabilities of layers
%% Propagation matrices
k = 2*pi/lambda; % wave vector of incident ray
kx = [k*cos(theta)]; % x component of wave vector k of incident ray
kz = k*sin(theta); % z component of wave vector k of incident ray
p = [0]; %
r = [0]; % reflection coefficients
for i = 1:(N_layers + 1)
    kx(i+1) = sqrt(k^2*c^2*eps(i+1)*mu(i+1) - kz^2); % extending array of kx in layers
    p(i) = kx(i)*eps(i+1)/(kx(i+1)*eps(i)); % auxiliary parameter
    r(i) = (1 - p(i))/(1 + p(i)); % reflection coefficients
    % propagation matrices
    V(:,:,i) = 0.5*(1 + p(i))*[exp(-I*(kx(i+1) - kx(i))*d(i)) r(i)*exp(-I*(kx(i+1) + kx(i))*d(i)); ...
                               r(i)*exp(I*(kx(i+1) + kx(i))*d(i)) exp(I*(kx(i+1) - kx(i))*d(i))];
end
%% Propagation
% Calculation of amplitudes A, B is done with propagation matrices:
% (A(l + 1), B(l + 1)) = V(l + 1)(A(l), B(l)), l - layer number.
% Here we suppose R = 0.01 to check the propagation, but it should be
% calculated exactly for consistency.
R = 0.01;
% set coordinates (x, z) for evaluation of field 
x = 0.01*[-50:550];
z = 0.0;
H = zeros(601);
for j=1:601
    i = 1;
    A = transpose([H0, R]);
    while (i-1)*dl < x(j)
        A = V(:,:,i)*A;
        i = i + 1;
    end
    H(j) = (A(1)*exp(I*kx(i)*x(j)) + A(2)*exp(-I*kx(i)*x(j)))*exp(I*kz*z);
end
plot(x, abs(H));