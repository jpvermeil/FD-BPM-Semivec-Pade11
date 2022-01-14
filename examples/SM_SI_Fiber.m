%% Create the concentration profile of an open MMI-Element

% Parameter definition
clear
close all
clc

PROGRESS = 'cl';

lambda              = 1330e-9;       % Wavelength
alpha               = 0.5;           % Crank-Nicolson-Constant
delta               = 1e-6;
alpha               = alpha + delta;
solver_tolerance    = 1e-12;         % Tolerance of BiCGM solver

POLARIZATION = 'TE';
FIELDCOMPONENT = 'Ex';
BC = 'ABC';
ABSORBER = 1e-3;                     % Small value needed to prevent diverging field. This can only be avoided by proper implementation of TBC, which is yet to be done 
    
dx_fine   = .1e-6;                   % Step size for fine step
dy_fine   = .1e-6;                   % Step size for fine step
dx_coarse = 2e-6;                    % Step coarse for fine step               
dy_coarse = 2e-6;                    % Step coarse for fine step
dz        = 1e-6;                    % Propagationstep

Excitation.type = 'modal';           % Type of Excitation
Excitation.visualize = 1;            % Plot Excitation field
Excitation.sigma_x = 5e-6;           % Only needed for fieldtype 'gauss'
Excitation.sigma_y = 5e-6;           % Only needed for fieldtype 'gauss'
Excitation.threshold = 1-1/exp(1);   % Only needed for fieldtype 'full'
Excitation.numberInterpolation = 0;  % Only needed for fieldtype 'modal'
Excitation.mode = 'beta_z';          % Only needed for fieldtype 'modal'

% Some possible manipulations of index profile
append = 1;
append_dist = 10e-6;  % um
prepend = 1;
prepend_dist = 10e-6; % um

%% Generation of the refractive index profile

out = 'Generating index profile...';
disp(out)
tic
    
% Definition of SI-Fiber
n1 = 1.45;                  % Index of core 
n2 = 1.448;                 % Index of cladding
r_core = 5e-6;              % Radius
beta_0 = 2*pi/lambda;       % Wave number
NA = sqrt(n1^2-n2^2);       % Numerical aperture
R = 2*pi*r_core*NA/lambda;  % Fiber parameter

% Grid
xfl = 15e-6; % Absolute limits for fine mesh
yfl = 15e-6;
x = [-50e-6:dx_coarse:-xfl -xfl+dx_fine:dx_fine:xfl-dx_fine xfl:dx_coarse:50e-6];   % Base vector
y = [-50e-6:dy_coarse:-yfl -yfl+dy_fine:dy_fine:yfl-dy_fine yfl:dy_coarse:50e-6];   % Base vector
z = [0:1e-6:100e-6];                                                % Base vector
[xg,yg,zg] = meshgrid(x,y,z);                                       % Base grid

xgb = squeeze(xg(:,:,1));  % Transversal grid of first slice
ygb = squeeze(yg(:,:,1));  % Transversal grid of first slice

% Refractive index profile
n = n2*ones(length(x),length(y));
n(sqrt(xgb.^2+ygb.^2) < r_core) = n1;
nc = repmat(n,[1 1 length(z)]);

out = ['  Index profile generated: ' num2str(toc) 's.'];
disp(out)
tic
    
%% Append exit
% Can be used to replicate the last slice

out = 'Prepending and/or appending first and last slice...';
disp(out)
tic

if append == 1

    append_steps = floor(append_dist/dz);

    nc_end  = repmat(squeeze(nc(:,:,end)),[1 1 append_steps]);
    nc      = cat(3,nc,nc_end); 

    z_end  = z(end)+dz:dz:z(end)+append_dist;
    z      = [z z_end]; 

    [xg,yg,zg] = meshgrid(x,y,z);

end

%% Prepend input
% Can be used to replicate the first slice

if prepend == 1

    prepend_steps = floor(prepend_dist/dz);

    nc_input = repmat(squeeze(nc(:,:,1)),[1 1 prepend_steps]);
    nc       = cat(3,nc_input,nc);

    z_input = 0:dz:prepend_dist-dz;
    z = [z_input z+prepend_dist];

    [xg,yg,zg] = meshgrid(x,y,z);

end

out = ['  Prepending/appending done: ' num2str(toc) 's.'];
disp(out)
tic

%% BPM flags and parameters

out = 'Defining variables and BPM flags...';
disp(out)
tic

dim_y   = size(nc,1); % Global dimension
dim_x   = size(nc,2); 

dim_yl   = dim_y - 2; % Local dimensions (without boundary values)
dim_xl   = dim_x - 2;

grid.localIndex = zeros(size(nc,1),size(nc,2));   % Grid variables
grid.globalIndex = zeros(size(nc,1),size(nc,2));

grid.boundaryFlag = zeros(size(nc,1),size(nc,2)); % Boundary flag
grid.boundaryFlag([1 end],[1:end 1:end]) = 1;
grid.boundaryFlag([1:end 1:end],[1 end]) = 1;

grid.localIndex(2:end-1,2:end-1) = reshape(1:1:dim_xl*dim_yl',dim_yl,dim_xl);
grid.globalIndex(1:end) = 1:1:length(grid.globalIndex(1:end));

out = ['  Defined necessary parameters: ' num2str(toc) 's.'];
disp(out)

%% Execute BPM

out = [POLARIZATION '-BPM with ' BC ' boundary condition started for SM-SI-Fibre with core radius of 5um at 1330nm wavelength.'];
disp(out);

[Ex,~,~,~] = FDBPMPade11Semivec(nc,lambda,n1,alpha,solver_tolerance,xg,yg,dz,Excitation,'TE',{'Ex';'Hy'},BC,ABSORBER,PROGRESS);

%% Visualization

Ex_int = squeeze(sum(abs(Ex) .* abs(yg)));
figure
surf(z(2:end),y,squeeze(Ex_int(:,2:end))/max(max(Ex_int(:,2:end))))
shading flat
xlabel('z [\mum]')
ylabel('x [\mum]')
title('Normalized absolute E-field (after y-integration) [a.u.]')