%% Create the concentration profile of an open MMI-Element
% Exemplary BPM calculation for a rectangular step-index MMI-element

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
ABSORBER = 0;                        % Small value needed to prevent diverging field, 0 to deactivate.
    
dx_fine   = .5e-6;                   % Step size for fine step
dy_fine   = 1e-6;                    % Step size for fine step
dx_coarse = 1e-6;                    % Step coarse for fine step               
dy_coarse = 1e-6;                    % Step coarse for fine step
dz        = 2e-6;                    % Propagationstep

Excitation.type = 'modal';           % Type of Excitation
Excitation.visualize = 1;            % Plot Excitation field
Excitation.sigma_x = 5e-6;           % Only needed for fieldtype 'gauss'
Excitation.sigma_y = 5e-6;           % Only needed for fieldtype 'gauss'
Excitation.threshold = 1-1/exp(1);   % Only needed for fieldtype 'full'
Excitation.numberInterpolation = 0;  % Only needed for fieldtype 'modal'
Excitation.mode = 'beta_z';          % Only needed for fieldtype 'modal'

% Some possible manipulations of index profile
append = 1;
append_dist = 4200e-6;  % um
prepend = 1;
prepend_dist = 10e-6;   % um

%% Generation of the refractive index profile

out = 'Generating index profile...';
disp(out)
tic
    
% Definition of fibers
n1 = 1.9;                  % Index of core 
n2 = 1.2;                  % Index of cladding
wInput = 5e-6;             % Width of rectangular input waveguide
wMMI = 20e-6;              % Width of rectangular input waveguide
beta_0 = 2*pi/lambda;      % Wave number

% Grid
xfl = 35e-6; % Absolute limits for fine mesh
yfl = 35e-6;
x = [-50e-6:dx_coarse:-xfl -xfl+dx_fine:dx_fine:xfl-dx_fine xfl:dx_coarse:50e-6];   % Base vector
y = [-50e-6:dy_coarse:-yfl -yfl+dy_fine:dy_fine:yfl-dy_fine yfl:dy_coarse:50e-6];   % Base vector
z = [0:1e-6:1e-6];                                                                  % Base vector
[xg,yg,zg] = meshgrid(x,y,z);                                                       % Base grid

xgb = squeeze(xg(:,:,1));  % Transversal grid of first slice
ygb = squeeze(yg(:,:,1));  % Transversal grid of first slice

% Refractive index profile input waveguide
nInput = n2*ones(length(y),length(x));
nInput(abs(xgb) <= wInput) = n1;
nInput(abs(ygb) >= wInput) = n2;
nMMI = n2*ones(length(y),length(x));
nMMI(abs(xgb) <= wMMI) = n1;
nMMI(abs(ygb) >= wMMI) = n2;
nc = cat(3,nInput,nMMI);

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

out = [POLARIZATION '-BPM with ' BC ' boundary condition started for step-index multimode interference based waveguide component at 1330nm wavelength.'];
disp(out);

[Ex,~,~,~] = FDBPMPade11Semivec(nc,lambda,n1,alpha,solver_tolerance,xg,yg,dz,Excitation,'TE',{'Ex';'Hy'},BC,ABSORBER,PROGRESS);

%% Visualization

Ex_int = squeeze(sum(abs(Ex) .* abs(yg)));
figure
surf(z(2:end),x,squeeze(Ex_int(:,2:end))/max(max(Ex_int(:,2:end))))
shading flat
xlabel('z [\mum]')
ylabel('x [\mum]')
title('Normalized absolute E-field (after y-integration) [a.u.]')
view([0 90])
axis([z(2) z(end) x(1) x(end) 0 1])