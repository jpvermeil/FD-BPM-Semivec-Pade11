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
ABSORBER = 0;                     % Small value needed to prevent diverging field, 0 to deactivate.
    
dx_fine   = .5e-6;                   % Step size for fine step
dy_fine   = .5e-6;                   % Step size for fine step
dx_coarse = 1e-6;                    % Step coarse for fine step               
dy_coarse = 1e-6;                    % Step coarse for fine step
dz        = 1e-6;                    % Propagationstep

Excitation.type = 'modal';           % Type of Excitation
Excitation.visualize = 1;            % Plot Excitation field
Excitation.sigma_x = 5e-6;           % Only needed for fieldtype 'gauss'
Excitation.sigma_y = 3e-6;           % Only needed for fieldtype 'gauss'
Excitation.threshold = 1-1/exp(1);   % Only needed for fieldtype 'full'
Excitation.numberInterpolation = 0;  % Only needed for fieldtype 'modal'
Excitation.mode = 'k_bar';           % Only needed for fieldtype 'modal'

% Some possible manipulations of index profile
append = 0;
append_dist = 0e-6;     % um
prepend = 0;
prepend_dist = 0e-6;    % um
upper_cladding = 1;

%% Generation of the refractive index profile

out = 'Loading index profile...';
disp(out)
tic
    
% Loading index profile and base vectors
load('Index_GI_MM_Diffused_Waveguide.mat')
n_max = max(max(n));
n_min = min(min(n));
nc = repmat(n,[1 1 101]);
beta_0 = 2*pi/lambda;       % Wave number
k_bar  = (n_max - n_min)/2; % Mean wave number

% Grid
z = [0:1e-6:100e-6];                                                                % Base vector
[xg,yg,zg] = meshgrid(x,y,z);                                                       % Base grid

xgb = squeeze(xg(:,:,1));  % Transversal grid of first slice
ygb = squeeze(yg(:,:,1));  % Transversal grid of first slice

out = ['  Index profile loaded: ' num2str(toc) 's.'];
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

%% Anbringung eines Claddingmaterials oberhalb des Substratblocks
    
if upper_cladding == 1
    
    out = 'Attaching upper cladding...';
    disp(out)
    tic

    % Generating new vectors and grids
    yCladd = y(end)+dy_fine:dy_fine:y(end)+10e-6;
    y      = [y yCladd];
    [xg,yg,zg] = meshgrid(x,y,z);

    % Removing round-off errors
    xg = 1e-12*round(xg*1e12);
    yg = 1e-12*round(yg*1e12);
    zg = 1e-12*round(zg*1e12);
     
    % Attaching upper cladding to existing index profile
    nCladd = 1.522*ones(length(yCladd),length(x),length(z));

    nc = cat(1,nc,nCladd);

    out = ['  Upper Cladding attached: ' num2str(toc) 's.'];
    disp(out)

end
 
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

out = [POLARIZATION '-BPM with ' BC ' boundary condition started for graded-index multimode diffused waveguide at 1330nm wavelength.'];
disp(out);

[Ex,~,~,~] = FDBPMPade11Semivec(nc,lambda,k_bar,alpha,solver_tolerance,xg,yg,dz,Excitation,'TE',{'Ex';'Hy'},BC,ABSORBER,PROGRESS);

%% Visualization

Ex_int = squeeze(sum(abs(Ex) .* abs(yg)));
figure
surf(z(2:end),x,squeeze(Ex_int(:,2:end))/max(max(Ex_int(:,2:end))))
shading flat
xlabel('z [\mum]')
ylabel('x [\mum]')
title('Normalized absolute E-field (after y-integration) [a.u.]')
axis([z(2) z(end) x(1) x(end) 0 1])