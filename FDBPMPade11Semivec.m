function [phi,globalAdrSlgs,dim_xl,dim_yl] = FDBPMPade11Semivec(n,lambda,neff,alpha,solverTolerance,xg,yg,dz,EXCITATION,POLARIZATION,FIELDCOMPONENTS,BC,ABSORBER,PROGRESS)
% Semi vectorial wide angle Padé(1,1) finite difference BPM for TE/TM E-
% and/or H-fields in 3D structures.
% 
% SYNOPSIS
%
% FDBPMPade11Semivec(n,lambda,neff,alpha,solverTolerance,xg,yg,dz,EXCITATION,POLARIZATION,FIELDCOMPONENTS,BC,ABSORBER,PROGRESS)
% 
% VARIABLES
%
%   n               Refractive index profile
%   lambda          Optical wavelength
%   n_eff_Ex        Propagation constant (For wide angle Pade BPM this is usually (n1 + n2) / 2)
%   alpha           Discretization parameter in propagation direction
%   solver_tol      Tolerance of bicgstab solver
%   xg              Grid of x dimension. Has to match dimensions of n. 
%                   Has to be X output of Matlab meshgrid function ([X,Y,Z] = meshgrid(x,y,z))
%   yg              Grid of y dimension. Has to match dimensions of n.
%                   Has to be Y output of Matlab meshgrid function ([X,Y,Z] = meshgrid(x,y,z))
%   EXCITATION      Structural variable for the definition of exciting field distribution 
%       .fieldtype              Can be 'gauss', 'full', 'modal' or 'external'
%       .visualize_excitation   Numerical flat for visualization of excitation
%       (gaussian)
%           .sigma_x            Sigma of gaussian beam in x direction in [um]
%           .sigma_y            Sigma of gaussian beam in y direction in [um]
%       (full)
%           .threshold          Numerical threshold value for core illumination: 0 < threhold < 1
%                               E.g.: 0.9 -> 90% of maximum refractive index contrast illuminated
%                                     0.2 -> 20% of maximum refractive index contrast illuminated
%       (modal)
%           .nb_interpolation   Number of interpolations of n for determining propagation constants.
%                               This causes n with dimension j x k to be interpolated. Resulting 
%                               profile ni will have dimensions:
%                               ((j-1)^nb_interpolation + 1) x ((k-1)^nb_interpolation + 1)
%                               Note that this feature is not supported yet. This means, that no
%                               other value than '0' is supported. 
%           .mode               String value that can either be 'beta_z' or 'k_bar'.
%                               beta_z will use propagation constant of calculated fundamental mode
%                               k_bar will use those of wide angle Padé approximation
%       (external)  
%           .field              A field distribution can be provided at this point. Its size has to
%                               match the defined grid sizes.
%   POLARIZATION    String value that can either be 'TE' or 'TM'
%   FIELDCOMPONENTS Can be 'Ex', 'Ey', 'Hx' or 'Hy'
%   BC              String value speficying boundary condition. Currently only 'ABC' supported
%   ABSORBER        Numerical value for the application of an delta_n-dependent absorber. 
%                   The valid range for this value is: 0 <= absorber < 1.
%                   The value specifies the absorption range of delta_n. 
%                   While 0 no absorber is applied.
%                   Example: 0.1 -> Absorbing where n < 0.1 * delta_n (lowest 10%)
%                   Example: 0.9 -> Absorbing where n < 0.9 * delta_n (lowest 90%)
%   PROGRESS        Defines type of progress information.
%                   String value that can either be:
%                   'bar':  Visual waitbar
%                   'cl':   Command line (useful for batch mode)
%                   'off':  No progress indicator (useful for step mode)

format long

%% Definition of variables and constants

beta_0 = 2*pi/lambda;           % wave number
beta_z = beta_0 * neff;         % propagation constant as defined by neff
n_max = max(max(n(:,:,1)));     % maximum value of n
n_min = min(min(n(:,:,1)));     % minimum value of n
delta_n = n_max - n_min;        % maximum refractive index contrast

% Check field components input
[FX,~] = checkFieldComponents(FIELDCOMPONENTS); % Function checks which field component is to be evaluated
        
%% BPM specific parameters

dim_y   = size(n,1); % Number of global elements in y direction
dim_x   = size(n,2); % Number of global elements in x direction

dim_xl  = size(n,2) - 2; % Number of local elements in y direction (without the boundary of the structure)
dim_yl  = size(n,1) - 2; % Number of local elements in x direction (without the boundary of the structure)

phi = zeros(size(n,1),size(n,2),size(n,3)); % field matrix

grid.localIndex = zeros(size(n,1),size(n,2)); % Grid that speficies lokal element numbers
grid.localIndex(2:end-1,2:end-1) = reshape(1:1:dim_xl*dim_yl',dim_yl,dim_xl);

grid.globalIndex = zeros(size(n,1),size(n,2)); % Grid that speficies global element numbers. Note that in the global case (looking at the outer dimensions of the grid that include the boundary) there is no difference between 'Index' and 'Adress'. For this reason 'Index' is used in matrix form and 'Adr' is used as vector.
grid.globalIndex(1:end) = 1:1:length(grid.globalIndex(1:end));

grid.boundaryFlag = zeros(size(n,1),size(n,2));  % Grid that flags boundary elements of the structure
grid.boundaryFlag([1 end],[1:end 1:end]) = 1;
grid.boundaryFlag([1:end 1:end],[1 end]) = 1;

globalIndexSlgs   = grid.globalIndex(2:end-1,2:end-1); % Vector for global adresses of all lokal elements
globalAdrSlgs   = reshape(globalIndexSlgs,size(globalIndexSlgs,1)*size(globalIndexSlgs,2),1);
     
%% Definition of initial fields depending on chosen excitation

fieldtype = EXCITATION.type;
VisualizeExcitation = EXCITATION.visualize;

% Check for fieldtype 'gauss', 'full' or 'modal'

if strcmp(fieldtype,'gauss') % Excitation with a gaussian beam
    
    sigma_x = EXCITATION.sigma_x; % sigma of gaussian distribution as measure for the extention of the gaussian beam in x-direction in [um]
    sigma_y = EXCITATION.sigma_y; % sigma of gaussian distribution as measure for the extention of the gaussian beam in y-direction in [um]
    
    if isnumeric(sigma_x) && isnumeric(sigma_y)

        xg1 = squeeze(xg(:,:,1));
        yg1 = squeeze(yg(:,:,1));

        [r_max,c_max] = find(squeeze(n(:,:,1)) == max(max(n(:,:,1)))); % Those parameters assure, that the gaussian beam hits the profile at its maximum preventing mismatch

        xg1 = xg1 - xg1(r_max(1),c_max(1)); % x-grid for the definition of the gaussian beam
        yg1 = yg1 - yg1(r_max(1),c_max(1)); % y-grid for the definition of the gaussian beam

        phiInput = 1*exp(-xg1.^2/(2*(sigma_x))^2 -yg1.^2/(2*(sigma_y))^2);  % Definition of gaussian beam
        phi(:,:,1) = phiInput;                                              % Merging gaussian beam in global field matrix 
        
    else
       
        out = 'Invalid specification of excitation. EXCITATION.sigma_x and EXCITATION.sigma_y have to be numeric values for field type ''gauss''.';
        disp(out)
        return
        
    end
    
elseif strcmp(fieldtype,'full')        % Excitation with a constant illumination pattern
    
    threshold = EXCITATION.threshold;   % threshold value determines how much of the waveguide's core is excited by full illuminating source. It determines how big the part of delta_n is which is illuminated.
    % Expample: threhold = 0.9; -> Maximum 90% of delta_n is illuminated
    % Expample: threhold = 0.6; -> Maximum 60% of delta_n is illuminated
    
    if isnumeric(threshold) && (threshold < 1) && (threshold > 0) 
        
        nThreshold = min(min(n(:,:,1))) + delta_n * (1 - threshold);
        phiInput = zeros(size(n,1),size(n,2));
        phiInput(squeeze(n(:,:,1)) >= nThreshold) = 1;   % Illuminating the part of delta_n that lies above the given threshold

        phi(:,:,1) = phiInput; % Merging illumination pattern in global field matrix 
        
    else
        
        out = 'Invalid specification of threshold parameter. EXCITATION.threshold has to be numeric and 0 < threshold < 1.';
        disp(out)
        return
        
    end

elseif strcmp(fieldtype,'modal') % Excitation with the fundamental mode
    
    numberInterpolation = EXCITATION.numberInterpolation;   % Specify number of interpolations of refractive index profile for the determination of propagation constants 
    propagationConstant = EXCITATION.mode;                  % Mode parameter determines if beta_z of the fundamental mode or k_bar of Padé approximation should be used for propagation
    
    nInterpolated = interpn(squeeze(n(:,:,1)),numberInterpolation);   % This causes n with dimension j x k to be interpolated. Resulting profile ni will have dimensions ((j-1)^numberInterpolation + 1) x ((k-1)^numberInterpolation + 1)
    
    % Calculate propagation constants of fundamental mode
    [neff_pc,modenfeld] = FDPropagationconstantsSemivec(nInterpolated,beta_0,squeeze(xg(:,:,1)),squeeze(yg(:,:,1)),dim_y,dim_xl,dim_yl,grid.globalIndex,grid.localIndex,POLARIZATION,FX,1); 

    phiInput = zeros(size(nInterpolated,1),size(nInterpolated,2));
    phiInput(2:end-1,2:end-1) = reshape(modenfeld(:,1),dim_yl,dim_xl);  % Reshaping modefield vector to match size of global field matrix
    phiInput = phiInput/max(max(abs((phiInput))));            
    
    phi(:,:,1) = phiInput;  % Merging mode field in global field matrix 
    
    if strcmp(propagationConstant,'beta_z')     % Use beta_z of calculated modefield for propagation
        
        neff = neff_pc; 
        
    elseif ~strcmp(propagationConstant,'k_bar') % Use k_bar of Padé approximation for propagation
        
        out = 'Invalid specification of excitation. Possible choices for EXCITATION.mode for modal propagation are: ''beta_z'' and ''k_bar''.';
        disp(out)
        return
        
    end
    
elseif strcmp(fieldtype,'external') % Excitation with an externally defined field
    
    if isnumeric(EXCITATION.field) && size(EXCITATION.field,1) == size(n,1) && size(EXCITATION.field,2) == size(n,2)    % Dimension check
    
        phiInput = EXCITATION.field;
        phiInput = phiInput/max(max(abs((phiInput))));
        phi(:,:,1) = phiInput;  % Merging mode field in global field matrix
        
    else
        
        out = 'Invalid specification of external field distribution. EXCITATION.field has to be a matrix matching the dimesions of the outer grids.';
        disp(out)
        return
        
    end
    
else
    
    out = 'Invalid excitation field. Possible choices for EXCITATION.field are: ''gauss'', ''full'', ''modal'' or ''external''.';
    disp(out)
    return
    
end

if VisualizeExcitation == 1 % Visualize excitation
    
    figure
    surf(xg(:,:,1),yg(:,:,1),abs(phi(:,:,1)))
    shading flat
    
end

%% Defining gamma for handling different boundary conditions

if strcmp(BC,'ABC')
    
    gammaBoundaryCondition = 0;
    
elseif strcmp(BC,'TBC')
    
    out = 'Transparent boundary condition not yet supported. Please use absorbing boundary condition: ''ABC''.';
    disp(out)
    return
    
else
    
    out = 'Invalid specification of boundary condition. Possible choices are ''ABC'' or ''TBC''.';
    disp(out)
    return
    
end

%% Generate multistep paramters

% Usually those parameters do not change for the whole propagation
% distance. The are however explicitely depending on step size, alpha and
% propagation constant.

[UX,VX] = genMultistepVars11(dz,alpha,beta_z);  

tic
c = 1; % Global progress counter (waitbar)

if strcmp(PROGRESS,'bar') == 1 
    
    h = waitbar(0,'','Name',['Computing Padé ' POLARIZATION '-BPM with ' BC ' boundary condition...']);

end

%% Propagation in z direction

for kz = 1:1:size(n,3)-1
    
    % Extract known field values for current propagation step
    
    pb = phi(:,:,kz);
 
    %% Call functions for calculation of diagonals
    
    [ diagC,diagN,diagS,diagE,diagW ]      = diagonalsPade(beta_0,neff,n(:,:,kz+1),xg(:,:,kz+1),yg(:,:,kz+1),dim_y,dim_xl,dim_yl,grid,gammaBoundaryCondition,POLARIZATION,FX,BC);
    [ diagBC,diagBN,diagBS,diagBE,diagBW ] = diagonalsPade(beta_0,neff,n(:,:,kz),xg(:,:,kz),yg(:,:,kz),dim_y,dim_xl,dim_yl,grid,gammaBoundaryCondition,POLARIZATION,FX,BC);
    
    %% Apply multistep method
    
    for ii = 1:1:1  % Padé(1,1) only requires one multistep
        
        % Merge diagonals in system matrix
        
        A = sparse((size(n,1)-2)*(size(n,2)-2),(size(n,1)-2)*(size(n,2)-2));
        A = spdiags(1 + VX(ii)*diagC,0,A);
        A = spdiags([VX(ii)*diagN(2:end); 0],-1,A);
        A = spdiags([0; VX(ii)*diagS(1:end-1)],1,A);
        A = spdiags([zeros(dim_yl,1); VX(ii)*diagE(1:end-dim_yl)],dim_yl,A);
        A = spdiags([VX(ii)*diagW(dim_yl+1:end); zeros(dim_yl,1)],-dim_yl,A);
        
        % Compute right side of equation
        
        diagBC = 1 + UX(ii)*diagBC;
        diagBN = [UX(ii)*diagBN(2:end); 0];
        diagBS = [0; UX(ii)*diagBS(1:end-1)];
        diagBE = [zeros(dim_yl,1); UX(ii)*diagBE(1:end-dim_yl)];
        diagBW = [UX(ii)*diagBW(dim_yl+1:end); zeros(dim_yl,1)];
        
        % Multiply right side with known field values at 'kz'
        
        diagBC = pb(globalAdrSlgs)          .* diagBC;
        diagBN = pb(globalAdrSlgs - 1)      .* diagBN;
        diagBS = pb(globalAdrSlgs + 1)      .* diagBS;
        diagBE = pb(globalAdrSlgs + dim_y)  .* diagBE;
        diagBW = pb(globalAdrSlgs - dim_y)  .* diagBW;
        
        % Form vector for right side
        
        b = sparse(diagBC + diagBN + diagBS + diagBE + diagBW);

        %% Solve SoLE
        
        [phiSolution,~] = bicgstab(A,b,solverTolerance);
       
        % Reshape solution to match field matrix
        
        pb = zeros(dim_y,dim_x); 
        pb(2:end-1,2:end-1) = reshape(phiSolution,dim_yl,dim_xl);
        
        %% Apply absorber if specified
        
        if isnumeric(ABSORBER) && (ABSORBER > 0) && (ABSORBER < 1) && (ABSORBER ~= 0)
        
            nThreshold             = n_min + delta_n * ABSORBER;
            adr_n_threshold         = squeeze(n(:,:,kz)) <= nThreshold;
            pb(adr_n_threshold)    = 0;
            
        elseif ~isnumeric(ABSORBER) || (ABSORBER < 0) || (ABSORBER > 1)
            
            out = 'Invalid specification of absorber: 0 <= absorber < 1.';
            disp(out)
            return
            
        end

    end
    
    %% Merge solution in global field matrix
    
    phi(:,:,kz+1) = pb;
    
    %% Update progress bar
    
    if strcmp(PROGRESS,'bar') == 1
    
        stepPercent = 1;  % Defines step size for progress bar update

        if floor(100*kz/(size(n,3)-1)) > c*stepPercent 

          if ishandle(h)    % Check if waitbar has not been closed to prevent error

            minutesPassed      = floor(toc/60);
            secondsPassed      = toc - minutesPassed*60;
            minutesRemaining   = floor((toc/c)*((100/stepPercent)-c)/60);
            secondsRemaining   = (toc/c)*((100/stepPercent)-c) - minutesRemaining*60;

            waitbar(kz/(size(n,3)-1),h,['Progress: ' num2str(c*stepPercent) '%. Time remaining: ' num2str(floor(minutesRemaining)) 'm ' num2str(ceil(secondsRemaining)) 's.']) 

            c = c + 1;

          else 

             disp('Process interrupted') 
             break

          end 
        end 
   
    elseif strcmp(PROGRESS,'cl') == 1
        
        stepPercent = 10;  % Defines step size for progress bar update

        if floor(100*kz/(size(n,3)-1)) > c*stepPercent 

            minutesPassed      = floor(toc/60);
            secondsPassed      = toc - minutesPassed*60;
            minutesRemaining   = floor((toc/c)*((100/stepPercent)-c)/60);
            secondsRemaining   = (toc/c)*((100/stepPercent)-c) - minutesRemaining*60;

            out = ['   Progress: ' num2str(c*stepPercent) '%. Time remaining: ' num2str(floor(minutesRemaining)) 'm ' num2str(ceil(secondsRemaining)) 's.'];
            disp(out)

            c = c + 1;
            
        end    
        
    end
        
end

try
    close(h)
end

minutesPassed     = floor(toc/60);
secondsPassed    = toc - minutesPassed*60;

if strcmp(PROGRESS,'off') == 0

    out = ['BPM-calculation finished after: ' num2str(minutesPassed) 'm ' num2str(secondsPassed) 's.'];
    disp(out);
    
end

end
