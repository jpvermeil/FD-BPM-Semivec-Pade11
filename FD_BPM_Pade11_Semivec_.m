function [Px,glob_adr_slgs,dim_xl,dim_yl] = FD_BPM_Pade11_Semivec_(n,lambda,n_eff_Px,alpha,solver_tol,xg,yg,dz,EXCITATION,POLARIZATION,FIELDCOMPONENTS,BC,ABSORBER,PROGRESS)
% Semi vectorial wide angle Padé(1,1) finite difference BPM for TE/TM E-
% and/or H-fields in 3D structures.
% 
% SYNOPSIS
%
% FD_BPM_Pade11_Semivec_(n,lambda,n_eff_Px,alpha,solver_tol,xg,yg,dz,EXCITATION,POLARIZATION,FIELDCOMPONENTS,BC,ABSORBER,PROGRESS)
% 
% VARIABLES
%
%   n               Refractive index profile
%   lambda          Optical wavelength
%   alpha           Discretization parameter in propagation direction
%   solver_tol      Tolerance of bicgstab solver
%   xg              Grid of x dimension. Has to match dimensions of n. 
%                   Has to be X output of Matlab meshgrid function ([X,Y,Z] = meshgrid(x,y,z))
%   yg              Grid of y dimension. Has to match dimensions of n.
%                   Has to be Y output of Matlab meshgrid function ([X,Y,Z] = meshgrid(x,y,z))
%   EXCITATION      Structural variable for the definition of exciting field distribution 
%       .fieldtype              Can be 'gauss', 'full' or 'modal'
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
%                               other value than »0« is supported. 
%           .mode               String value that can either be 'beta_z' or 'k_bar'.
%                               beta_z will use propagation constant of calculated fundamental mode
%                               k_bar will use those of wide angle Padé approximation  
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

format long

%% Definition of variables and constants

beta_0 = 2*pi/lambda;
beta_z = beta_0 * n_eff_Px;
n_max = max(max(squeeze(n(:,:,1))));
n_min = min(min(squeeze(n(:,:,1))));
delta_n = n_max - n_min;                % refractive index contrast

% Check field components input
[FX,FY] = check_fieldcomponents_(FIELDCOMPONENTS);
        
%% BPM specific parameters

dim_y   = size(n,1); % Number of global elements in y direction
dim_x   = size(n,2); % Number of global elements in x direction

dim_xl  = size(n,2) - 2; % Number of local elements in y direction (without the boundary of the structure)
dim_yl  = size(n,1) - 2; % Number of local elements in x direction (without the boundary of the structure)

Px = zeros(size(n,1),size(n,2),size(n,3));

dG.cl = zeros(size(n,1),size(n,2)); % Grid that speficies lokal element numbers
dG.cl(2:end-1,2:end-1) = reshape(1:1:dim_xl*dim_yl',dim_yl,dim_xl);

dG.cg = zeros(size(n,1),size(n,2)); % Grid that speficies global element numbers
dG.cg(1:end) = 1:1:length(dG.cg(1:end));

dG.f = zeros(size(n,1),size(n,2));  % Grid that flags boundary elements of the structure
dG.f([1 end],[1:end 1:end]) = 1;
dG.f([1:end 1:end],[1 end]) = 1;

glob_adr_slgs   = dG.cg(2:end-1,2:end-1); % Vector for global adresses of all lokal elements
glob_adr_slgs   = reshape(glob_adr_slgs,size(glob_adr_slgs,1)*size(glob_adr_slgs,2),1);
     
%% Definition of initial fields depending on chosen excitation

excitation = EXCITATION.fieldtype;
visualize_excitation = EXCITATION.visualize;

% Check for fieldtype 'gauss', 'full' or 'modal'

if strcmp(excitation,'gauss') % Excitation with a gaussian beam
    
    sigma_x = EXCITATION.sigma_x; % sigma of gaussian distribution as measure for the extention of the gaussian beam in x-direction in [um]
    sigma_y = EXCITATION.sigma_y; % sigma of gaussian distribution as measure for the extention of the gaussian beam in y-direction in [um]
    
    if isnumeric(sigma_x) && isnumeric(sigma_y)

        xg1 = squeeze(xg(:,:,1));
        yg1 = squeeze(yg(:,:,1));

        [r_max,c_max] = find(squeeze(n(:,:,1)) == max(max(squeeze(n(:,:,1))))); % Those parameters assure, that the gaussian beam hits the profile at its maximum preventing mismatch

        xg1 = xg1 - xg1(r_max(1),c_max(1)); % x-grid for the definition of the gaussian beam
        yg1 = yg1 - yg1(r_max(1),c_max(1)); % y-grid for the definition of the gaussian beam

        Px_a = 1*exp(-xg1.^2/(2*(sigma_x))^2 -yg1.^2/(2*(sigma_y))^2);  % Definition of gaussian beam
        Px(:,:,1) = Px_a;                                               % Merging gaussian beam in global field matrix 
        
    else
       
        out = 'Invalid specification of excitation. EXCITATION.sigma_x and EXCITATION.sigma_y have to be numeric values for field type ''gauss''.';
        disp(out)
        return
        
    end
    
elseif strcmp(excitation,'full')        % Excitation with a constant illumination pattern
    
    threshold = EXCITATION.threshold;   % threshold value determines how much of the waveguide's core is excited by full illuminating source. It determines how big the part of delta_n is which is illuminated.
    % Expample: threhold = 0.9; -> Maximum 90% of delta_n are illuminated
    % Expample: threhold = 0.6; -> Maximum 60% of delta_n are illuminated
    
    if isnumeric(threshold) && (threshold < 1) && (threshold > 0) 
        
        n_threshold = min(min(squeeze(n(:,:,1)))) + delta_n * (1 - threshold);
        Px_a = zeros(size(n,1),size(n,2));
        Px_a(find(squeeze(n(:,:,1)) >= n_threshold)) = 1;   % Illuminating the part of delta_n that lies above the given threshold

        Px(:,:,1) = Px_a; % Merging illumination pattern in global field matrix 
        
    else
        
        out = 'Invalid specification of threshold parameter. EXCITATION.threshold has to be numeric and 0 < threshold < 1.';
        disp(out)
        return
        
    end

elseif strcmp(excitation,'modal') % Excitation with the fundamental mode
    
    nb_interpolation = EXCITATION.nb_interpolation;     % Specify number of interpolations of refractive index profile for the determination of propagation constants 
    propagation_constant = EXCITATION.mode;             % Mode parameter determines if beta_z of the fundamental mode or k_bar of Padé approximation should be used for propagation
    
    ni = interpn(squeeze(n(:,:,1)),nb_interpolation);   % This causes n with dimension j x k to be interpolated. Resulting profile ni will have dimensions ((j-1)^nb_interpolation + 1) x ((k-1)^nb_interpolation + 1)
    
    % Calculate propagation constants of fundamental mode

    [n_eff_pc_Px,Modenfeld_Px] = FD_propagationconstants_Semivec_(ni,beta_0,squeeze(xg(:,:,1)),squeeze(yg(:,:,1)),dim_y,dim_xl,dim_yl,dG.cg,dG.cl,POLARIZATION,FX,1); 

    Px_a = zeros(size(ni,1),size(ni,2));

    Px_a(2:end-1,2:end-1) = reshape(Modenfeld_Px(:,1),dim_yl,dim_xl);   % Reshaping modefield vector to match size of global field matrix

    Px_a = Px_a/max(max(abs((Px_a))));
    
    Px(:,:,1) = Px_a; % Merging mode field in global field matrix 
    
    if strcmp(propagation_constant,'beta_z') % Use beta_z of calculated modefield for propagation
        
        n_eff_Px = n_eff_pc_Px; 
        
    elseif ~strcmp(propagation_constant,'k_bar') % Use k_bar of Padé approximation for propagation
        
        out = 'Invalid specification of excitation. Possible choices for EXCITATION.mode for modal propagation are: ''beta_z'' and ''k_bar''.';
        disp(out)
        return
        
    end
    
else
    
    out = 'Invalid excitation field. Possible choices for EXCITATION.field are: ''gauss'', ''full'' or ''modal''.';
    disp(out)
    return
    
end

if visualize_excitation == 1 % Visualize excitation
    
    figure
    surf(xg(:,:,1),yg(:,:,1),Px(:,:,1))
    shading flat
    
end

%% Defining gamma for handling different boundary conditions

if strcmp(BC,'ABC')
    
    gamma_bc = 0;
    
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

[ux,vx] = gen_multistep_vars_11_(dz,alpha,beta_z);  

tic
c = 1; % Global progress counter (waitbar)

if strcmp(PROGRESS,'bar') == 1 
    
    h = waitbar(0,'','Name',['Computing Padé ' POLARIZATION '-BPM with ' BC ' boundary condition...']);

end

%% Propagation in z direction

for kz = 1:1:size(n,3)-1
    
    % Extract known field values for current propagation step
    
    Pbx = Px(:,:,kz);
 
    %% Call functions for calculation of diagonals
    
    [ Cxx,Nxx,Sxx,Exx,Wxx ]      = diagonals_pade_(beta_0,n_eff_Px,n(:,:,kz+1),xg(:,:,kz+1),yg(:,:,kz+1),dim_y,dim_xl,dim_yl,dG,gamma_bc,POLARIZATION,FX,BC);
    [ Cbxx,Nbxx,Sbxx,Ebxx,Wbxx ] = diagonals_pade_(beta_0,n_eff_Px,n(:,:,kz),xg(:,:,kz),yg(:,:,kz),dim_y,dim_xl,dim_yl,dG,gamma_bc,POLARIZATION,FX,BC);
    
    %% Apply multistep method
    
    for ii = 1:1:1
        
        % Paste diagonals in system matrix
        
        Axx = sparse((size(n,1)-2)*(size(n,2)-2),(size(n,1)-2)*(size(n,2)-2));
        Axx = spdiags(1 + vx(ii)*Cxx,0,Axx);
        Axx = spdiags([vx(ii)*Nxx(2:end); 0],-1,Axx);
        Axx = spdiags([0; vx(ii)*Sxx(1:end-1)],1,Axx);
        Axx = spdiags([zeros(dim_yl,1); vx(ii)*Exx(1:end-dim_yl)],dim_yl,Axx);
        Axx = spdiags([vx(ii)*Wxx(dim_yl+1:end); zeros(dim_yl,1)],-dim_yl,Axx);
        
        % Compute right side of equation
        
        Cbxx = 1 + ux(ii)*Cbxx;
        Nbxx = [ux(ii)*Nbxx(2:end); 0];
        Sbxx = [0; ux(ii)*Sbxx(1:end-1)];
        Ebxx = [zeros(dim_yl,1); ux(ii)*Ebxx(1:end-dim_yl)];
        Wbxx = [ux(ii)*Wbxx(dim_yl+1:end); zeros(dim_yl,1)];
        
        % Multiply right side with known field values at »kz«
        
        Cbxx = Pbx(glob_adr_slgs)          .* Cbxx;
        Nbxx = Pbx(glob_adr_slgs - 1)      .* Nbxx;
        Sbxx = Pbx(glob_adr_slgs + 1)      .* Sbxx;
        Ebxx = Pbx(glob_adr_slgs + dim_y)  .* Ebxx;
        Wbxx = Pbx(glob_adr_slgs - dim_y)  .* Wbxx;
        
        % Form vector for right side
        
        bxx = sparse(Cbxx + Nbxx + Sbxx + Ebxx + Wbxx);

        %% Solve SoLE
        
        [Px_l,~] = bicgstab(Axx,bxx,solver_tol);
       
        % Reshape solution to match field matrix
        
        Pbx = zeros(dim_y,dim_x); 
        Pbx(2:end-1,2:end-1) = reshape(Px_l,dim_yl,dim_xl);
        
        %% Apply absorber if specified
        
        if isnumeric(ABSORBER) && (ABSORBER > 0) && (ABSORBER < 1) && (ABSORBER ~= 0)
        
            n_threshold             = n_min + delta_n * ABSORBER;
            adr_n_threshold         = find(squeeze(n(:,:,kz)) <= n_threshold);
            Pbx(adr_n_threshold)    = 0;
            
        elseif ~isnumeric(ABSORBER) || (ABSORBER < 0) || (ABSORBER > 1)
            
            out = 'Invalid specification of absorber: 0 <= absorber < 1.';
            disp(out)
            return
            
        end

    end
    
    %% Merge solution in global field matrix
    
    Px(:,:,kz+1) = Pbx;
    
    % Update progress bar
    
    if strcmp(PROGRESS,'bar') == 1
    
        step_percent = 1;  % Defines step size for progress bar update

        if floor(100*kz/(size(n,3)-1)) > c*step_percent 

          if ishandle(h)    % Check if waitbar has not been closed to prevent error

            minutes_passed      = floor(toc/60);
            seconds_passed      = toc - minutes_passed*60;
            minutes_remaining   = floor((toc/c)*((100/step_percent)-c)/60);
            seconds_remaining   = (toc/c)*((100/step_percent)-c) - minutes_remaining*60;

            waitbar(kz/(size(n,3)-1),h,['Progress: ' num2str(c*step_percent) '%. Time remaining: ' num2str(floor(minutes_remaining)) 'm ' num2str(ceil(seconds_remaining)) 's.']) 

            c = c + 1;

          else 

             disp('Process interrupted') 
             break

          end 
        end 
   
    elseif strcmp(PROGRESS,'cl') == 1
        
        step_percent = 10;  % Defines step size for progress bar update

        if floor(100*kz/(size(n,3)-1)) > c*step_percent 

            minutes_passed      = floor(toc/60);
            seconds_passed      = toc - minutes_passed*60;
            minutes_remaining   = floor((toc/c)*((100/step_percent)-c)/60);
            seconds_remaining   = (toc/c)*((100/step_percent)-c) - minutes_remaining*60;

            out = ['   Progress: ' num2str(c*step_percent) '%. Time remaining: ' num2str(floor(minutes_remaining)) 'm ' num2str(ceil(seconds_remaining)) 's.'];
            disp(out)

            c = c + 1;
            
        end
        
        
    end
    
end

try
    close(h)
end

minutes_passed     = floor(toc/60);
seconds_passed    = toc - minutes_passed*60;

out = ['BPM-calculation finished after: ' num2str(minutes_passed) 'm ' num2str(seconds_passed) 's.'];
disp(out);

end
