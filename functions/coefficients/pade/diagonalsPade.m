function [ diagC,diagN,diagS,diagE,diagW ] = diagonalsPade( beta_0,neff,n,xg,yg,dim_y,dim_xl,dim_yl,grid,gammaBoundaryCondition,POLARIZATION,FIELDCOMPONENT,BC)
% This function populaltes the vectors that are merged into the system
% matrix A as the respective diagonals of for each step in z-direction.
%
% SYNOPSIS
% 
% [diagC,diagN,diagS,diagE,diagW] = diagonalsPade(beta_0,neff,n,xg,yg,dim_y,dim_xl,dim_yl,grid,gammaBoundaryCondition,POLARIZATION,FIELDCOMPONENT,BC)
% 
% VARIABLES
% 
%   beta_0          wave number
%   neff            Propagation constant (For wide angle Pade BPM this is usually (n1 + n2) / 2)
%   n               Refractive index profile
%   xg              Grid of x dimension. Has to match dimensions of n. 
%                   Has to be X output of Matlab meshgrid function ([X,Y,Z] = meshgrid(x,y,z))
%   yg              Grid of y dimension. Has to match dimensions of n.
%                   Has to be Y output of Matlab meshgrid function ([X,Y,Z] = meshgrid(x,y,z))
%   dim_y           Size of the global matrix in y-direction
%   dim_xl          Size of the local matrix in x-direction
%   dim_yl          Size of the local matrix in y-direction
%   grid            Contains the information of element numbers and boundaries
%       .boundaryFlag           Matrix that flags boundary elements with a '1'
%       .globalIndex            Global element numbers
%       .localIndex             Local element numbers
%   gammaBoundaryCondition
%                   Defines how the boundary condition is handled
%                   For 'ABC' this is '0'
%                   For 'TBC' it contains 4 Vectors of quotients of the field values at the respective boundaries
%   POLARIZATION    String value that can either be 'TE' or 'TM'
%   FIELDCOMPONENTS Can be 'Ex', 'Ey', 'Hx' or 'Hy'
%   BC              String value speficying boundary condition. Currently only 'ABC' supported

% This function populaltes the vectors that are merged into the system
% matrix A as the respective diagonals of for each step in z-direction. It
% belongs to the program "Semi vectorial wide angle Padé(1,1) finite
% difference BPM for TE/TM E- and/or H-fields in 3D structures".
% Copyright (C) 2017 Jan-Philipp Roth (JanPhilipp.Roth@gmail.com) 
% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 3 of the License, or (at your
% option) any later version. This program is distributed in the hope that
% it will be useful, but WITHOUT ANY WARRANTY; without even the implied
% warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details. You should have received a
% copy of the GNU General Public License along with this program; if not,
% write to the Free Software Foundation, Inc., 51 Franklin Street, Fifth
% Floor, Boston, MA 02110-1301  USA

if ((strcmp('TE',POLARIZATION) || strcmp('TM',POLARIZATION)) && ((strcmp('Ex',FIELDCOMPONENT) || strcmp('Ey',FIELDCOMPONENT) || strcmp('Hx',FIELDCOMPONENT) || strcmp('Hy',FIELDCOMPONENT)))) == 0
        
    out = 'WARNING: Unknown field component or polarization. Possible choices are: ''Ex'', ''Ey'', ''Hx'', ''Hy'' or ''TE'' and ''TM'' respectively.';
    disp(out);
    return

end

if ((strcmp(BC,'ABC') || (strcmp(BC,'TBC')))) == 0
    
    out = 'WARNING: Unknown boundary condition. Possible choices are: ''ABC'' or ''TBC''';
    disp(out);
    return
    
end   

epsilon = n.*n;

%% Address space initialization
    
    globalIndexSlgs   = grid.globalIndex(2:end-1,2:end-1);
    globalAdrSlgs   = reshape(globalIndexSlgs,size(globalIndexSlgs,1)*size(globalIndexSlgs,2),1);

    globalIndexN = grid.globalIndex(3:end-1,2:end-1); 
    globalAdrN = reshape(globalIndexN,size(globalIndexN,1)*size(globalIndexN,2),1);
    
    globalIndexS = grid.globalIndex(2:end-2,2:end-1); 
    globalAdrS = reshape(globalIndexS,size(globalIndexS,1)*size(globalIndexS,2),1);

    globalIndexE = grid.globalIndex(2:end-1,2:end-2);
    globalAdrE = reshape(globalIndexE,size(globalIndexE,1)*size(globalIndexE,2),1);
    
    globalIndexW = grid.globalIndex(2:end-1,3:end-1);
    globalAdrW = reshape(globalIndexW,size(globalIndexW,1)*size(globalIndexW,2),1);
    
    localIndexSlgs    = grid.localIndex(2:end-1,2:end-1);
    localAdrSlgs    = reshape(localIndexSlgs,size(localIndexSlgs,1)*size(localIndexSlgs,2),1);
    
    localIndexN  = grid.localIndex(3:end-1,2:end-1); 
    localAdrN  = reshape(localIndexN,size(localIndexN,1)*size(localIndexN,2),1);
    
    localIndexS  = grid.localIndex(2:end-2,2:end-1); 
    localAdrS  = reshape(localIndexS,size(localIndexS,1)*size(localIndexS,2),1);
    
    localIndexE = grid.localIndex(2:end-1,2:end-2);
    localAdrE = reshape(localIndexE,size(localIndexE,1)*size(localIndexE,2),1);
    
    localIndexW = grid.localIndex(2:end-1,3:end-1);
    localAdrW = reshape(localIndexW,size(localIndexW,1)*size(localIndexW,2),1);
    
    if strcmp(BC,'TBC') == 1
        
        kx = 2:1:size(n,2)-1;
        ky = 2:1:size(n,1)-1;
        
        localAdrNTBC = grid.localIndex(2,kx)';
        localAdrETBC = grid.localIndex(ky,end-1)';
        localAdrSTBC = grid.localIndex(end-1,kx)';
        localAdrWTBC = grid.localIndex(ky,2)';

        globalAdrNTBC = grid.globalIndex(2,kx)';
        globalAdrETBC = grid.globalIndex(ky,end-1)';
        globalAdrSTBC = grid.globalIndex(end-1,kx)';
        globalAdrWTBC = grid.globalIndex(ky,2)';

        gammaN = ones(dim_xl*dim_yl,1);
        gammaE = ones(dim_xl*dim_yl,1);
        gammaS = ones(dim_xl*dim_yl,1);
        gammaW = ones(dim_xl*dim_yl,1);
        
        gammaN(localAdrNTBC) = gammaBoundaryCondition{1};
        gammaE(localAdrETBC) = gammaBoundaryCondition{2};
        gammaS(localAdrSTBC) = gammaBoundaryCondition{3};
        gammaW(localAdrWTBC) = gammaBoundaryCondition{4};
        
    end
    
%% Generation of diagonals
    
% The functions pade*() generate vectors matching the passed element
% addresses. The position of the calculated coefficients in the diagonal
% vector is contained in the according address vector.

    if strcmp(BC,'ABC')
    
        diagC = padeX(n,xg,dim_xl,dim_yl,dim_y,localAdrSlgs,globalAdrSlgs,POLARIZATION,FIELDCOMPONENT) + padeY(n,yg,dim_xl,dim_yl,localAdrSlgs,globalAdrSlgs,POLARIZATION,FIELDCOMPONENT) + beta_0^2 .* (epsilon(globalAdrSlgs) - neff^2);
        diagN = padeN(n,yg,dim_xl,dim_yl,localAdrN,globalAdrN,POLARIZATION,FIELDCOMPONENT);
        diagS = padeS(n,yg,dim_xl,dim_yl,localAdrS,globalAdrS,POLARIZATION,FIELDCOMPONENT);
        diagE = padeE(n,xg,dim_y,dim_xl,dim_yl,localAdrE,globalAdrE,POLARIZATION,FIELDCOMPONENT);
        diagW = padeW(n,xg,dim_y,dim_xl,dim_yl,localAdrW,globalAdrW,POLARIZATION,FIELDCOMPONENT);
        
    elseif strcmp(BC,'TBC')
        
        diagNTBC = gammaN .* padeN(n,yg,dim_xl,dim_yl,localAdrNTBC,globalAdrNTBC,POLARIZATION,FIELDCOMPONENT); 
        diagSTBC = gammaS .* padeS(n,yg,dim_xl,dim_yl,localAdrSTBC,globalAdrSTBC,POLARIZATION,FIELDCOMPONENT); 
        diagETBC = gammaE .* padeE(n,xg,dim_y,dim_xl,dim_yl,localAdrETBC,globalAdrETBC,POLARIZATION,FIELDCOMPONENT); 
        diagWTBC = gammaW .* padeW(n,xg,dim_y,dim_xl,dim_yl,localAdrWTBC,globalAdrWTBC,POLARIZATION,FIELDCOMPONENT); 
        
        diagTBC = (diagNTBC + diagSTBC + diagETBC + diagWTBC);
        
        diagC = padeX(n,xg,dim_xl,dim_yl,dim_y,localAdrSlgs,globalAdrSlgs,POLARIZATION,FIELDCOMPONENT) + padeY(n,yg,dim_xl,dim_yl,localAdrSlgs,globalAdrSlgs,POLARIZATION,FIELDCOMPONENT) + beta_0^2 .* (epsilon(globalAdrSlgs) - neff^2);
        diagN = padeN(n,yg,dim_xl,dim_yl,localAdrN,globalAdrN,POLARIZATION,FIELDCOMPONENT);
        diagS = padeS(n,yg,dim_xl,dim_yl,localAdrS,globalAdrS,POLARIZATION,FIELDCOMPONENT);
        diagE = padeE(n,xg,dim_y,dim_xl,dim_yl,localAdrE,globalAdrE,POLARIZATION,FIELDCOMPONENT);
        diagW = padeW(n,xg,dim_y,dim_xl,dim_yl,localAdrW,globalAdrW,POLARIZATION,FIELDCOMPONENT);
        
        diagC = (diagC + diagTBC);
        
    end
        
end

