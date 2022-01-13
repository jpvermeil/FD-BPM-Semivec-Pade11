function [ aY ] = padeY( n,yg,dim_xl,dim_yl,localAdrSlgs,globalAdrSlgs,POLARIZATION,FIELDCOMPONENT)
% This function populaltes the y coefficient vector that is merged into the
% system matrix A as part of the respective main diagonal of for each step
% in z-direction.
%
% SYNOPSIS
% 
% [ aY ] = padeY(n,xg,dim_y,dim_xl,dim_yl,localAdrE,globalAdrE,POLARIZATION,FIELDCOMPONENT)
%
% VARIABLES
%
%   n               Refractive index profile
%   xg              Grid of x dimension. Has to match dimensions of n. 
%                   Has to be X output of Matlab meshgrid function ([X,Y,Z] = meshgrid(x,y,z))
%   yg              Grid of y dimension. Has to match dimensions of n.
%                   Has to be Y output of Matlab meshgrid function ([X,Y,Z] = meshgrid(x,y,z))
%   dim_y           Size of the global matrix in y-direction
%   dim_xl          Size of the local matrix in x-direction
%   dim_yl          Size of the local matrix in y-direction
%   localAdrSlgs    Contains the local address information of all elements
%   globalAdrSlgs   Contains the global address information of all elements
%   POLARIZATION    String value that can either be 'TE' or 'TM'
%   FIELDCOMPONENTS Can be 'Ex', 'Ey', 'Hx' or 'Hy'

% This function populaltes the coefficient vector that is merged into the
% system matrix A as part of the respective main diagonal of for each step
% in z-direction. It belongs to the program "Semi vectorial wide angle
% Padé(1,1) finite difference BPM for TE/TM E- and/or H-fields in 3D
% structures".
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
    
    dn = yg(globalAdrSlgs - 1) - yg(globalAdrSlgs);
    ds = yg(globalAdrSlgs) - yg(globalAdrSlgs + 1);
    
    % Round off error
    dn = 1e-12*round(dn*1e12);
    ds = 1e-12*round(ds*1e12);

    aNorth = padeN(n,yg,dim_xl,dim_yl,localAdrSlgs,globalAdrSlgs,POLARIZATION,FIELDCOMPONENT);
    aSouth = padeS(n,yg,dim_xl,dim_yl,localAdrSlgs,globalAdrSlgs,POLARIZATION,FIELDCOMPONENT);
    
    if (strcmp('TE',POLARIZATION)) == 1
        
        aY = - aNorth - aSouth; 
        
    elseif (strcmp('TM',POLARIZATION) && (strcmp('Hx',FIELDCOMPONENT) || strcmp('Hy',FIELDCOMPONENT))) == 1
        
        aY = - aNorth - aSouth; 
    
    elseif (strcmp('TM',POLARIZATION) && (strcmp('Ex',FIELDCOMPONENT) || strcmp('Ey',FIELDCOMPONENT))) == 1
        
        aY = -4./(dn.*ds) + aNorth + aSouth;
        
    else
        
        out = ['WARNING: Unknown field component or polarization. Possible choices are: ''Ex'', ''Ey'', ''Hx'', ''Hy'' or ''TE'' and ''TM'' respectively.'];
        disp(out);
       
    end
   
end

