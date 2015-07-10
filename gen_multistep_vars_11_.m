function [u,v]=gen_multistep_vars_11_(dz, alpha, beta_z)

u = zeros(1,1);
v = zeros(1,1);

R=2*beta_z;

C=1/R^2-1j*dz*(1-alpha)*1/R;
D=1/R^2+1j*dz*( alpha )*1/R;

u(1)= C;

v(1)= D;

end