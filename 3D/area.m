function [a] = area(H,e)

h = size(H);

Nx = h(1);
Ny = h(2);
Nz = h(3);

a = 3*sum(sum(sum(abs(gradient(H)).^2)))/(4*sqrt(2)*e*Nx*Ny*Nz);
