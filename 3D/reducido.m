function [v] = reducido(H,e)

a = area(H,e);
vol = volumen(H);

Ro = sqrt(a/4*pi);
v = vol/(4*pi/3)*Ro^3;