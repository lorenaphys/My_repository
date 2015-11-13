function [V] = volumen(H)

V = sum(sum(sum(H > .99)));
