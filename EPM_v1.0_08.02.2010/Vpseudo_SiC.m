
function [V]=Vpseudo_SiC(G)
%This function calculates the pseudopotential
%corresponding to a certain G vector
%This is only in the case of i~=j (not equal)
%The function argument is the reciprocal lattice vector G
%The Form Factors are in this file.

DotProduct=dot(G(:),G(:));
if DotProduct==3  
    Vs = (-0.4280)*13.6059;  %13.6 is the Rydberg energy in eV's
    Va = 0.0010*13.6059;
elseif DotProduct==4
    Vs = 0;
    Va = 0.0800*13.6059;
elseif DotProduct==8  
    Vs = 0.1017*13.6059;
    Va = 0;
elseif DotProduct==11 
    Vs = 0.1108*13.6059;
    Va = 0.0277*13.6059;
else
    Vs = 0;
    Va = 0;
end;


T=[1/8 1/8 1/8];

V=Vs*cos(2*pi*dot(G(:),T(:))) + Va*j*sin(2*pi*dot(G(:),T(:)));


