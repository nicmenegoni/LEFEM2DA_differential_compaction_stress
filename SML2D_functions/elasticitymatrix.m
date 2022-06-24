function [D] = elasticitymatrix(E,nu, plainstrain)
% DESCRIPTION:
% This fucntion define the elasticity matrix;
if plainstrain==1 %Plain strain formulation
    D = (E / ((1+nu)*(1 - 2*nu)) )* ...
        [ 1-nu, nu,            0;
        nu,  1-nu,            0;
        0,  0, (1 - 2*nu)/2];
elseif plainstrain==2 %PlainSTRESS formulation
    D = (E / (1-nu^2) )* ...
        [ 1-nu, nu,            0;
        nu,  1-nu,            0;
        0,  0, (1 - nu)/2];
end
end