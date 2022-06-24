function [Sp, Taumax,teta2p, teta2s, tetap] = principalstresses_eigen(S)
% DESCRIPTION:
% This function calculates for each node of the mesh the principal stresses
% and their orientation using the stress components xx, yy and xy defined
% on each node of the 2D mesh.
% INPUT VARIABLES:
% S = matrix containg the stress components xx, yy and xy for each
% node, it is a N x 3 matrix.
% OUTPUT VARIABLES:
% SP = matrix containg principal stress 1 and 2 for each node, it is N x 2 
% dimension;
% Taumax = array containg the maximum values of the shear stress component
% for each node;
% teta2p, teta2s = array containg the angle 2theta that defines the
% orientation of the principal stress components (teta2p) or the maximum
% shear stress (teta2s) for each node of the mesh;
% tetap = array containg the angle 2theta that defines the orientation of 
% the principal stress components resect to the horizontal (Sxx) for each
% node.

for i =1 : size(S, 1)
    
    [V,D]=eig([S(i,1),S(i,3);S(i,3),S(i,2)]);%Magitude
    Sp(i,1) = D(2,2)
    
    Sp(i,2) = D(1,1)
    
    Taumax(i,1)= sqrt(((S(i,1)-S(i,2))/2 )^2 +....
        S(i,3)^2);%Tau max
    
    %Orientation
    teta2p(i,1) = atan(2*S(i,3)/((S(i,1)-S(i,2)) ));%possible orientation of Sigma1 in radians expressed as 2teta (see page 209 of book 'Stuctural Geology 4th ed.' of Ragan D.M.)
    teta2s(i,1) = atan(-((S(i,1)-S(i,2)) )/2*S(i,3));%posible orientation of Tau max in radians expressed as 2teta
    %the above angles correspon to 2teta. The oriantion could be teta or
    %teta+pi (180°). To understand the orientation we must rotate the sigma x
    %and according the two angle (teta and teta+pi). The angle that gives a
    %rotated sigma equals to the magintude of sigma 1 is the correct one.
    %%%Srot(i,1)= (S(i,1)+S(i,2))/2 + ((S(i,1)-S(i,2))/2)*cos(teta2p(i,1)) +S(i,3)* sin(teta2p(i,1));
    %Srot(i,1) = Sp(i,1) * cos(teta2p(i,1)/2)^2 + Sp(i,2) * sin(teta2p(i,1)/2)^2;
    %Srot(i,1) = S(i,1) * cos(teta2p(i,1)/2)^2 + S(i,2) * sin(teta2p(i,1)/2)^2 + S(i,3)* cos(teta2p(i,1)/2) * sin(teta2p(i,1)/2);
    %if Sp(i,1) ==((S(i,1)+S(i,2))/2 + cos(teta2p(i,1))*(S(i,1)-S(i,2))/2 + sin(teta2p(i,1))*S(i,3))
    Srot(i,1)= (S(i,1)+S(i,2))/2 + ((S(i,1)-S(i,2))/2)*cos(teta2p(i,1)) +S(i,3)* sin(teta2p(i,1));
    if abs(Sp(i,1)-Srot(i,1))<abs(Sp(i,2)-Srot(i,1))
        
        tetap(i,1) = teta2p(i,1)/2;
        tetap(i,2) = teta2p(i,1)/2 + pi/2;
    else
        tetap(i,1) = teta2p(i,1)/2 + pi/2;
        tetap(i,2) = teta2p(i,1)/2;
    end
end