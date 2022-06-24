function [F] = bodyforce(t, p, g, rho, material, F, APPLY2MATERIAL)
% DESCRIPTION:
% This function defines the body forces for each 

for element = 1 : size(t, 1)%in this for cicle the force are applied to selected nodes
        %Body forces (e.g. gravity or inertial force)
        if sum(material(element) == APPLY2MATERIAL) > 0
        nodes = t(element, :);
        dofs = reshape([2 * nodes - 1; 2 * nodes], 1, 2 * numel(nodes));
        area_of_element = abs(det([ones(1, 3); p(:, nodes)]))/2;
        Hb = 1/3 * [1,0; 0,1;1,0; 0,1;1,0; 0,1];
        gForce= -g*rho(material(element));%
        Fb = Hb * [0; gForce] * area_of_element;
        F(dofs) = F(dofs) + Fb;
        end
end
end