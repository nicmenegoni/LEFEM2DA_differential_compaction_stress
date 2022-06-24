function [K,F] = boundarycondition (node, K, F, type)
% DESCRIPTION:
%This function definesthe boundary condition for the analyzed node.
% node = analyzed node;
% K = stifness matrix;
% F = forces matrix;
% type = type of boundary condition:
% if type = 1, is a fixed node; 
% if type = 2, is x-displacement is blocked;
% if type = 3, is y-displacement is blocked.
if type ==1
        dofs = [2 * node - 1; 2 * node];
        K(dofs,:) = 0;
        K(dofs, dofs) = eye(numel(dofs));
        F(dofs) = 0;
elseif type ==2
            dofs = [2 * node - 1; 2 * node];
        K(dofs(1,:), :) = 0;
        K(dofs(1,:), dofs(1,:)) = eye(numel(dofs(1,:)));
        F(dofs(1,:)) = 0;
elseif type ==3
        dofs = [2 * node - 1; 2 * node];
        K(dofs(2,:),:) = 0;
        K(dofs(2,:), dofs(2,:)) = eye(numel(dofs(2,:)));
        F(dofs(2,:)) = 0;
end