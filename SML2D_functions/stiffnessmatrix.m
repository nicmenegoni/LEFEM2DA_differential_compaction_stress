function [K,B] = stiffnessmatrix (t,p, material, D,K)
% DESCRIPTION:
% Calculation of Ke, Fe & assembly of K and F
% Ke = is element stiffness matrix.
% phi = the vector of primary unknown quantities normally expressed by the 
% quantities at the nodes.
% B = is the boundary nodal quantities
% K = is the model stiffness matrix
% B = strain displacement matrix
% see page 32 of https://www.sjsu.edu/me/docs/hsu-Chapter%2011%20Finite%20element%20analysis_04-25-19.pdf
for element = 1 : size(t,1)
    nodes = t(element, :);
    P = [ones(1, 3); p(:, nodes)];
    C = inv(P);
    area_of_element = abs(det(P))/2;
    diff_Phi = C(:, 2:3);
    
    B{element} = [];%initializing the strain displacement matrix
    for i = 1 : 3
        b_e = [
            diff_Phi(i, 1)  0;
            0               diff_Phi(i, 2);
            diff_Phi(i, 2)  diff_Phi(i, 1)];%element strain displacement matrix
        
        B{element} = [B{element}, b_e];
    end
    
    
    Ke = B{element}' * cell2mat(D(material(element))) * B{element} * area_of_element;%element stifness matrix
    dofs = reshape([2 * nodes - 1; 2 * nodes], 1, 2 * numel(nodes));
    K(dofs, dofs) = K(dofs, dofs) + Ke;%Model stiffness matrix
    
end
end