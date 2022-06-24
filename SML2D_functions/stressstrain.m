function [S,Eps] = stressstrain(number_of_nodes, t,displacements, D, B, material)
% t =trangular element
% displacements = nodal displacement
% D = Flexural rigidity
% B = strain displacement matrix
% material = defined material for each element
%% Calculate stress xx, yy and zz components
S = zeros(1 * number_of_nodes, 3);
node_occurences = zeros(1 * number_of_nodes, 1);
for element = 1 : size(t,1)
    nodes = t(element, :);
    displacement_nodes_of_element = displacements(:, nodes);
    U_e = displacement_nodes_of_element(:);
    s = cell2mat(D(material(element)))* B{element} * U_e;
    S(nodes', :) = S(nodes', :) + repmat(s', 3, 1);
end
%% Calculate strain xx, yy and zz components
Eps = zeros(1 * number_of_nodes, 3);
node_occurences = zeros(1 * number_of_nodes, 1);
for element = 1 : size(t,1)
    nodes = t(element, :);
    displacement_nodes_of_element = displacements(:, nodes);
    U_e = displacement_nodes_of_element(:);
    eps = B{element} * U_e;
    Eps(nodes', :) = Eps(nodes', :) + repmat(eps', 3, 1);
end
%Determine node frequencies (occurrences)
for i = 1 : number_of_nodes
    node_occurences(i) = numel(find(t == i));
end
S = S ./ repmat(node_occurences, 1, 3);
Eps = Eps ./ repmat(node_occurences, 1, 3);
end