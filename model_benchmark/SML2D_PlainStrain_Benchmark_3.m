%%%%%close all; clear; clc
% FEM Analysis in Plain Strain -> it works only for thick plate
%% Mesh reading
%%%%edge = e(1, :);
%%%%t = t(1:3, :)';
t=tria;
p=vert';
number_of_nodes = size(p, 2);
number_of_elements = size(t, 1);

%% Parameters
consider_gravity=1;%1 yes, 0 no
g=9.81;%[m/s^2] gravitational acceleration
material=tnum;

E1 = 0.1e11;%GPa
nu1 = 0.3;
D1 = (E1 / ((1+nu1)*(1 - 2*nu1)) )* ...
    [ 1-nu1 nu1            0;
    nu1  1-nu1            0;
    0  0 (1 - 2*nu1) / 2];%Flexural rigidity or elasticity matrix
rho1=0.271e4;%Kg/m^3

rho=[rho1];
D={D1};
%% Initialization of K and F
K = zeros(2 * number_of_nodes);
F = zeros(2 * number_of_nodes, 1);

%% Calculation of Ke, Fe & assembly of K and F
%Ke is the element coefficient matrix?
%phi the vector of primary unknown quantities normally expressed by the quantities at the nodes?
%B is the boundary nodal quantities
%see page 32 of https://www.sjsu.edu/me/docs/hsu-Chapter%2011%20Finite%20element%20analysis_04-25-19.pdf
for element = 1 : number_of_elements
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
            diff_Phi(i, 2)  diff_Phi(i, 1)];
        
        B{element} = [B{element}, b_e];
    end
    
    
    Ke = B{element}' * cell2mat(D(material(element))) * B{element} * area_of_element;%element stifness matrix
    dofs = reshape([2 * nodes - 1; 2 * nodes], 1, 2 * numel(nodes));
    K(dofs, dofs) = K(dofs, dofs) + Ke;%Stiffness matrix?
    
end

%% Neumann boundary (where the force is applied?)
t_Neumann = [];
if consider_gravity==0
    for e = 1 : number_of_elements
        nodes = t(e, :);
        
        I = p(1, nodes) == 0 | p(1, nodes) == 1000 ;%This should apply the force on the last x value?
        if( sum(I) == 2)
            t_Neumann = [t_Neumann; nodes(I)];
        end
        
    end
end
% if consider_gravity==0
%     for nodes = 1 : number_of_nodes
%         if p(1, i) == 0 || p(1, i) == 1000
%             t_Neumann = [t_Neumann; nodes];
%         end
%     end
% end
%% Forces definition

if consider_gravity == 1
    for element = 1 : size(t, 1)%in this for cicle the force are applied to selected nodes
        %Body forces (e.g. gravity or inertial force)
        nodes = t(element, :);
        P = [ones(1, 3); p(:, nodes)];
        C = inv(P);
        area_of_element = abs(det(P))/2;
        diff_Phi = C(:, 2:3);
        dofs = reshape([2 * nodes - 1; 2 * nodes], 1, 2 * numel(nodes));
        area_of_element = abs(det([ones(1, 3); p(:, nodes)]))/2;
        Hb = 1/3 * [1,0; 0,1;1,0; 0,1;1,0; 0,1];
        gForce= -g*rho(material(element));%
        Fb = Hb * [0; gForce] * area_of_element;
        F(dofs) = F(dofs) + Fb;
        %Fb = B{element}' * [0; gForce] * B{element} * area_of_element;
        
    end
end

if size(t_Neumann, 1)>0
    for element = 1 : size(t_Neumann, 1)
        nodes = t_Neumann(element, :);
        dofs_Neumann = reshape([2 * nodes - 1; 2 * nodes], 1, 2 * numel(nodes));
        P = p(:, nodes);
        length_of_element = norm(diff(P, 1, 2));
        H_mean = 1/2 * repmat(eye(2), 1, 2);
        if P(1,:)==0
        q = [-5e7; 0];
        elseif P(1,:)==1000
        q = [5e7; 0];
        end
        Fe = H_mean' * q * length_of_element;
        F(dofs_Neumann) = F(dofs_Neumann) + Fe;
    end
end
% % % % % % % %% Dirichlet boundary conditions (This define the points with no movement?)
% % % % % % % for nodes = 1 : size(vert, 1)
% % % % % % %     %nodes = etri(element, :);
% % % % % % %     
% % % % % % %     if p(2, nodes) == 0 && p(1,nodes) ==500 %Fixed Boundary (blocked onto y)
% % % % % % %         dofs = [2 * nodes - 1; 2 * nodes];
% % % % % % %         K(dofs,:) = 0;
% % % % % % %         K(dofs, dofs) = eye(numel(dofs));
% % % % % % %         F(dofs) = 0;
% % % % % % %     elseif p(2, nodes) == 0  %Roller Boundary (blocked onto y)
% % % % % % %         dofs = [2 * nodes - 1; 2 * nodes];
% % % % % % %         K(dofs(2,:),:) = 0;
% % % % % % %         K(dofs(2,:), dofs(2,:)) = eye(numel(dofs(2,:)));
% % % % % % %         F(dofs(2,:)) = 0;
% % % % % % % %     elseif p(1,nodes) == 2000%Roller Boundary (blocked onto x)
% % % % % % % %         dofs = [2 * nodes - 1; 2 * nodes];
% % % % % % % %         K(dofs(1,:), :) = 0;
% % % % % % % %         K(dofs(1,:), dofs(1,:)) = eye(numel(dofs(1,:)));
% % % % % % % %         F(dofs(1,:)) = 0;
% % % % % % %      end
% % % % % % % end

%% Solve
U = K \ F;
displacements = [U(1 : 2 : end), U(2 : 2 : end)]';
magnification =1;%5e2;
p_new = p + magnification * displacements;

%% Calculate stress xx, yy and zz components
S = zeros(1 * number_of_nodes, 3);
node_occurences = zeros(1 * number_of_nodes, 1);
for element = 1 : number_of_elements
    nodes = t(element, :);
    displacement_nodes_of_element = displacements(:, nodes);
    U_e = displacement_nodes_of_element(:);
    s = cell2mat(D(material(element)))* B{element} * U_e;
    S(nodes', :) = S(nodes', :) + repmat(s', 3, 1);
end

%% Calculate strain xx, yy and zz components
Eps = zeros(1 * number_of_nodes, 3);
node_occurences = zeros(1 * number_of_nodes, 1);
for element = 1 : number_of_elements
    nodes = t(element, :);
    displacement_nodes_of_element = displacements(:, nodes);
    U_e = displacement_nodes_of_element(:);
    eps = B{element} * U_e;
    Eps(nodes', :) = Eps(nodes', :) + repmat(eps', 3, 1);
end
%% Determine node frequencies (occurrences)
for i = 1 : number_of_nodes
    node_occurences(i) = numel(find(t == i));
end
S = S ./ repmat(node_occurences, 1, 3);
Eps = Eps ./ repmat(node_occurences, 1, 3);
%% Calculate princiapl components of stress and strain
%Calculate the magnitude and orientation of the principal stresses
[Sp, Taumax,teta2p, teta2s, tetap] = principalstresses(S)
%% Plotting
%Stress (Sigma)
figure(1)
subplot(3,1,1)
STRESS_COMPONENT = 1;
splot = S(:, STRESS_COMPONENT);
set(gcf, 'color', 'w');colormap parula; hold on
trisurf(t, p_new(1, :), p_new(2, :), zeros(1, number_of_nodes), splot, 'EdgeColor', 'k', 'FaceColor', 'interp');
trisurf(t, p(1, :), p(2, :), zeros(1, number_of_nodes), 'EdgeColor',  0.75 * [1 1 1], 'FaceColor', 'none');
colorbar;view(2);axis equal;axis tight
title('\sigma XX')

subplot(3,1,2)
STRESS_COMPONENT = 2;
splot = S(:, STRESS_COMPONENT);
set(gcf, 'color', 'w');colormap parula; hold on
trisurf(t, p_new(1, :), p_new(2, :), zeros(1, number_of_nodes), splot, 'EdgeColor', 'k', 'FaceColor', 'interp');
trisurf(t, p(1, :), p(2, :), zeros(1, number_of_nodes), 'EdgeColor',  0.75 * [1 1 1], 'FaceColor', 'none');
colorbar;view(2);axis equal;axis tight
title('\sigma YY')

subplot(3,1,3)
STRESS_COMPONENT = 3;
splot = S(:, STRESS_COMPONENT);
set(gcf, 'color', 'w');colormap parula; hold on
trisurf(t, p_new(1, :), p_new(2, :), zeros(1, number_of_nodes), splot, 'EdgeColor', 'k', 'FaceColor', 'interp');
trisurf(t, p(1, :), p(2, :), zeros(1, number_of_nodes), 'EdgeColor',  0.75 * [1 1 1], 'FaceColor', 'none');
colorbar;view(2);axis equal;axis tight
title('\sigma XY')

%Strain (Epsilon)
figure(2)
subplot(3,1,1)
STRAIN_COMPONENT = 1;
epsplot = Eps(:, STRAIN_COMPONENT);
set(gcf, 'color', 'w');colormap parula; hold on
trisurf(t, p_new(1, :), p_new(2, :), zeros(1, number_of_nodes), epsplot, 'EdgeColor', 'k', 'FaceColor', 'interp');
trisurf(t, p(1, :), p(2, :), zeros(1, number_of_nodes), 'EdgeColor',  0.75 * [1 1 1], 'FaceColor', 'none');
colorbar;view(2);axis equal;axis tight
title('\epsilon XX')

subplot(3,1,2)
STRAIN_COMPONENT = 2;
epsplot = Eps(:, STRAIN_COMPONENT);
set(gcf, 'color', 'w');colormap parula; hold on
trisurf(t, p_new(1, :), p_new(2, :), zeros(1, number_of_nodes), epsplot, 'EdgeColor', 'k', 'FaceColor', 'interp');
trisurf(t, p(1, :), p(2, :), zeros(1, number_of_nodes), 'EdgeColor',  0.75 * [1 1 1], 'FaceColor', 'none');
colorbar;view(2);axis equal;axis tight
title('\epsilon YY')

subplot(3,1,3)
STRAIN_COMPONENT = 3;
epsplot = Eps(:, STRAIN_COMPONENT);
set(gcf, 'color', 'w');colormap parula; hold on
trisurf(t, p_new(1, :), p_new(2, :), zeros(1, number_of_nodes), epsplot, 'EdgeColor', 'k', 'FaceColor', 'interp');
trisurf(t, p(1, :), p(2, :), zeros(1, number_of_nodes), 'EdgeColor',  0.75 * [1 1 1], 'FaceColor', 'none');
colorbar;view(2);axis equal;axis tight
title('\epsilon XY')

%Principal stresses plot
figure(3)
subplot(3,1,1)
STRESS_COMPONENT = 1;
splot = Sp(:, STRESS_COMPONENT);
set(gcf, 'color', 'w');colormap parula; hold on
trisurf(t, p_new(1, :), p_new(2, :), zeros(1, number_of_nodes), splot, 'EdgeColor', 'k', 'FaceColor', 'interp');
trisurf(t, p(1, :), p(2, :), zeros(1, number_of_nodes), 'EdgeColor',  0.75 * [1 1 1], 'FaceColor', 'none');
colorbar;view(2);axis equal;axis tight
title('\sigma 1')

subplot(3,1,2)
STRESS_COMPONENT = 2;
splot = Sp(:, STRESS_COMPONENT);
set(gcf, 'color', 'w');colormap parula; hold on
trisurf(t, p_new(1, :), p_new(2, :), zeros(1, number_of_nodes), splot, 'EdgeColor', 'k', 'FaceColor', 'interp');
trisurf(t, p(1, :), p(2, :), zeros(1, number_of_nodes), 'EdgeColor',  0.75 * [1 1 1], 'FaceColor', 'none');
colorbar;view(2);axis equal;axis tight
title('\sigma 2')

subplot(3,1,3)
splot = Taumax(:, 1);
set(gcf, 'color', 'w');colormap parula; hold on
trisurf(t, p_new(1, :), p_new(2, :), zeros(1, number_of_nodes), splot, 'EdgeColor', 'k', 'FaceColor', 'interp');
trisurf(t, p(1, :), p(2, :), zeros(1, number_of_nodes), 'EdgeColor',  0.75 * [1 1 1], 'FaceColor', 'none');
colorbar;view(2);axis equal;axis tight
title('\tau max')

%Principal stresses plot
figure(4)
subplot(3,1,1)
STRESS_COMPONENT = 1;
splot = Sp(:, STRESS_COMPONENT)*1e-6;;
set(gcf, 'color', 'w');colormap parula; hold on
trisurf(t, p_new(1, :), p_new(2, :), zeros(1, number_of_nodes), splot, 'EdgeColor', 'k', 'FaceColor', 'interp');
quiver(p_new(1, :)', p_new(2, :)', cos(tetap(:,1)), sin(tetap(:,1)),0.3,'k')
c=colorbar;view(2);axis equal;axis tight
c.Label.String = 'MPa';
title('\sigma 1')

subplot(3,1,2)
STRESS_COMPONENT = 2;
splot = Sp(:, STRESS_COMPONENT)*1e-6;
set(gcf, 'color', 'w');colormap parula; hold on
trisurf(t, p_new(1, :), p_new(2, :), zeros(1, number_of_nodes), splot, 'EdgeColor', 'k', 'FaceColor', 'interp');
quiver(p_new(1, :)', p_new(2, :)', cos(tetap(:,2)), sin(tetap(:,2)),0.3,'k')
c=colorbar;view(2);axis equal;axis tight
c.Label.String = 'MPa';
title('\sigma 2')

subplot(3,1,3)
splot = Taumax(:, 1)*1e-6;;
set(gcf, 'color', 'w');colormap parula; hold on
trisurf(t, p_new(1, :), p_new(2, :), zeros(1, number_of_nodes), splot, 'EdgeColor', 'k', 'FaceColor', 'interp');
quiver(p_new(1, :)', p_new(2, :)', cos(teta2s(:,1)/2), sin(teta2s(:,1)/2),0.3,'k')
c=colorbar;view(2);axis equal;axis tight
c.Label.String = 'MPa';
title('\tau max')