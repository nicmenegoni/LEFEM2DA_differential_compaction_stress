%%%%%close all; clear; clc
% FEM Analysis in Plain Strain -> it works only for thick plate
%% Mesh reading
%%%%edge = e(1, :);
%%%%t = t(1:3, :)';
tic
t=tria;% triplets of nodes defining the triangular elements (N.elements x3 array)
p=vert';% x and y coordinates of the model nodes (2xN.nodes array)
material=tnum;% element material (N.elements x 1 array)
number_of_nodes = size(p, 2);%
number_of_elements = size(t, 1);

%% Calculation Parameters
%General parameters
g=9.81;%[m/s^2] gravitational acceleration
consider_bodyforce = 1; % 0 = no body forces, 1 = yes body forces
APPLY2MAT=2;%material where the body force is applied; 1= Basin, 2=1st Platform; 3 =2nd platform
apply_distribforce = 0;% 0= no distribuited forces, 1 = yes distribuited forces
apply_nodalforce = 0;% 0= no nodal/concentrated forces, 1 = yes nodal/concentrated forces

%% Material parameters
E1 = 0.1e9;%GPa
nu1 = 0.15;
D1 = elasticitymatrix(E1, nu1,1);%Flexural rigidity or elasticity matrix
rho1=0.22e4;%Kg/m^3

%Material parameters
E2 = 0.1e11;%GPa
nu2 = 0.3;
D2 = elasticitymatrix(E2, nu2,1);%Flexural rigidity or elasticity matrix
rho2=0.271e4;%Kg/m^3

%Material parameters
E3 = 01e11;%GPa
nu3 = 0.3;
D3 = elasticitymatrix(E3, nu3,1);%Flexural rigidity or elasticity matrix
rho3=0.271e4;%Kg/m^3

rho=[rho1,rho2,rho3];
D={D1,D2,D3};

%% Initialization of K and F
K = zeros(2 * number_of_nodes);
F = zeros(2 * number_of_nodes, 1);

%% Calculation of Ke, Fe & assembly of K and F
[K,B] = stiffnessmatrix (t,p, material, D, K);

%% Forces definition
if consider_bodyforce == 1 %Applying body forces
[F] = bodyforce(t, p, g, rho, material, F, APPLY2MAT);
end

t_Neumann = [];
if apply_distribforce ==1 %Applying distribuited force
%1) Neumann boundary (where the force is applied?)
    for e = 1 : number_of_elements
        nodes = t(e, :);
        
        I = p(1, nodes) == 0 | p(1, nodes) == 1000 ;%This should apply the force on the last x value?
        if( sum(I) == 2)
            t_Neumann = [t_Neumann; nodes(I)];
        end
        
    end
%2) Neumann 'forces'
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

if apply_nodalforce == 1
    nodes=ismember(p',[500,35],rows)
dofs_Neumann = reshape([2 * nodes - 1; 2 * nodes], 1, 2 * numel(nodes));
end

%% Dirichlet boundary conditions (This define the points with no movement?)
for nodes = 1 : size(vert, 1)
    % type 1 --> fixed
    % type 2 --> roller on y-axis (fixed x-dir)
    % type 3 --> roller on x-axis (fixed y-dir)
    if p(1, nodes) == 0 && p(2, nodes) == 0 
        [K,F] = boundarycondition (nodes, K, F, 1);
    elseif p(2, nodes) == 0  %Roller Boundary (blocked onto y)
        [K,F] = boundarycondition (nodes, K, F, 3);
    elseif p(1,nodes) == 0 %Roller Boundary (blocked onto x)
        [K,F] = boundarycondition (nodes, K, F, 2);
    end
end

%% Solve
U = K \ F;
displacements = [U(1 : 2 : end), U(2 : 2 : end)]';%defined at the node
magnification =1;%5e2;
p_new = p + magnification * displacements;

[S,Eps] = stressstrain(number_of_nodes, t,displacements, D, B, material);

%% Calculate principal components of stress and strain
%Calculate the magnitude and orientation of the principal stresses
[Sp, Taumax,teta2p, teta2s, tetap] = principalstresses(S);

disp('Jub runned in ', num2str(toc/60), 'minutes')


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


save('slm2d_1plat_mesh','vert','p','p_new','tria','tnum','material','etri','number_of_nodes')
save('slm2d_1plat_results','p_new','S','Sp', 'Taumax','teta2p', 'teta2s', 'tetap','Eps','displacements')