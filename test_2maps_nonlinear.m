clc 
clear all

% DEFINE THE FAMILY OF CONTINUOUS, HOMOGENEOUS and ORDER-PRESERVING MAPS

% 1ST EXAMPLE
alpha_1=0.3;
alpha_2=0.5;
A=[0.8,0.1;0.1,0.8];
B=[0.8,0.2;0,0.8];
C=[0.9,0;0.8,0.7];
D=[0.8,0;0.1,0.8];

F_1=@(x) (C*((D*x).^(alpha_1))).^(1/alpha_1);
F_2=@(x) (B*((A*x).^(alpha_2))).^(1/alpha_2);
F_family = {F_1,F_2};


% 2ND EXAMPLE
% alpha_1=3;
% alpha_2=5;
% A=[0.8,0.1;0,0.8];
% B=[0.8,0.2;0,0.8];
% C=[0.9,0;0.8,0.7];
% D=[0.8,0;0.1,0.8]; 
% F_1=@(x) (C*((D*x).^(alpha_1))).^(1/alpha_1);
% F_2=@(x) (B*((A*x).^(alpha_2))).^(1/alpha_2);
% F_family = {F_1,F_2};


% FIND DOMINANT EIGENVECTOR AND EIGENVALUE OF THE PRESUMED SPECTRAL MAXIMIZING
% MAP and NORMALIZE THE EIGENVALUE TO THE PRESUMED JSR

% Compute spectral radius of F_1
[x,r]= pw_method(F_1,2, 1000, 1e-15);
% Proposed JSR
rho = r^(1/1);

% Run polytopal prenorm algorithm
[vertices,rho,info] = polytope_prenorm(F_family,x,rho,30);

% Display results
disp('Computed JSR:');
disp(rho);
disp("Number of vertices:");
disp(size(vertices,2))
disp('Information:');
disp(info);

% DISPLAY STEP BY STEP CONSTRUCTION OF THE EXTREMAL PRENORM
% normalize maps in the family
F_norm = cell(1,size(F_family,2));
for q = 1:size(F_norm,2)
    F_norm{q} = @(x) F_family{q}(x) / rho;
end

% number of vertices of the extremal prenorm
N_vertices = size(vertices,2);

%%% DISPLAY THE EXTREMAL PRENORM
% [sortisci,Isort]=sort(Y(1,:));
% Y=Y(:,Isort);
% G=zeros(2,2*(N_vertices)+2);
% G(:,2:2:2*(N_vertices))=Y(:,1:N_vertices);
% G(2,1)=Y(2,1);
% G(1,2*(N_vertices)+1)=Y(1,N_vertices);
% for i=3:2:2*(N_vertices)-1
%    G(:,i)=min(G(:,i-1),G(:,i+1));
% end
% pgon=polyshape(G(1,:),G(2,:));
% figure
% plot(pgon)


%%% DISPLAY THE CONSTRUCTION OF THE EXTREMAL PRENORM, STEP BY STEP
for k=1:N_vertices
   % for each k<=N_vertices display a new figure with the extremal prenorm
   % induced by the first k vertices
   Y=vertices(:,1:k); 
   x_0=Y(:,k);
   [sortisci,Isort]=sort(Y(1,:));
   Y=Y(:,Isort);
   G=zeros(2,2*(k)+2);
   G(:,2:2:2*(k))=Y(:,1:k);
   G(2,1)=Y(2,1);
   G(1,2*(k)+1)=Y(1,k);
   for i=3:2:2*(k)-1
       G(:,i)=min(G(:,i-1),G(:,i+1));
   end
   pgon=polyshape(G(1,:),G(2,:));
f = figure('Color','w','Position',[100 100 800 700]);

ax = axes(f,'Position',[0.15 0.15 0.75 0.75]);

plot(ax,pgon)
axis(ax,[0 0.3 0 1.2])
hold(ax,'on')
% Display the k-th vertex
eps = 1e-2;
plot(ax,x_0(1),x_0(2), 'o', 'Color','red', 'MarkerSize',12, 'MarkerFaceColor','#F77100')
text(ax,x_0(1)+eps,x_0(2)+eps,'x', ...
    'FontSize',18, 'Color','k', 'HorizontalAlignment','left', 'VerticalAlignment','bottom')
title(ax,sprintf('{%d}-th step of polytopal algorithm',k), 'FontSize',18, 'Color','k')

% Display the image of the k-th vertex through the action
% of the maps in the family
for j = 1:size(F_norm,2)
    x = F_norm{j}(x_0);
    plot(ax,x(1),x(2), ...
        'hexagram', ...
        'Color','black', 'MarkerSize',12, 'MarkerFaceColor','#00728D')
    text(ax,x(1)+eps,x(2)+eps, ...
        sprintf('f_{%d}(x)',j), ...
        'FontSize',18, 'Color','k', 'HorizontalAlignment','left','VerticalAlignment','bottom')
end

% Make sure there is enough space for axes, tick labels and title
set(ax,'LooseInset',get(ax,'TightInset'))
drawnow
% Save the entire figure, including axes and title
exportgraphics(f, fullfile("figures",sprintf("%d-step_prenorm.png",k)), 'Resolution',300);
end



