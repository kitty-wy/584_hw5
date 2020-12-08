clc; clear; close all

%% (a) ground truth
A = zeros(10,10);

for i = 1:10
    A(i,i) = rand*randi(10);
    
    for j = (i+1):10
        A(i,j) = rand*randi(10);
        A(j,i) = A(i,j);
    end
end

[V,D] = eig(A);

%% (b) max eigenvalue
v = ones(10,1)./sqrt(10); %????
max_itr = 1000;
tol = 1e-12;

% k = 1
w = A*v;
v = w./norm(w);
l_itr(1) = v'*A*v;

% k = 2
w = A*v;
v = w./norm(w);
l_itr(2) = v'*A*v;

k = 2;

while abs(l_itr(k)-l_itr(k-1)) >= tol && k <= max_itr
    k = k+1;
    w = A*v;
    v = w./norm(w);
    l_itr(k) = v'*A*v;
end

figure(1)
subplot(1,2,1)
plot(1:length(l_itr), l_itr)
xlim([1,length(l_itr)])
xlabel('# iteration')
ylabel('$\lambda^{(k)}$', 'interpreter','latex')
title(['Power Iteration: max eigenvalue'])

subplot(1,2,2)
plot(1:length(l_itr), abs(l_itr-max(diag(D))))
xlim([1,length(l_itr)])
xlabel('# iteration')
ylabel('$|\lambda^{(k)}-\lambda|$', 'interpreter','latex')
title(['Power Iteration: absolute error'])

%% (c) all 10 eigenvalues
Vrq = [v zeros(10,9)];
%Vrq = zeros(10,10);
%Drq = zeros(10,10);
Drq(1,1) = l_itr(end);
max_itr = 1000;
tol = 1e-12;
eps = 0.2;

fig2 = figure(2);
fig3 = figure(3);
set(fig2,'Position',[100 100 1500 1000])
set(fig3,'Position',[100 100 1500 1000])
nRows = 3; nCols = 4;

for i = 1:10
    clear l_itr v
    
    v = V(:,i) + eps.*(2.*rand(10,1)-1);
    
    l0 = v'*A*v;
        
    % k = 1
    w = (A-l0.*eye(10))\v;
    v = w./norm(w);
    l_itr(1) = v'*A*v;
    
    % k = 2
    w = (A-l_itr(1).*eye(10))\v;
    v = w./norm(w);
    l_itr(2) = v'*A*v;
    
    k = 2;
    
    while l_itr(k)-l_itr(k-1) >= tol && k <= max_itr
        k = k+1;
        w = (A-l_itr(k-1).*eye(10))\v;
        v = w./norm(w);
        l_itr(k) = v'*A*v;
    end
    
    Drq(i,i) = l_itr(k);
    Vrq(:,i) = v;
    
    figure(2)
    subplot(nRows,nCols,i)
    plot(1:length(l_itr), l_itr)
    xlim([1,length(l_itr)])
    ylim([min(l_itr), max(l_itr)])
    xlabel('# iteration')
    ylabel('$\lambda^{(k)}$', 'interpreter','latex')
    title(['eigenvalue ', num2str(i)])
    
    figure(3)
    subplot(nRows,nCols,i)
    plot(1:length(l_itr), abs(l_itr-l_itr(end)))
    xlim([1,length(l_itr)])
    xlabel('# iteration')
    ylabel('$|\lambda^{(k)}-\lambda|$', 'interpreter','latex')
    title(['eigenvalue ', num2str(i)])
    
end
