clc; clear; close all


%% ground truth
A = rand(10,10).*randi(100) + 1i.*rand(10,10).*randi(100);

[V,D] = eig(A);

%% max eigenvalue
v = (1+1i).*ones(10,1);
v = v./norm(v);
max_itr = 10000;
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

while norm(l_itr(k)-l_itr(k-1)) >= tol && k <= max_itr
    k = k+1;
    w = A*v;
    v = w./norm(w);
    l_itr(k) = v'*A*v;
end

figure(1)

plot((real(l_itr)), (imag(l_itr)), '-s')
hold on
plot([real(l_itr(end))-5 real(l_itr(end))+5], [imag(l_itr(end)) imag(l_itr(end))], 'k')
plot([real(l_itr(end)) real(l_itr(end))], [imag(l_itr(end))-5 imag(l_itr(end))+5], 'k')
scatter(real(l_itr(end)), imag(l_itr(end)))
xlim([min(real(l_itr)), max(real(l_itr))])
ylim([min(imag(l_itr)), max(imag(l_itr))])
xlabel('Re\{$\lambda$\}', 'interpreter', 'latex')
ylabel('Im\{$\lambda$\}', 'interpreter', 'latex')
title(['Power Iteration'])

%% all 10 eigenvalues
%Vrq = [v zeros(10,9)];
Vrq = zeros(10,10);
Drq = zeros(10,10);
%Drq(1,1) = l_itr(end);
max_itr = 1000;
tol = 1e-12;
eps = 0.1;

fig2 = figure(2);
set(fig2,'Position',[100 100 1000 1000])
nRows = 3; nCols = 4;

for i = 1:10
    clear l_itr v
    
    %{
    v = rand(10,1);
    
    for j = 1:i-1
        v = v-(v'*Vrq(:,i-1)).*v;
    end
    
    v = v./norm(v);
    %}
    
    v = V(:,i) + eps.*(2.*rand(10,1)-1);
    
    v'*Vrq(:,1:i-1)
    
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
    
    while norm(l_itr(k)-l_itr(k-1)) >= tol && k <= max_itr
        k = k+1;
        w = (A-l_itr(k-1).*eye(10))\v;
        v = w./norm(w);
        l_itr(k) = v'*A*v;
    end
    
    Drq(i,i) = l_itr(k);
    Vrq(:,i) = v;
    
    figure(2)
    subplot(nRows,nCols,i)
    plot((real(l_itr)), (imag(l_itr)), '-s')
    hold on
    plot([real(l_itr(end))-5 real(l_itr(end))+5], [imag(l_itr(end)) imag(l_itr(end))], 'k')
    plot([real(l_itr(end)) real(l_itr(end))], [imag(l_itr(end))-5 imag(l_itr(end))+5], 'k')
    scatter(real(l_itr(end)), imag(l_itr(end)))
    xlim([min(real(l_itr)), max(real(l_itr))])
    xlabel('Re\{$\lambda^k$\} - Re\{$\lambda$\}', 'interpreter', 'latex')
    ylabel('Im\{$\lambda^k$\} - Im\{$\lambda$\}', 'interpreter', 'latex')
    title(['eigenvalue ', num2str(i)])
end
