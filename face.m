clc; clear; close all

dn1 = './CroppedYale/';
dn2 = './yalefaces_uncropped/yalefaces/';
dn_fig = './figs/';

dir1 = dir(fullfile(dn1,'yaleB*'));

%% read + reshape (CROPPED)

% 1. For face space of ONE PERSON, set 'i_subdir' and uncomment 
%    'for ss = i_subdir:i_subdir' for the loop condition
% 2. For face space of ALL PEOPLE, uncomment 'for ss = 1:length(dir1)'


A = []; % data matrix
i_subdir = 1;

for ss = 1:length(dir1)
%for ss = i_subdir:i_subdir
    disp(ss)
    sdn = dir1(ss).name;
    sdir = dir(fullfile([dn1,sdn],'*.pgm'));
    
    % compile all images in each subdir
    for ff = 1:length(sdir)
        fn = sdir(ff).name;
        full_fn = [dn1,sdn,'/',fn];
        
        A_ff = double(imread(full_fn, 'pgm')); % image matrix
        A_ff_c = reshape(A_ff, length(A_ff(:)), 1); % reshape into column
        A = [A, A_ff_c]; % append to A
        clear A_ff_c
    end
    
    % display average image
    figure(1)
    imshow(uint8(reshape(mean(A, 2), size(A_ff))))
    title(gca, 'average face of selected images')
    
    clear A_ff
end

disp('A compiled');

%% dominant evs & svd
[U,S,V] = svd(A, 'econ');

m = size(A,1); n = size(A,2);
C = A*A';

v = ones(m,1);
v = v./norm(v);
max_itr = 1000;
tol = 1e-6;

% k = 1
w = C*v;
v = w./norm(w);
l_itr(1) = v'*C*v;

% k = 2
w = C*v;
v = w./norm(w);
l_itr(2) = v'*C*v;

k = 2;

while abs(l_itr(k)-l_itr(k-1)) >= tol && k <= max_itr
    k = k+1;
    w = C*v;
    v = w./norm(w);
    l_itr(k) = v'*C*v;
end

figure
plot(1:length(l_itr), sqrt(l_itr), [1, length(l_itr)], [max(diag(S)) max(diag(S))])
xlabel('# iteration')
ylabel('$\lambda$', 'interpreter', 'latex')
legend('square root of eigenvalue of $A^T A$', 'leading order of S', 'interpreter', 'latex')

figure
subplot(1,2,1)
imshow(reshape(v.*100,192,168))
title('power iteration dominant mode')

subplot(1,2,2)
imshow(reshape(abs(U(:,1)).*100,192,168))
title('SVD dominant mode')

%% random
k = 20;
Omega = randn(n,k);

Y = A*Omega;

[Q,R] = qr(Y, 0);

B = Q'*A;
[U_, S_rec, V_rec] = svd(B, 'econ');

U_rec = Q*U_;

%% modes
fig1=figure(1)
set(fig1,'Position',[100 100 1000 1000])

nRows = 3; nCols = 3;

for pp = 1:(nRows*nCols-1)
    subplot(nRows,nCols,pp)
    imshow([reshape(abs(U_rec(:,pp)).*100,192,168), ...
        reshape(abs(U(:,pp)).*100,192,168)])
end

subplot(nRows,nCols,nRows*nCols)
semilogy(1:nRows*nCols, diag(S_rec(1:nRows*nCols,1:nRows*nCols)));
hold on
semilogy(1:nRows*nCols, diag(S(1:nRows*nCols,1:nRows*nCols)));
xlim([1,nRows*nCols])
xlabel('eigenvalue index')
ylabel('eigenvalue')
legend('random','true modes')

%% singular value decay vs k
K = floor(logspace(1, 2, 5));
legend_str = cell(1,5);

figure(2)

for kk = 1:5
    k = K(kk);
    legend_str(kk) = {['k=', num2str(K(kk))]};
    
    Omega = randn(n,k);
    
    Y = A*Omega;
    
    [Q,R] = qr(Y, 0);
    
    B = Q'*A;
    [~, S_rec, ~] = svd(B, 'econ');
    
    figure(2)
    hold on
    semilogy(1:nRows*nCols, diag(S_rec(1:nRows*nCols,1:nRows*nCols)));
end

hold on
semilogy(1:nRows*nCols, diag(S(1:nRows*nCols,1:nRows*nCols)), 'k', 'linewidth', 1.5);

legend({legend_str{:}, 'SVD'})
xlabel('singular value index')
ylabel('singular value')
