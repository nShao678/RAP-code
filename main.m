warning('off');
addpath(genpath(pwd))
sympref('FloatingPointOutput',true);
global T
iterMax = 20000;
err = cell(4,8);
N = 2;
for n = 1:8
    [K,M,T,x0] = genDD(n,N);
    lambda = eigs(K,M,1,'smallestabs');

    tol = 1e-6;
    tic
    mu = 8;
    L = 50;
    err{1,n} = rap(K,M,T,x0,mu,L,iterMax,lambda,tol);
    [~,~,~,rho] = lobpcg(x0,K,M,'Tinv',tol);
    rho = (rho-lambda)/lambda;
    idx = find(rho<tol,1);
    err{2,n} = rho(1:idx);
    err{3,n} = pinvit(K,M,T,x0,iterMax,lambda,tol);
    err{4,n} = lanczos_generalized(K,M,x0,iterMax,lambda,tol);
    toc

end
IterNum = zeros(size(err));
for ii = 1:size(IterNum,1)
    for jj = 1:size(IterNum,2)
        IterNum(ii,jj) = length(err{ii,jj});
    end
end
latex(sym(IterNum))

