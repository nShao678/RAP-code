function err = lanczos_generalized(A,M,v,iterMax,lambda,tol)

    L = ichol(M);
    T = zeros(iterMax+1, iterMax+1);
    vv = zeros(size(v));
    v = v / sqrt(v'*M*v); % Normalize the initial vector
    err = zeros(1,iterMax);
    beta = 0;
    for j = 1:iterMax
        w = A*v;
        alpha = v'*w;
        T(j,j) = alpha;
        err(j) = (eigs(T(1:j,1:j),1,'smallestabs')-lambda)/lambda;
        if err(j)<tol
            err = err(1:j);
            break
        end
        ww = w;
        [w,~] = pcg(M,w,1e-12,20,L,L',v);
        w = w-alpha*v-beta*vv;
        beta = sqrt(w'*ww);
        if ~isreal(beta)
            1
        end
        vv = v;
        v = w/beta;
        T(j+1,j) = beta;
        T(j,j+1) = beta;
    end
end
