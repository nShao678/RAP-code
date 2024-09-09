function err = pinvit(A,M,T,x,iterMax,lambda,tol)
    err = zeros(1,iterMax);
    x = x/normm(x);
    rho = (x'*A*x);
    for iter = 1:iterMax
        err(iter) = (rho-lambda)/lambda;
        if err(iter)<tol
            err = err(1:iter);
            break
        end
        r = (A*x-rho*M*x);
        g = T(r);
        V = [x,g];
        V = orth(V);
        Av = V'*A*V;
        Mv = V'*M*V;
        [c,~] = eigs(Av,Mv,1,'smallestabs');
        x = V*c;
        x = x/normm(x);
        rho = x'*A*x;
    end
    function y = normm(x)
        y = sqrt(x'*M*x);
    end
end
