function err = rap(A,M,B,x0,mu,L,iterMax,lambda,tol)
err = zeros(1,iterMax);
rq = @(x) x'*A*x/(x'*M*x);
kappa = sqrt(mu/L); 
beta = 3*kappa/(2-4*kappa); 
alpha = (sqrt(beta^2+4*(beta+1)*kappa^2)-beta)/2;
gamma = alpha*mu/(alpha+beta);
c1 = alpha/(alpha+beta+1); 
c2 = (1-alpha)/alpha; 
c3 = alpha/(1+beta)/gamma;
x = B(x0);
eta = sqrt(x'*x0); x = x/eta; x0 = x0/eta;
v = x; v0 = x0;

for iter = 1:iterMax
    
    err(iter) = (rq(x)-lambda)/lambda;
    if err(iter) < tol
        break
    end
    % update y
    eta = (x'*v0+v'*x0)/2;
    w = v-x*eta;
    if norm(w) < 100*eps
        sigma = rq(x);
        g0 = (A*x-sigma*(M*x))*2/(x'*M*x);
        g = B(g0);
        q = -c3*g;
        q0 = -c3*g0;
        eta = sqrt(q'*q0);
        q = q/eta;
        q0 = q0/eta;
        v = v*cos(eta)+q*sin(eta);
        v0 = v0*cos(eta)+q0*sin(eta);
        [U,R] = qr([x,g],0);
        [c,~] = eigs(U'*A*U,U'*M*U,1,'smallestabs');
        x = U*c;
        c = R\c;
        x0 = x0*c(1)+g0*c(2);
    else
        theta = c1*acos(eta);
        if ~isreal(theta)
            theta = 0;
        end
        w0 = v0-x0*eta;
        eta = sqrt(w'*w0);
        w = w/eta;
        w0 = w0/eta;
        y = x*cos(theta)+w*sin(theta);
        y0 = x0*cos(theta)+w0*sin(theta);
    
        % update v
        sigma = rq(y);
        eta = (y'*v0+v'*y0)/2;
        p = v-y*eta;
        p0 = v0-y0*eta;
        eta = sqrt(p'*p0);
        p = p/eta;
        p0 = p0/eta;
        g0 = (A*y-sigma*(M*y))*2/(y'*M*y);
        g = B(g0);
        q = c2*theta*p-c3*g;
        q0 = c2*theta*p0-c3*g0;
        eta = sqrt(q'*q0);
        q = q/eta;
        q0 = q0/eta;
        v = y*cos(eta)+q*sin(eta);
        v0 = y0*cos(eta)+q0*sin(eta);
        % update x
        [U,R] = qr([x,y,g],0);
        [c,~] = eigs(U'*A*U,U'*M*U,1,'smallestabs');
        x = U*c;
        c = R\c;
        x0 = x0*c(1)+y0*c(2)+g0*c(3);
    
    end
    eta = sqrt(x'*x0);
    x = x/eta;
    x0 = x0/eta;
end
err = err(1:iter);


end

