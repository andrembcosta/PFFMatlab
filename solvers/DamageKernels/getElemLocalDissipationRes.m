function r_d = getElemLocalDissipationRes(msh,mat,gauss,local_dissipation,e)

    %coeficients
    Gc = mat.Gc;
    l = mat.l;
    c0 = local_dissipation.c0;
    
    % 2d quad element residual vector routine
    r_d   = zeros(1,4);
    one  = ones(1,4);
    psiJ = [-1, +1, +1, -1]; etaJ = [-1, -1, +1, +1];
    
    % get coordinates of element nodes
    xe = zeros(1,4); ye = zeros(1,4); de = zeros(1,4); 
    for j=1:4
      je = msh.ien(e,j+1); xe(j) = msh.xnod(je,2); ye(j) = msh.xnod(je,3);
      de(j) = d(je);
    end
    
    % compute element residual
    for i=1:2
      for j=1:2

        eta = gauss(i); psi = gauss(j);

        % compute derivatives of shape functions in reference coordinates
        NJpsi = 0.25*psiJ.*(one + eta*etaJ);
        NJeta = 0.25*etaJ.*(one + psi*psiJ);
        NJ    = 0.25*(one + psi*psiJ).*(one + eta*etaJ);

        % compute derivatives of x and y wrt psi and eta
        xpsi = NJpsi*xe'; ypsi = NJpsi*ye';
        xeta = NJeta*xe'; yeta = NJeta*ye';
        jcob = xpsi*yeta - xeta*ypsi;
        
        %get d at quadrature point
        d_qp = NJ * de'; %u at quadrature point

        % apply quadrature (with unit weights)
        r_d = r_d + (Gc/(c0*l)) * local_dissipation.first_derivative(d_qp) * NJ' * jcob;

      end
    end
    
end