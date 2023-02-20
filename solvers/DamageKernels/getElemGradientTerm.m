function [r_grad, ke_grad] = getElemGradientResTerm(msh,mat,gauss,local_dissipation,d,e)
    
    %coeficients
    Gc = mat.Gc;
    l = mat.l;
    c0 = local_dissipation.c0;
    
    % 2d quad element residual vector routine
    ke_grad   = zeros(4,4);
    r_grad   = zeros(4,1);
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

        eta = gauss(i);  psi = gauss(j);

        % compute derivatives of shape functions in reference coordinates
        NJpsi = 0.25*psiJ.*(one + eta*etaJ);
        NJeta = 0.25*etaJ.*(one + psi*psiJ);

        % compute derivatives of x and y wrt psi and eta
        xpsi = NJpsi*xe'; ypsi = NJpsi*ye';
        xeta = NJeta*xe'; yeta = NJeta*ye';
        Jinv = [yeta, -xeta; -ypsi, xpsi]';
        jcob = xpsi*yeta - xeta*ypsi;

        % compute derivatives of shape functions in element coordinates
        NJdpsieta = [NJpsi; NJeta];
        NJdxy = Jinv*NJdpsieta;

        % assemble B matrix
        BJ = NJdxy;
        
        %get grad_d at quadrature point
        grad_d_qp = BJ * de'; %u at quadrature point

        % assemble ke
        r_grad = r_grad + (2*Gc*l/c0) * (BJ' * grad_d_qp) / jcob;
        ke_grad = ke_grad + (2*Gc*l/c0) * (BJ'*BJ) / jcob;
        
      end
    end

end
              