function r_grad = getElemGradientResTerm(msh,mat,gauss,local_dissipation,d,e)
    
    %coeficients
    Gc = mat.Gc;
    l = mat.l;
    c0 = local_dissipation.c0;
    
    % 2d quad element residual vector routine
    r_grad   = zeros(1,4);
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

        % assemble ke
        r_grad = r_grad + (Gc*l/c0) * (BJ'*BJ) * de / jcob;
        
      end
    end

end
              