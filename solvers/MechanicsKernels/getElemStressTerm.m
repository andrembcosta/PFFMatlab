function [r, J] = getElemStressTerm(msh,gauss,e,constitutiveModel,u)
    
    % 2d quad element residual vector routine
    J   = zeros(4,4);
    r   = zeros(4,1);
    one  = ones(1,4);
    psiJ = [-1, +1, +1, -1]; etaJ = [-1, -1, +1, +1];
    
    % get coordinates of element nodes
    xe = zeros(1,4); ye = zeros(1,4); 
    for j=1:4
      je = msh.ien(e,j+1); xe(j) = msh.xnod(je,2); ye(j) = msh.xnod(je,3);
    end
    
    % compute element stiffness
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
        BJ = zeros(3,8);
        BJ(1,1:2:7) = NJdxy(1,1:4);  BJ(2,2:2:8) = NJdxy(2,1:4);
        BJ(3,1:2:7) = NJdxy(2,1:4);  BJ(3,2:2:8) = NJdxy(1,1:4);

        % residual 
        r = r + BJ' * constitutiveModel.getStress(e,i,j,u);
        J = r + BJ' * constitutiveModel.getStiffness(e,i,j,u) * BJ / jcob;

      end
    end


end
              