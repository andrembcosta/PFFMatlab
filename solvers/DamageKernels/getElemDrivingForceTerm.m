function [r_driv, ke_driv] = getElemDrivingForceTerm(msh, gauss, degradation, active_energy ,d, e)

    % 2d quad element mass matrix routine
    r_driv = zeros(4,1);
    ke_driv  = zeros(4,4);
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
        
        %convert qp from 1:2 x 1:2 to 1:4
        qp_ind = 2*(i-1) + j;

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
        %residual term
        r_driv = r_driv + degradation.first_derivative(d_qp) * active_energy(qp_ind) * NJ' * jcob;
        %jacobian term
        ke_driv = ke_driv + degradation.second_derivative(d_qp) * active_energy(qp_ind) * (NJ' * NJ) * jcob;

      end
    end


end