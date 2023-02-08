%class that implements the consitutive equations of a simple, linear
%elastic solid
classdef solidModel < handle 
    
    properties
        mesh
        youngsModulus
        poissonsRatio
    end
    
    methods
        %constructor 
        function obj = solidModel()
        end

        %compute strain
        function strain = getStrain(elem,qp,u) 
            % 2d quad element stiffness matrix routine
            one  = ones(1,4);
            psiJ = [-1, +1, +1, -1]; etaJ = [-1, -1, +1, +1];

            % get coordinates of element nodes
            for j=1:4
              je = obj.mesh.ien(elem,j+1); xe(j) = obj.mesh.xnod(je,2); ye(j) = obj.mesh.xnod(je,3);
              ue(2*(j-1)+1) = u(2*(je-1)+1);
              ue(2*(j-1)+2) = u(2*(je-1)+2);
            end

            % compute element stiffness
            %go from qp=2*(i-1)+j to i,j
            i = ceil(qp/2);
            j = qp - 2*(i-1);

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

            % compute strain 
            strain = BJ*ue'/jcob;
            
        end 
        
        %compute stress
        function stress = getStress(elem,qp,u)
            E = obj.youngsModulus.getValue(elem);
            nu = obj.poissonsRatio.getValue(elem);
            D = (E/(1+nu)*(1-2*nu))*[1-nu nu 0; nu 1-nu 0; 0 0 (1-2*nu)/2];
            stress = D*getStrain(elem,qp,u);
            %D =  lambda * [1 1 0; 1 1 0; 0 0 0] + mu * [2 0 0; 0 2 0; 0 0 1]; 
        end
           
        %compute stiffness
        function stiffness = getStiffness(~,~,~) 
            E = obj.youngsModulus.getValue(elem);
            nu = obj.poissonsRatio.getValue(elem);
            stiffness = (E/(1+nu)*(1-2*nu))*[1-nu nu 0; nu 1-nu 0; 0 0 (1-2*nu)/2];
        end
        
    end
    
end