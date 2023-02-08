%class that implements the consitutive equations of a simple, linear
%elastic solid
classdef damageSolidModel < solidModel
    
    properties
        mesh
        criticalFractureEnergy
        regularizationLength
        decomposition
        degradationFunction
        localDissipation
    end
    
    methods
        %constructor
        function obj = damageSolidModel()
        end
        
        %compute stress
        function stress = getStress(elem, qp, u, d)
            if obj.decomposition=='none'
                E = obj.youngsModulus.getValue(elem);
                nu = obj.poissonsRatio.getValue(elem);
                D = (E/(1+nu)*(1-2*nu))*[1-nu nu 0; nu 1-nu 0; 0 0 (1-2*nu)/2];
                stress = obj.degradationFunction(d)*D*getStrain(elem,qp,u);            
            end
            if obj.decomposition=='spectral'
                %TODO: implement this
                E = obj.youngsModulus.getValue(elem);
                nu = obj.poissonsRatio.getValue(elem);
                D = (E/(1+nu)*(1-2*nu))*[1-nu nu 0; nu 1-nu 0; 0 0 (1-2*nu)/2];
                stress = D*getStrain(elem,qp,u);            
            end
        end
           
        %compute stiffness
        function stiffness = getStiffness(~,~,~,d) 
            if obj.decomposition == 'none'
                E = obj.youngsModulus.getValue(elem);
                nu = obj.poissonsRatio.getValue(elem);
                stiffness = obj.degradationFunction(d)*(E/(1+nu)*(1-2*nu))*[1-nu nu 0; nu 1-nu 0; 0 0 (1-2*nu)/2];
            end
            if obj.decomposition == 'spectral'
                %TODO: implement this
            end
        end
        
        %compute active energy
        function activeEnergy = getActiveEnergy(elem, qp, u, d)
            activeEnergy = 0;
            %TODO: inneficient repeated calls to getStrain and getStrain
            if obj.decomposition == 'none'
                activeEnergy = 0.5*obj.getStress(elem, qp, u, d)*obj.getStrain(elem, qp, u);
            end
            if obj.decomposition == 'spectral'
                %TODO: implement this
            end
        end
        
        %compute nonlocal microforce
        function nonLocal = getNonlocalMicroforce(elem, qp, d)
            ell = obj.regularizationLegth.getValue(elem);
            Gc = obj.criticalFractureEnergy.getValue(elem);
            c0 = obj.localDissipation.getC0();
            nonLocal = Gc*ell/c0 * grad_d(elem, qp, d);
        end
        
        %compute internal microforce
        function local = getLocalMicroforce(elem, qp, u, d)
            ell = obj.regularizationLegth.getValue(elem);
            Gc = obj.criticalFractureEnergy.getValue(elem);
            c0 = obj.localDissipation.getC0();
            local = obj.degradation.firstDerivative(d) * getActiveEnergy(elem, qp, u, d) ...
                    + Gc/(ell*c0)*obj.localDissipation.firstDerivative(d);
        end
        
    end
    
end