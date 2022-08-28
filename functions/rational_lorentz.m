classdef rational_lorentz
    properties
        p;
        m;
    end
    methods
        function obj = rational_lorentz(p,m)
            obj.p = p;
            obj.m = m;
        end
        function y = value(obj,x)
            y = (1-x)^2/( (1-x)^2 + obj.m*x*(1+obj.p*x) );
        end
        function y = first_derivative(obj,x)
            num = - ( (1-x)^2*obj.m*(1+2*obj.p*x) + 2*(1-x)*obj.m*x*(1+obj.p*x) );
            den = ( (1-x)^2 + obj.m*x*(1+obj.p*x) )^2;
            y = num/den;
        end
        function y = second_derivative(obj, x)
            den = ( (1-x)^2 + obj.m*x*(1+obj.p*x) )^4;
            num_part_1 = - ( -2*(1-x)*obj.m*(1+2*obj.p*x) + 2*(1-x)^2*obj.m*obj.p + 2*(1-2*x)*obj.m*(1+obj.p*x) + 2*obj.p*obj.m*x*(1-x) ) * ( (1-x)^2 + obj.m*x*(1+obj.p*x) )^2;
            num_part_2 = ( (1-x)^2*obj.m*(1+2*obj.p*x) + 2*(1-x)*obj.m*x*(1+obj.p*x) ) * 2 * ( (1-x)^2 + obj.m*x*(1+obj.p*x) ) * (-2*(1-x) + obj.m*(1+obj.p*x)+obj.p*obj.m*x);
            num = num_part_1 + num_part_2;
            y = num/den;
        end
    end
end