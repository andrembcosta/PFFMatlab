classdef quadratic
    properties
        c0 = 2;
    end
    methods
        function y = value(obj,x)
            y = x^2;
        end
        function y = first_derivative(obj,x)
            y = 2*x;
        end
        function y = second_derivative(obj, x)
            y = 2;
        end
    end
end