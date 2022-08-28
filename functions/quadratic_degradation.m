classdef quadratic_degradation
    methods
        function y = value(obj,x)
            y = (1-x)^2;
        end
        function y = first_derivative(obj,x)
            y = 2*(x-1);
        end
        function y = second_derivative(obj, x)
            y = 2;
        end
    end
end