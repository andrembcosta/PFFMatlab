classdef identity
    properties
        c0 = 8/3;
    end
    methods
        function y = value(obj,x)
            y = x;
        end
        function y = first_derivative(obj,x)
            y = 1;
        end
        function y = second_derivative(obj, x)
            y = 0;
        end
    end
end