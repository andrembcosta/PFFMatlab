classdef quadratureField < handle
    
    properties
       names
       values
    end
    
    methods
        function obj = quadratureField()
        end
        
        function field = getField(field_name)
            field = [];
            for i = 1:size(obj.names,1)
               if obj.names(i)==field_name
                    field = obj.values(i);
               end
            end      
        end
    
    end

    
end
