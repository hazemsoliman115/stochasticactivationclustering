classdef accesspoint
    
    properties
        x_pos
        y_pos
    end
    
    methods
        function obj = accesspoint(x,y)
            if nargin > 0
                obj.x_pos = x;
                obj.y_pos = y;
            end
        end
    end
end