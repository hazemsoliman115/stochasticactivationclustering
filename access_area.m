classdef access_area
    
    properties
        x_pos
        y_pos
        queue_load
        qos_factor
    end
    
    methods
        function obj = access_area(x,y,t,q)
            if nargin > 0
                obj.x_pos = x;
                obj.y_pos = y;
                obj.queue_load = t;
                obj.qos_factor = q;
            end
        end
    end
end