classdef LensOptimizer < Lens
    % Class that helps with optimization
    
    properties
        z1
        z2
        antenna_r
    end
    
    methods
        function obj = LensOptimizer(L, lambda, dx, z1, z2, antenna_r)
            obj@Lens(L, lambda, dx)
            % Constructor
            obj.z1 = z1;
            obj.z2 = z2;
            obj.antenna_r = antenna_r;
        end
        
        function u1 = idealcpu1(obj, xangle, yangle)
            source = obj.pso(0,0 , obj.z1);
            lens = obj.makecplens(obj.z1, obj.z2, obj.antenna_r, 0);
            center = obj.lenspropagate(source, lens, obj.z1, obj.z2);
            xshiftblocks = fix(obj.z2*sind(xangle)/obj.dx);
            yshiftblocks = fix(obj.z2*sind(yangle)/obj.dx);
            u1 = circshift(center, [yshiftblocks, xshiftblocks]);            
        end
        
        function u1 = normdb(~, u1)
            u1 = abs(u1);
            u1 = mag2db(u1) - max(max(u1));
        end
        
        function error = field2error(obj, u1, ideal, dblim)
            u1 = obj.normdb(u1);
            u1(u1 < dblim) = dblim;
            error = immse(ideal, u1);
        end
    end
end

