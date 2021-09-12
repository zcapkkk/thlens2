classdef Lens < Propagator
    properties
        r
    end
    methods
        function obj = Lens(L, lambda, dx)
            obj@Propagator(L, lambda, dx);
            obj.r = sqrt(obj.X.^2 + obj.Y.^2);
        end
        
        function phase0 = phaseprofile(obj, coeffs, antenna_r)
            phase0 = 0;
            r1 = obj.r/antenna_r;
            for i = 1:length(coeffs)
                phase0 = phase0 + coeffs(i)*(r1.^(2*i));
            end
            phase0(obj.r > antenna_r) = 0;
        end
            
        function lens = makephaselens(obj, coeffs, antenna_r, padding)
            phase0 = obj.phaseprofile(coeffs, antenna_r);
            lens = exp(1i*phase0);
            if padding == 0
                lens(obj.r > antenna_r) = 0;
            end
        end
        
        function u1 = lenspropagate(obj, u0, lens, z1, z2)
            b4lens = obj.prop(u0, z1);
            afterlens = b4lens.*lens;
            u1 =  obj.prop(afterlens, z2);
        end
        
         function lens = makecplens(obj, z1, z2, antenna_r, padding)
            R1 = sqrt(obj.X.^2+obj.Y.^2+z1^2);
            R2 = sqrt(obj.X.^2+obj.Y.^2+z2^2);
            lens = exp(1i*(2*pi/obj.lambda)*(R1+R2));
            if padding == 0
                lens(obj.r > antenna_r) = 0;
            elseif padding == 1
                lens(obj.r > antenna_r) = 1;
            end
        end
        
        
    end
end
