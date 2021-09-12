
classdef Propagator
     properties
        L
        lambda
        M
        dx
        x
        y
        X
        Y
     end
     
     methods
         function obj = Propagator(L, lambda, dx)
             % Constructor
             obj.L = L;
             obj.lambda = lambda;
             obj.dx = dx;
             obj.M = L/dx+1;
             obj.x = -L/2:obj.dx:L/2;
             obj.y = -L/2:obj.dx:L/2;
             [obj.X, obj.Y] = meshgrid(obj.x, obj.y);
         end
         
         function u0 = pso(obj, xangle, yangle, R1)
             % make point source
             xshift = R1*sind(xangle);
             xshiftblocks = fix(xshift/obj.dx);
             yshift = R1*sind(yangle);
             yshiftblocks = fix(yshift/obj.dx);
             u0 = zeros(obj.M);
             u0(ceil(obj.M/2)-yshiftblocks,ceil(obj.M/2)-xshiftblocks) = 1;
         end
         
         function u1 = prop(obj, u0, z)

           k = 2*pi/obj.lambda;

           fx = -1/(2*obj.dx):1/obj.L:1/(2*obj.dx);
           fy = -1/(2*obj.dx):1/obj.L:1/(2*obj.dx);
           [Fx, Fy] = meshgrid(fx,fy);
           H = exp(-1i*k*z*sqrt(1 - (obj.lambda*Fx).^2 - (obj.lambda*Fy).^2));
           H(Fx.^2 + Fy.^2 > 1/obj.lambda^2) = 0;
           U0 = fftshift(fft2(fftshift(u0)));
           U1 = H.*U0;
           u1 = ifftshift(ifft2(ifftshift(U1)));
         end
         
         function u0 = backprop(obj, u0, z)
           u0 = obj.prop(obj, u0, -z);
         end        
     end

end