function [Call] = makeVerticalCoil(opt)
% From a structure array which describes the coils from a human
% perspective, this code creates a vector array which model each coil as a
% set of simple current loops. 
%
% Each coil loop is just
%
% [zposition radius current]
%
% This code initializes the current to be 1A, but that is modified later.

% Coil radial thickness
thickness = (opt.OD - opt.ID)/2;

Nr = opt.Nrad;
Nax = opt.Nax;

% This is the individual wire thickness, this should match the first
% dimension of the wire (it is rectangular).
w = thickness/Nr;

w = opt.WireDim(1);
h = opt.WireDim(2);

% Radius of inner-most coil
R0 = (opt.ID/2) + w/2;

% Z position of bottom layer
Z0 = opt.zbot + h/2;

% X Location
X0 = 0;

Call = zeros(Nax*Nr,4);
ii = 1;
for cc = 0:(Nax-1)
    for rr = 0:(Nr-1)
        C1 = [X0 +(Z0+cc*h) R0+rr*w 1];
        Call(ii,:) = C1;ii=ii+1;
    end
end


end

