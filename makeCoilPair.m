function Call = makeCoilPair(opt)

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

% Z Position of closet coil (assumed to have a pair)
Z0 = (opt.InnerSeparation/2) + h/2;

% Axis Locatio
X0 = opt.Position(1);

Call = zeros(Nax*Nr*2,4);
ii = 1;
for cc = 0:(Nax-1)
    for rr = 0:(Nr-1)
        C1 = [X0 +(Z0+cc*h) R0+rr*w +1];
        C2 = [X0 -(Z0+cc*h) R0+rr*w -1];
        Call(ii,:) = C1;ii=ii+1;
        Call(ii,:) = C2;ii=ii+1;
    end
end

end

