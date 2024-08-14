function cpush = makePushCoil
%{
The design of the push coil called for windings that had a variable number of radial layers
(Nrad). In addition to the constraint imposed by the MOT beams, as can be seen in Figure 5.7a,
the front of the push coil must fit in between the horizontal transport mounting plates. The
simplest way was to make a push coil with variable Nrad for different layers. In total, the push
coil had 14 axial layers (Nax), and it was resolved that the easiest way to construct this was to
wind 7 subcoils each with two axial layers. Originally we tried winding the subcoils so that one
layer of the subcoil had a different Nrad than the other layer (ex. subcoil 5 having Nrad = 24, 25).
However, when winding the coils, we found that the coil would not hold its position well. In the
end we resolved to wrap subcoils with fixed Nrad. The dimensions of the subcoils are shown in
Figure 5.5.
%}



% Push coil OD corresponds to the last segment, and P inner separation corresponds to the distance from the first segment to
% the MOT center.

%% Coil 1 ("Push")
cpush = struct;
cpush.Name = 'Push';
cpush.Name2 = 'Push';
cpush.ID = 15*1e-3;
cpush.OD = NaN;
cpush.Nrad = [30 30 30 30 28 28 25 24 11 11 11 11 11 11];
cpush.Nrad = flip([30 30 30 30 28 28 25 24 11 11 11 11 11 11]);

cpush.Nax = length(cpush.Nrad);

cpush.WireDim = [1.02 2.29]*1e-3;
cpush.Position = -48*1e-3; % separation from front of push to the MOT center

% cpush.Position = -40*1e-3; % separation from front of push to the MOT center

cpush.InnerSeparation = NaN;
%% Create Coil Array

w = cpush.WireDim(1);
h = cpush.WireDim(2);

% inner radius
R0 = (cpush.ID/2);

Call = zeros(sum(cpush.Nrad),4);
ii = 1;
for ax = 1:length(cpush.Nrad)
    X0 = cpush.Position - (ax-0.5)*h;
    for rr = 1:cpush.Nrad(ax)
        C1 = [0 X0 R0+(rr-0.5)*w 1];
        Call(ii,:) = C1;ii=ii+1;
    end
end
cpush.Coil = Call;

end

