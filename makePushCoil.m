function coil = makePushCoil
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
cpush.OD = 75.6*1e-3; % NOT ACTUALLY THE OD
cpush.Nrad = [30 30 30 30 28 28 25 24 11 11 11 11 11 11];

cpush.WireDim = [1.02 2.29]*1e-3;
cpush.InnerSeparation = 42.5*1e-3; % NOT ACTUALLY THE SEPARATION, UNDEFIND
cpush.Position = [NaN NaN]*1e-3;

end

