function coils = makeHorizontalCoils

%% Coil 2 ("MOT")
c2 = struct;
c2.Name = 'MOT';
c2.Name2 = 'MOT';
c2.ID = 45.7*1e-3;
c2.OD = 96*1e-3;
c2.Nrad = 33;
c2.Nax = 2;
c2.WireDim = [0.73 2.63]*1e-3;
c2.InnerSeparation = 50*1e-3;
c2.Position = [0 0]*1e-3;

%% Coil 3 ("Horiz 1")
c3 = struct;
c3.Name = 'Coil 3';
c3.Name2 = 'Horiz 1';
c3.ID = 24.8*1e-3;
c3.OD = 61*1e-3;
c3.Nrad = 16;
c3.Nax = 2;
c3.WireDim = [1.13 2.63]*1e-3;
c3.InnerSeparation = 63.48*1e-3;
c3.Position = [53.3 0]*1e-3;

%% Coil 4 ("Horiz 2")
c4 = struct;
c4.Name = 'Coil 4';
c4.Name2 = 'Horiz 2';
c4.ID = 24.8*1e-3;
c4.OD = 61*1e-3;
c4.Nrad = 16;
c4.Nax = 2;
c4.WireDim = [1.13 2.63]*1e-3;
c4.InnerSeparation = 51.48*1e-3;
c4.Position = [84.8 0]*1e-3;

%% Coil 5 ("Horiz 3")
c5 = struct;
c5.Name = 'Coil 3';
c5.Name2 = 'Horiz 3';
c5.ID = 24.8*1e-3;
c5.OD = 61*1e-3;
c5.Nrad = 16;
c5.Nax = 2;
c5.WireDim = [1.13 2.63]*1e-3;
c5.InnerSeparation = 63.48*1e-3;
c5.Position = [116.3 0]*1e-3;

%% Coil 6 ("Horiz 4")
c6 = struct;
c6.Name = 'Coil 6';
c6.Name2 = 'Horiz 4';
c6.ID = 24.8*1e-3;
c6.OD = 61*1e-3;
c6.Nrad = 16;
c6.Nax = 2;
c6.WireDim = [1.13 2.63]*1e-3;
c6.InnerSeparation = 51.48*1e-3;
c6.Position = [147.8 0]*1e-3;

%% Coil 7 ("Horiz 5")
c7 = struct;
c7.Name = 'Coil 7';
c7.Name2 = 'Horiz 5';
c7.ID = 24.8*1e-3;
c7.OD = 61*1e-3;
c7.Nrad = 16;
c7.Nax = 2;
c7.WireDim = [1.13 2.63]*1e-3;
c7.InnerSeparation = 63.48*1e-3;
c7.Position = [179.3 0]*1e-3;

%% Coil 8 ("Horiz 6")
c7 = struct;
c7.Name = 'Coil 8';
c7.Name2 = 'Horiz 6';
c7.ID = 24.8*1e-3;
c7.OD = 61*1e-3;
c7.Nrad = 16;
c7.Nax = 2;
c7.WireDim = [1.13 2.63]*1e-3;
c7.InnerSeparation = 51.48*1e-3;
c7.Position = [210.8 0]*1e-3;

%% Coil 9 ("Horiz 7")
c8 = struct;
c8.Name = 'Coil 8';
c8.Name2 = 'Horiz 7';
c8.ID = 24.8*1e-3;
c8.OD = 61*1e-3;
c8.Nrad = 16;
c8.Nax = 2;
c8.WireDim = [1.13 2.63]*1e-3;
c8.InnerSeparation = 63.48*1e-3;
c8.Position = [242.3 0]*1e-3;

%% Coil 9 ("Horiz 8")
c9 = struct;
c9.Name = 'Coil 9';
c9.Name2 = 'Horiz 8';
c9.ID = 24.8*1e-3;
c9.OD = 61*1e-3;
c9.Nrad = 16;
c9.Nax = 2;
c9.WireDim = [1.13 2.63]*1e-3;
c9.InnerSeparation = 51.48*1e-3;
c9.Position = [273.8 0]*1e-3;

%% Coil 10 ("Horiz 9")
c10 = struct;
c10.Name = 'Coil 10';
c10.Name2 = 'Horiz 9';
c10.ID = 24.8*1e-3;
c10.OD = 61*1e-3;
c10.Nrad = 16;
c10.Nax = 2;
c10.WireDim = [1.13 2.63]*1e-3;
c10.InnerSeparation = 63.48*1e-3;
c10.Position = [305.3 0]*1e-3;

%% Assemble Coils

coils = [c3 c4 c5 c6 c7 c8 c9 c10];

for kk=1:length(coils)
coils(kk).Coil = makeCoilPair(coils(kk));
end

end


