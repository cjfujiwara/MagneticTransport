%% Specify coils parameters
% Coil radius
R=1;

% First set
C1=[0 -R/2 R];
C2=[0 R/2 R];

% Second set
C3=[R/2 -R/2 R];
C4=[R/2 R/2 R];

% Third set
C5=[R -R/2 R];
C6=[R R/2 R];

% Group the coil specifications together
Cmat=[C1;C2;C3;C4;C5;C6];
%% Specify the X and Y mesh

N=1000;
inds=[0 1]+N;
xVec=linspace(-1.25*R,1.25*R,N); % inspect over many x's
y=[-1 1]*R*1E-2;              % inspect over only the center y
z=[-1 1]*R*1E-2;               % inspect over only the center z


dZ=1E-2*R;
dR=1E-2*R;
for kk=1:length(xVec)
    x=xVec(kk);
        
    % Calculate the magnitude of the field at (x,y=0,z=0);
    [Bx1,By1,Bz1]=fieldCoil3(x,0,0,Cmat(ind,:));
    [Bx2,By2,Bz2]=fieldCoil3(x,0,0,Cmat(ind+1,:));
    
    % Superpose the fields    
    Bz=Bz1-Bz2;Bx=Bx1-Bx2;By=By1-By2;
        
    % Calculate the magnitude at the center point
    B=sqrt(Bz.^2+By.^2+Bz.^2);
    
    % Calculate the z field gradient
    [~,~,Bz1]=fieldCoil3(x,0,dZ/2,Cmat(ind,:));
    [~,~,Bz2]=fieldCoil3(x,0,dZ/2,Cmat(ind+1,:));    
    [~,~,Bz3]=fieldCoil3(x,0,-dZ/2,Cmat(ind,:));
    [~,~,Bz4]=fieldCoil3(x,0,-dZ/2,Cmat(ind+1,:));    
    Gz=((Bz3-Bz4)-(Bz1-Bz2))/dZ;
    
    % Calculate the y field gradient
    [~,By1,~]=fieldCoil3(x,dR/2,0,Cmat(ind,:));
    [~,By2,~]=fieldCoil3(x,dR/2,0,Cmat(ind+1,:));    
    [~,By3,~]=fieldCoil3(x,-dR/2,0,Cmat(ind,:));
    [~,By4,~]=fieldCoil3(x,-dR/2,0,Cmat(ind+1,:));    
    Gy=((By3-By4)-(By1-By2))/dR;
    
    % Calculate the y field gradient
    [~,By1,~]=fieldCoil3(x,dR/2,0,Cmat(ind,:));
    [~,By2,~]=fieldCoil3(x,dR/2,0,Cmat(ind+1,:));    
    [~,By3,~]=fieldCoil3(x,-dR/2,0,Cmat(ind,:));
    [~,By4,~]=fieldCoil3(x,-dR/2,0,Cmat(ind+1,:));    
    Gy=((By3-By4)-(By1-By2))/dR;
    
%     Gz=Bz(2,:,
    
    % Calculate the gradients

    
    


end

