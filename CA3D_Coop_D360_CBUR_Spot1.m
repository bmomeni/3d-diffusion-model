%% Diffusion in 3D
% 3D diffusion for transfer of nutrients
% nonhomogeneous media
% space occupied by dead cells
% continuous uptake
% v2: removed agar from cell domain to occupy less memory
% v3: decoupled lateral size of cell domain from that of diffusion to save
% search space when dealing with cells
% v4: checks all neighboring cells (instead of eight attempts) before
% budding up; cells uniformly distributed over the surface; cells drop down
% if there is empty space underneath
% PBC: periodic boundary condition
% LRA: live lys-requiring cells release adenine
% UU1: Unrestricted uptake (not limited to available nutrient within a
% single spatial unit; lowered diffusion time-step for better accuracy

clear
rand('twister',5489)

savename = 'CA3DCOOP_CBUR072_Spot_SC10_1335N500_1340N500_t120_rg5_N160Nz60Ns1000';

%% Parameters
SSF = 5e-4; % spatial scaling factor; grid size in cm
r1 = 0.51; % 1335 reproduction rate, per hour
r2 = 0.44; % 1340 reproduction rate, per hour
d1 = 0.0021; % death rate, per hour
d2 = 0.015; % death rate, per hour
alphaA = 3.1; % reproduction nutrient consumption factor; fmole per cell
alphaL = 5.4; % reproduction nutrient consumption factor; fmole per cell
gammaA = 0.26; % hyp release rate; fmole per cell per hour
% gammaL = 0.4; % death nutrient release factor; fmole per cell
HgL = 1e6*[0 0.58 0.65 0.73 0.80 1.07 1000]; % Hyp levels (in fmole/ml) at which release rate was measured
gammaLH = [0.52 0.83 0.78 0.53 0.38 0.08 0.08]; % measured Lys release rate (fmole/cell/hr) at given [Hyp]
nL = 3.2; % Moser's power for Lys
nA = 1.5; % Moser's power for Hyp
KA = 1.3e6; % Monod constant; fmole/ml
KL = 2.1e6; % Monod constant; fmole/ml

vmL = alphaL*r1; % maximum uptake rate, fmole/hour
vmA = alphaA*r2; % maximum uptake rate, fmole/hour

KmmL = KL; % Michaelis-Menten constant
KmmA = KA; % Michaelis-Menten constant

D0 = 360; % diffusion constant in agarose/agar, microns2/s
D1 = 360; % diffusion constant in cells, microns2/s

%% Definition of solution domain
%% Cell simulation domain
Nc = 160;
Nz = 60;
nc = round(Nc/2);
% [xx,yy,zz] = meshgrid(1:Nc,1:Nc,1:Nz);

%% Diffusion simulation domain
SC = 10; % ratio of diffusion grid to cell grid
Nd = 16; % domain size
Nsh = round(0.5*(Nd - floor(Nc/SC)));
Ns = 100; % thickness of agarose substrate
Ndz = floor(Ns+Nz/SC); % height
nd = round(Nd/2);
ndz = round(Ndz/2);

SA = zeros(Nd,Nd,Ndz);
SL = zeros(Nd,Nd,Ndz);
SL(:,:,1:Ns) = 0.7e6; % initial lys conc. in agar, fmole/ml
SAde = zeros(Nc,Nc,Nz);
SLde = zeros(Nc,Nc,Nz);

De = zeros(Nd+2,Nd+2,Ndz+2);

D = zeros(Nd,Nd,Ndz);
D(:,:,1:Ns) = D0;

% Length scale
g = 5; % grid size for cells, microns
d = g*SC; % grid size for diffusion, microns
dVc = g^3 * 1e-12; % volume of each grid cube for cells, ml
dV = d^3 * 1e-12; % volume of each grid cube for diffusion, ml

rg0 = 5; % neghborhood radius for expansion to sides

%% Initial cell distrbution
CSD = 500; % initial cell surface density per type, 1/mm2
Ncs = Nc/4;
Nc1 = round(CSD*g*Ncs*g*Ncs/1e6); % initial number of type 1 cells, (R)
Nc2 = round(CSD*g*Ncs*g*Ncs/1e6); % initial number of type 2 cells, (Y)
X = zeros(Nc,Nc,Nz); % Cells
U = zeros(Nc,Nc,Nz); % Nutrients taken up by cells
NB = zeros(Nc,Nc,Nz); % Number of budding events at each location
nn = 0;
while nn < Nc1
    i = floor(nc+Ncs*(rand(1)-0.5)+1);
    j = floor(nc+Ncs*(rand(1)-0.5)+1);
    k = 1;
    if X(i,j,k) == 0
        nn = nn+1;
        X(i,j,k) = 1;
        U(i,j,k) = rand(1)*alphaL; %to account for the asynchronous nature of cells - how much nutrient accumulated up to t=0
    end
end
nn = 0;
while nn < Nc2
    i = floor(nc+Ncs*(rand(1)-0.5)+1);
    j = floor(nc+Ncs*(rand(1)-0.5)+1);
    k = 1;
    if X(i,j,k) == 0
        nn = nn+1;
        X(i,j,k) = 2;
        U(i,j,k) = rand(1)*alphaA; %to account for the asynchronous nature of cells - how much nutrient accumulated up to t=0
    end
end

%% Cell-growth time-course
tau0 = 0; % in hours
tauf = 120; % in hours
dtau = 0.1; % in hours
ts = 1;

r1e = r1*dtau; % probability of reproduction in a time step
d1e = d1*dtau; % probability of death in a time step
r2e = r2*dtau; % probability of reproduction in a time step
d2e = d2*dtau; % probability of death in a time step

taurng = tau0:dtau:tauf;
% tu = 0.1*KmmL*dV/vmL; % uptake update time-scale
ctu = 0; % counter for uptake timing

X1lm = zeros(size(taurng)); % number of live cells of type 1
X2lm = zeros(size(taurng)); % number of live cells of type 1
X1m = zeros(size(taurng)); % number of live cells of type 1
X2m = zeros(size(taurng)); % number of live cells of type 1
UAcc1 = zeros(size(taurng)); % Accumulated nutrients in live cells
UAcc2 = zeros(size(taurng)); % Accumulated nutrients in live cells
UW1 = zeros(size(taurng)); % Wasted nutrients in dead cells
UW2 = zeros(size(taurng)); % Wasted nutrients in dead cells
SLAccA = zeros(size(taurng)); % Total nutrients in the agar region
SLAccC = zeros(size(taurng)); % Total nutrients in the cell region
SAAccA = zeros(size(taurng)); % Total nutrients in the agar region
SAAccC = zeros(size(taurng)); % Total nutrients in the cell region

ct = 0;
cS = 0;

Q = [1 0; -1 0; 0 1; 0 -1; 1 -1; -1 1; 1 1; -1 -1]; % locations of eight neighbor grids
Qc = [1 0; -1 0; 0 1; 0 -1]; % locations of four neighbor grids

for tau = taurng
    ct = ct+1;
    tic
    
    for z = 1:Nz-1
        zd = floor((z-1)/SC)+Ns+1;
        [I,J] = find((X(:,:,z)==1)+(X(:,:,z)==2));
        Ncc = length(I);
        [SS,OS] = sort(rand(1,Ncc));
        I = I(OS);
        J = J(OS);
        for cc = 1:Ncc
            
            Id = floor((I(cc)-1)/SC)+Nsh+1;
            Jd = floor((J(cc)-1)/SC)+Nsh+1;
            xd = I(cc); % location of cell along x
            yd = J(cc); % location of cell along y
            
            live = 1;
            % Is the cell dead?
            cd = rand(1);
            if (X(xd,yd,z) == 1)&&(cd < d1e)&&(sum(sum(sum(X==1)))>1)
                X(xd,yd,z) = 0.5;
                live = 0;
            end
            if (X(xd,yd,z) == 2)&&(cd < d2e)&&(sum(sum(sum(X==2)))>1)
                X(xd,yd,z) = 1.5;
                live = 0;
            end
            
            % Cell division and rearrangement
            if (X(xd,yd,z) == 1)
                alpha = alphaL;
                celltype = 1;
            end
            if (X(xd,yd,z) == 2)
                alpha = alphaA;
                celltype = 2;
            end
            
            if (live==1)&&(U(xd,yd,z) >= alpha)
                Budded = 0;
                Natt = 0;
                U(xd,yd,z) = U(xd,yd,z)-alpha;
                [qm1,qi1] = sort(rand(1,4)); % random index for neighboring grids
                [qm2,qi2] = sort(rand(1,4)); % random index for neighboring grids
                qi = [qi1 4+qi2];
                while (Budded==0)&&(Natt<8) % try immediate neighboring grids
                    Natt = Natt+1;
                    xb = I(cc) + Q(qi(Natt),1);
                    if xb > Nc, xb = xb-Nc; end
                    if xb < 1, xb = xb+Nc; end
                    yb = J(cc) + Q(qi(Natt),2);
                    if yb > Nc, yb = yb-Nc; end
                    if yb < 1, yb = yb+Nc; end
                    
                    if X(xb,yb,z)==0 % if available space, bud into it
                        zb = z;
                        while (zb>=2)&&(X(xb,yb,zb-1)==0)
                            zb = zb-1;
                        end
                        Budded = 1;
                        X(xb,yb,zb) = X(xd,yd,z);
                        NB(xd,yd,z) = NB(xd,yd,z)+1;
                        U(xb,yb,zb) = 0;
                    end
                end
                if (Budded==0) % try extended neighborhood in the same plane
                    rg = rg0;
                    xdp = xd;
                    ydp = yd;
                    if xd+rg > Nc
                        Xp = [X(:,:,z); X(1:xd+rg-Nc,:,z)];  % extend along +x; size(Xp)=xd+rg
                    else
                        if xd-rg < 1
                            Xp = [X(Nc+xd-rg:Nc,:,z); X(:,:,z)]; % extend along -x; size(Xp)=Nc-xd+rg+1
                            xdp = rg+1;
                        else
                            Xp = X(:,:,z);
                        end
                    end
                    if yd+rg > Nc
                        Xp = [Xp, Xp(:,1:yd+rg-Nc)];
                    else
                        if yd-rg < 1
                            Xp = [Xp(:,Nc+yd-rg:Nc), Xp];
                            ydp = rg+1;
                        end
                    end
                    Ng = Xp(xdp-rg:xdp+rg,ydp-rg:ydp+rg);
                    if sum(sum(Ng>0.25))<(2*rg+1)^2
                        [Ie,Je] = find(Ng==0);
                        [SSe,OSe] = min(sqrt((Ie-rg-1).^2+(Je-rg-1).^2));
                        if SSe <= rg
                            xe = Ie(OSe);
                            ye = Je(OSe);
                            TP = TracePath(rg+1,rg+1,xe,ye);
                            NT = size(TP,1);
                            TP(:,1) = TP(:,1) + xd-rg-1;
                            TP(:,1) = TP(:,1) + Nc*((TP(:,1)<1)-(TP(:,1)>Nc));
                            TP(:,2) = TP(:,2) + yd-rg-1;
                            TP(:,2) = TP(:,2) + Nc*((TP(:,2)<1)-(TP(:,2)>Nc));
                            zb = z;
                            while (zb>=2)&&(X(TP(NT,1),TP(NT,2),zb-1)==0)
                                zb = zb-1;
                            end
                            X(TP(NT,1),TP(NT,2),zb) = X(TP(NT-1,1),TP(NT-1,2),z);
                            U(TP(NT,1),TP(NT,2),zb) = U(TP(NT-1,1),TP(NT-1,2),z);
                            for ctr = NT-1:-1:2
                                X(TP(ctr,1),TP(ctr,2),z) = X(TP(ctr-1,1),TP(ctr-1,2),z);
                                U(TP(ctr,1),TP(ctr,2),z) = U(TP(ctr-1,1),TP(ctr-1,2),z);
                            end
                            X(TP(1,1),TP(1,2),z) = X(xd,yd,z);
                            U(TP(1,1),TP(1,2),z) = 0;
                            Budded = 1;
                        end
                    end
                end
                if (Budded == 0) % bud to top or sides and push cells above
                    cq = rand(1);
                    if cq < 0.7
                        X(xd,yd,z+1:Nz) = X(xd,yd,z:Nz-1);
                        U(xd,yd,z+1:Nz) = U(xd,yd,z:Nz-1);
                    else
                        [qm3,qi3] = sort(rand(1,4)); % random index for neighboring grids
                        xb = I(cc) + Qc(qi3(1),1);
                        if xb > Nc, xb = xb-Nc; end
                        if xb < 1, xb = xb+Nc; end
                        yb = J(cc) + Qc(qi3(1),2);
                        if yb > Nc, yb = yb-Nc; end
                        if yb < 1, yb = yb+Nc; end
                        X(xb,yb,z+1:Nz) = X(xb,yb,z:Nz-1);
                        U(xb,yb,z+1:Nz) = U(xb,yb,z:Nz-1);
                        X(xb,yb,z) = celltype;
                        U(xb,yb,z) = 0;
                        xd = xb;
                        yd = yb;
                    end
                    zt = Nz;
                    while (zt>z)&&(X(xd,yd,zt)==0)
                        zt = zt-1;
                    end
                    Xe = [X(Nc,:,zt-1); X(:,:,zt-1); X(1,:,zt-1)];
                    Xe = [Xe(:,Nc), Xe, Xe(:,1)];
                    if sum(sum(Xe(xd:xd+2,yd:yd+2)>0.1))<=6
                        [xe,ye] = find(Xe(xd:xd+2,yd:yd+2)==0);
                        et = floor(rand(1)*length(xe))+1;
                        xet = xd + xe(et) - 2;
                        xet = xet + Nc*((xet<1)-(xet>Nc));
                        yet = yd + ye(et) - 2;
                        yet = yet + Nc*((yet<1)-(yet>Nc));
                        zn = zt-1;
                        while (zn>=2)&&(X(xet,yet,zn-1)==0)
                            zn = zn-1;
                        end
                        X(xet,yet,zn) = X(xd,yd,zt);
                        U(xet,yet,zn) = U(xd,yd,zt);
                        X(xd,yd,zt) = 0;
                        U(xd,yd,zt) = 0;
                    end
                end
            end
            
        end
        
        X1lm(z,ct) = sum(sum(X(:,:,z)==1));
        X2lm(z,ct) = sum(sum(X(:,:,z)==2));
        X1m(z,ct) = sum(sum((X(:,:,z)>0.25).*(X(:,:,z)<1.25)));
        X2m(z,ct) = sum(sum((X(:,:,z)>1.25).*(X(:,:,z)<2.25)));
    end
    UAcc1(ct) = sum(sum(sum(U.*(X==1))));
    UAcc2(ct) = sum(sum(sum(U.*(X==2))));
    UW1(ct) = sum(sum(sum(U.*(X==0.5))));
    UW2(ct) = sum(sum(sum(U.*(X==1.5))));
    SLAccC(ct) = sum(sum(sum(SL(:,:,Ns+1:Ndz))));
    SLAccA(ct) = sum(sum(sum(SL(:,:,1:Ns))));
    SAAccC(ct) = sum(sum(sum(SA(:,:,Ns+1:Ndz))));
    SAAccA(ct) = sum(sum(sum(SA(:,:,1:Ns))));
    disp([tau  sum(X1lm(:,ct)) sum(X2lm(:,ct)) 1e-6*mean(mean(mean(SL))) 1e-6*mean(mean(mean(SA)))])
    
    %% Diffusion loop
    % Update diffusion constant
    D(:,:,Ns+1:Ndz) = zeros(Nd,Nd,Ndz-Ns);
    for i1 = 1:SC
        for i2 = 1:SC
            for i3 = 1:SC
                D(:,:,Ns+1:Ndz) = D(:,:,Ns+1:Ndz) + D1/SC^3*(X(i1:SC:Nc,i2:SC:Nc,i3:SC:Nz)>0.1);
            end
        end
    end
    De(2:Nd+1,2:Nd+1,2:Ndz+1) = D(1:Nd,1:Nd,1:Ndz);
    De(1,2:Nd+1,2:Ndz+1) = D(Nd,1:Nd,1:Ndz);
    De(Nd+2,2:Nd+1,2:Ndz+1) = D(1,1:Nd,1:Ndz);
    De(2:Nd+1,1,2:Ndz+1) = D(1:Nd,Nd,1:Ndz);
    De(2:Nd+1,Nd+2,2:Ndz+1) = D(1:Nd,1,1:Ndz);
    De(2:Nd+1,2:Nd+1,1) = D(1:Nd,1:Nd,1);
    De(2:Nd+1,2:Nd+1,Ndz+2) = D(1:Nd,1:Nd,Ndz);
    
    t0 = 0; % in seconds
    dt = 0.072*d^2/D0; % in seconds
    tf = dtau*3600; % in seconds
    trng = t0:dt:tf; % in seconds
    
    for t = trng
        
        SLe(2:Nd+1,2:Nd+1,2:Ndz+1) = SL(1:Nd,1:Nd,1:Ndz);
        SLe(1,2:Nd+1,2:Ndz+1) = SL(Nd,1:Nd,1:Ndz);
        SLe(Nd+2,2:Nd+1,2:Ndz+1) = SL(1,1:Nd,1:Ndz);
        SLe(2:Nd+1,1,2:Ndz+1) = SL(1:Nd,Nd,1:Ndz);
        SLe(2:Nd+1,Nd+2,2:Ndz+1) = SL(1:Nd,1,1:Ndz);
        SLe(2:Nd+1,2:Nd+1,1) = SL(1:Nd,1:Nd,1);
        SLe(2:Nd+1,2:Nd+1,Ndz+2) = SL(1:Nd,1:Nd,Ndz);
        
        dSL = dt/d^2 * De(2:Nd+1,2:Nd+1,2:Ndz+1).*(SLe(1:Nd,2:Nd+1,2:Ndz+1) + SLe(3:Nd+2,2:Nd+1,2:Ndz+1) + SLe(2:Nd+1,1:Nd,2:Ndz+1) + SLe(2:Nd+1,3:Nd+2,2:Ndz+1) + SLe(2:Nd+1,2:Nd+1,1:Ndz) + SLe(2:Nd+1,2:Nd+1,3:Ndz+2) - 6*SLe(2:Nd+1,2:Nd+1,2:Ndz+1)) + dt/d^2 * ((De(3:Nd+2,2:Nd+1,2:Ndz+1)-De(2:Nd+1,2:Nd+1,2:Ndz+1)).*(SLe(3:Nd+2,2:Nd+1,2:Ndz+1)-SLe(2:Nd+1,2:Nd+1,2:Ndz+1)) + (De(2:Nd+1,3:Nd+2,2:Ndz+1)-De(2:Nd+1,2:Nd+1,2:Ndz+1)).*(SLe(2:Nd+1,3:Nd+2,2:Ndz+1)-SLe(2:Nd+1,2:Nd+1,2:Ndz+1)) + (De(2:Nd+1,2:Nd+1,3:Ndz+2)-De(2:Nd+1,2:Nd+1,2:Ndz+1)).*(SLe(2:Nd+1,2:Nd+1,3:Ndz+2)-SLe(2:Nd+1,2:Nd+1,2:Ndz+1)));
        SL = SL + dSL;
        
        SAe(2:Nd+1,2:Nd+1,2:Ndz+1) = SA(1:Nd,1:Nd,1:Ndz);
        SAe(1,2:Nd+1,2:Ndz+1) = SA(Nd,1:Nd,1:Ndz);
        SAe(Nd+2,2:Nd+1,2:Ndz+1) = SA(1,1:Nd,1:Ndz);
        SAe(2:Nd+1,1,2:Ndz+1) = SA(1:Nd,Nd,1:Ndz);
        SAe(2:Nd+1,Nd+2,2:Ndz+1) = SA(1:Nd,1,1:Ndz);
        SAe(2:Nd+1,2:Nd+1,1) = SA(1:Nd,1:Nd,1);
        SAe(2:Nd+1,2:Nd+1,Ndz+2) = SA(1:Nd,1:Nd,Ndz);
        
        dSA = dt/d^2 * De(2:Nd+1,2:Nd+1,2:Ndz+1).*(SAe(1:Nd,2:Nd+1,2:Ndz+1) + SAe(3:Nd+2,2:Nd+1,2:Ndz+1) + SAe(2:Nd+1,1:Nd,2:Ndz+1) + SAe(2:Nd+1,3:Nd+2,2:Ndz+1) + SAe(2:Nd+1,2:Nd+1,1:Ndz) + SAe(2:Nd+1,2:Nd+1,3:Ndz+2) - 6*SAe(2:Nd+1,2:Nd+1,2:Ndz+1)) + dt/d^2 * ((De(3:Nd+2,2:Nd+1,2:Ndz+1)-De(2:Nd+1,2:Nd+1,2:Ndz+1)).*(SAe(3:Nd+2,2:Nd+1,2:Ndz+1)-SAe(2:Nd+1,2:Nd+1,2:Ndz+1)) + (De(2:Nd+1,3:Nd+2,2:Ndz+1)-De(2:Nd+1,2:Nd+1,2:Ndz+1)).*(SAe(2:Nd+1,3:Nd+2,2:Ndz+1)-SAe(2:Nd+1,2:Nd+1,2:Ndz+1)) + (De(2:Nd+1,2:Nd+1,3:Ndz+2)-De(2:Nd+1,2:Nd+1,2:Ndz+1)).*(SAe(2:Nd+1,2:Nd+1,3:Ndz+2)-SAe(2:Nd+1,2:Nd+1,2:Ndz+1)));
        SA = SA + dSA;
        
        %% Uptake
        dtu = dt/3600; % in hours
        gammaL = interp1(HgL,gammaLH,SA(:,:,Ns+1:Ndz));
        for i1 = 1:SC
            for i2 = 1:SC
                for i3 = 1:SC
                    SLde(i1:SC:Nc,i2:SC:Nc,i3:SC:Nz) = SL(:,:,Ns+1:Ndz) + dtu/dVc*gammaL.*(X(i1:SC:Nc,i2:SC:Nc,i3:SC:Nz)==2).*(1+1/alphaA*U(i1:SC:Nc,i2:SC:Nc,i3:SC:Nz));
                    SAde(i1:SC:Nc,i2:SC:Nc,i3:SC:Nz) = SA(:,:,Ns+1:Ndz) + dtu/dVc*gammaA*(X(i1:SC:Nc,i2:SC:Nc,i3:SC:Nz)==1).*(1+1/alphaL*U(i1:SC:Nc,i2:SC:Nc,i3:SC:Nz));
                end
            end
        end
        UcL = (X==1).*(1+1/alphaL*U).*(dtu*vmL*SLde.^nL./(SLde.^nL+KmmL^nL).*(SLde>0));
        UcA = (X==2).*(1+1/alphaA*U).*(dtu*vmA*SAde.^nA./(SAde.^nA+KmmA^nA).*(SAde>0));
        U = U+UcL+UcA;
        SLde = SLde-1/dVc*UcL;
        SAde = SAde-1/dVc*UcA;
        SL(:,:,Ns+1:Ndz) = zeros(Nd,Nd,floor(Nz/SC));
        SA(:,:,Ns+1:Ndz) = zeros(Nd,Nd,floor(Nz/SC));
        for ii = 1:SC
            for jj = 1:SC
                for kk = 1:SC
                    SL(:,:,Ns+1:Ndz) = SL(:,:,Ns+1:Ndz) + 1/SC^3*SLde(ii:SC:Nc,jj:SC:Nc,kk:SC:Nz);
                    SA(:,:,Ns+1:Ndz) = SA(:,:,Ns+1:Ndz) + 1/SC^3*SAde(ii:SC:Nc,jj:SC:Nc,kk:SC:Nz);
                end
            end
        end
        
    end
    if mod(tau+0.03,ts)<0.1
        cS = cS+1;
        save(strcat(savename,'_ts1_cs',num2str(cS)),'X','tau','SL','SA')
    end
    toc
    
end

save(strcat(savename,'.mat'))
