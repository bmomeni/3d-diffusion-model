%% Cellular Automata model for growth of yeast communities
%% Diffusion in 3D
% 3D diffusion for transfer of nutrients
% nonhomogeneous media
% space occupied by dead cells
% CU: continuous uptake
% v2: removed agar from cell domain to occupy less memory
% v3: decoupled lateral size of cell domain from that of diffusion to save
% search space when dealing with cells
% v4: checks all neighboring cells (instead of eight attempts) before
% budding up; cells uniformly distributed over the surface; cells drop down
% if there is empty space underneath
% v5: No Spatial: no limitation for cells confined among other cells,
% smaller diffusion constant and grid size for diffusion in cells
% v6: Shorter uptake timescale; track growth at different heights; random
% order in height when updating cells
% variation 2: uptake update in sync. with diffusion
% v7: expansion from sides with fixed radius
% FCU: fast uptake
% CI: confinement independent
% NoES2: Expansion to sides only in the first layer on the surface
% BT: bud to top, when confined from sides
% EB: Boundary condition (continuous flow) explicitly applied at the agar-community interface
% DR: Delayed release of lys
% EBP: Periodic boundary condition on the sides (in addition to EB)
% Conf2: each population physiaclly confined to a x-y region
% NoES2: no expansion to the sides, except for the very first cell layer (which rearranges according to ES2)
% GD: glucose-dependent growth
% CBUR: Both release and uptake are modulated by cell biomass

clear
rseed = 2417;
rand('twister',rseed)
savename = 'CA3DCOOP_D20_CBUR05_Spot_ECOPHYS_c5g15d60_1335N500_1340N500_t120_rg5Nc48Na16_Ndz30Naz80';

%% Parameters
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
D1 = 20; % diffusion constant in cells, microns2/s

%% Definition of solution domain
% Length scale
c = 5; % grid size for cells, microns
SC = 3; % ratio of diffusion to cell grids
SD = 4; % ratio of diffusion grid (agar to cells)
g = SC*c; % grid size for diffusion in cells, microns
d = SD*g; % grid size for diffusion in agar, microns
dV = c^3 * 1e-12; % volume of each grid cube for cells, ml
dVc = g^3 * 1e-12; % volume of each grid cube for diffusion in cells, ml
dVa = d^3 * 1e-12; % volume of each grid cube for diffusion in agar, ml
h = 0.5*(g+d); % grid size at the interface between the two grids

rg0 = 5; % neghborhood radius for expansion to sides

%% Diffusion simulation domain
Nax = 16; % agar domain size along x
Nay = 16; % agar domain size along y
Naz = 80; % agar domain height
Ndx = Nax*SD; % cell domain size for diffusion
Ndy = Nay*SD; % cell domain size for diffusion
Ndz = 30; % cell domain height for diffusion

%% Cell simulation domain
Ncx = Ndx*SC; % cell domain size for cells
Ncy = Ndy*SC; % cell domain size for cells
Ncz = Ndz*SC; % cell domain height for cells

Ncsx = round(Ncx/4); % initial spot for cells
Ncsy = round(Ncy/4); % initial spot for cells

SLa = 0.7e6*ones(Nax,Nay,Naz); %zeros(Nax,Nay,Naz); % fmole/ml initial concentration
SLd = zeros(Ndx,Ndy,Ndz);
SLde = zeros(Ncx,Ncy,Ncz);
SAa = zeros(Nax,Nay,Naz);
SAd = zeros(Ndx,Ndy,Ndz);
SAde = zeros(Ncx,Ncy,Ncz);

Dd = zeros(Ndx,Ndy,Ndz);
SLb = zeros(Ndx,Ndy); % concentration of lys at the boundary of agar and community; 0.5 or 1 or 0 makes no difference
SAb = zeros(Ndx,Ndy); % concentration of ade at the boundary of agar and community; 0.5 or 1 or 0 makes no difference

[xxd,yyd] = meshgrid((SD+1)/2:Ndx+(SD-1)/2,(SD+1)/2:Ndy+(SD-1)/2);

%% Initial cell distrbution
CSD = 500; % initial cell surface density per type, 1/mm2
N1 = round(CSD*c*Ncsx*c*Ncsy/1e6); % initial number of type 1 cells, (R)
N2 = round(CSD*c*Ncsx*c*Ncsy/1e6); % initial number of type 2 cells, (Y)
X = zeros(Ncx,Ncy,Ncz); % Cells
U = zeros(Ncx,Ncy,Ncz); % Nutrients taken up by cells
nn = 0;
while nn < N1
    i = floor(0.5*Ncx+Ncsx*(rand(1)-0.5)+1);
    j = floor(0.5*Ncy+Ncsy*(rand(1)-0.5)+1);
    k = 1;
    if X(i,j,k) == 0
        nn = nn+1;
        X(i,j,k) = 1;
        U(i,j,k) = rand(1)*alphaL;
    end
end
nn = 0;
while nn < N2
    i = floor(0.5*Ncx+Ncsx*(rand(1)-0.5)+1);
    j = floor(0.5*Ncy+Ncsy*(rand(1)-0.5)+1);
    k = 1;
    if X(i,j,k) == 0
        nn = nn+1;
        X(i,j,k) = 2;
        U(i,j,k) = rand(1)*alphaA;
    end
end

%% Cell-growth time-course
tau0 = 0; % in hours
tauf = 120; % in hours
dtau = 0.1; % in hours, cell growth update and uptake timescale
ts = 1; % in hours, sampling time for snapshots of sections

r1e = r1*dtau; % probability of reproduction in a time step
r2e = r2*dtau; % probability of reproduction in a time step

d1e = d1*dtau; % probability of death in a time step
d2e = d2*dtau; % probability of death in a time step

taurng = tau0:dtau:tauf;

X1m = zeros(size(taurng));
X2m = zeros(size(taurng));
X1lm = zeros(size(taurng));
X2lm = zeros(size(taurng));
SPLa = zeros(Naz,length(taurng));
SPLd = zeros(Ndz,length(taurng));
SPAa = zeros(Naz,length(taurng));
SPAd = zeros(Ndz,length(taurng));
UAccL = zeros(size(taurng)); % Accumulated nutrients in live cells
UAccA = zeros(size(taurng)); % Accumulated nutrients in live cells
UW1 = zeros(size(taurng)); % Wasted nutrients in dead cells
UW2 = zeros(size(taurng)); % Wasted nutrients in dead cells
SLAccA = zeros(size(taurng)); % Total lys nutrients in the agar region
SLAccC = zeros(size(taurng)); % Total lys nutrients in the cell region
SAAccA = zeros(size(taurng)); % Total ade nutrients in the agar region
SAAccC = zeros(size(taurng)); % Total ade nutrients in the cell region

ct = 0;
cS = 0;

Q = [1 0; -1 0; 0 1; 0 -1; 1 -1; -1 1; 1 1; -1 -1]; % locations of eight neighbor grids
Qc = [1 0; -1 0; 0 1; 0 -1]; % locations of bud for confined cells

for tau = taurng
    ct = ct+1;
    tic

    % Update cell activity
    [zm,zi] = sort(rand(1,Ncz-1)); % random order for cells at different heights
    for z = zi
        zd = floor((z-1)/SC)+1;
        [I,J] = find((X(:,:,z)==1)+(X(:,:,z)==2)+(X(:,:,z)==3)); % find all live cells at height z
        Ncc = length(I);
        [SS,OS] = sort(rand(1,Ncc));
        I = I(OS);
        J = J(OS);
        for cc = 1:Ncc

            Id = floor((I(cc)-1)/SC)+1; % diffusion index along I
            Jd = floor((J(cc)-1)/SC)+1; % diffusion index along J
            xd = I(cc); % location of cell along x
            yd = J(cc); % location of cell along y

            live = 1;
            cd = rand(1);
            if (X(xd,yd,z) == 1)&&(cd < d1e)
                X(xd,yd,z) = 0.5;
                live = 0;
            end
            if (X(xd,yd,z) == 2)&&(cd < d2e)
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
                    if xb > Ncx, xb = xb-Ncx; end
                    if xb < 1, xb = xb+Ncx; end
                    yb = J(cc) + Q(qi(Natt),2);
                    if yb > Ncy, yb = yb-Ncy; end
                    if yb < 1, yb = yb+Ncy; end

                    if X(xb,yb,z)==0 % if available space, bud into it
                        zb = z;
                        Budded = 1;
                        X(xb,yb,zb) = X(xd,yd,z);
                        U(xb,yb,zb) = 0;
                    end
                end
                if (Budded==0) % try extended neighborhood in the same plane
                    rg = rg0;
                    xdp = xd;
                    ydp = yd;
                    if xd+rg > Ncx
                        Xp = [X(:,:,z); X(1:xd+rg-Ncx,:,z)];
                    else
                        if xd-rg < 1
                            Xp = [X(Ncx+xd-rg:Ncx,:,z); X(:,:,z)];
                            xdp = rg+1;
                        else
                            Xp = X(:,:,z);
                        end
                    end
                    if yd+rg > Ncy
                        Xp = [Xp, Xp(:,1:yd+rg-Ncy)];
                    else
                        if yd-rg < 1
                            Xp = [Xp(:,Ncy+yd-rg:Ncy), Xp];
                            ydp = rg+1;
                        end
                    end
                    Ng = Xp(xdp-rg:xdp+rg,ydp-rg:ydp+rg);
                    if sum(sum(abs(Ng)>0.25))<(2*rg+1)^2
                        [Ie,Je] = find(Ng==0);
                        [SSe,OSe] = min(sqrt((Ie-rg-1).^2+(Je-rg-1).^2));
                        if SSe <= rg
                            xe = Ie(OSe);
                            ye = Je(OSe);
                            TP = TracePath(rg+1,rg+1,xe,ye);
                            NT = size(TP,1);
                            TP(:,1) = TP(:,1) + xd-rg-1;
                            TP(:,1) = TP(:,1) + Ncx*((TP(:,1)<1)-(TP(:,1)>Ncx));
                            TP(:,2) = TP(:,2) + yd-rg-1;
                            TP(:,2) = TP(:,2) + Ncy*((TP(:,2)<1)-(TP(:,2)>Ncy));
                            for ctr = NT-1:-1:2
                                X(TP(ctr,1),TP(ctr,2),z) = X(TP(ctr-1,1),TP(ctr-1,2),z);
                                U(TP(ctr,1),TP(ctr,2),z) = U(TP(ctr-1,1),TP(ctr-1,2),z);
                            end
                            zb = z;
                            X(TP(NT,1),TP(NT,2),zb) = X(TP(NT-1,1),TP(NT-1,2),z);
                            U(TP(NT,1),TP(NT,2),zb) = U(TP(NT-1,1),TP(NT-1,2),z);
                            X(TP(1,1),TP(1,2),z) = X(xd,yd,z);
                            U(TP(1,1),TP(1,2),z) = 0;
                            Budded = 1;
                        end
                    end
                end
                if (Budded == 0) % bud to top or sides and push cells above
                    cq = rand(1);
                    if cq < 0.7
                        X(xd,yd,z+1:Ncz) = X(xd,yd,z:Ncz-1);
                        U(xd,yd,z+1:Ncz) = U(xd,yd,z:Ncz-1);
                    else
                        [qm3,qi3] = sort(rand(1,4)); % random index for neighboring grids
                        xb = I(cc) + Qc(qi3(1),1);
                        if xb > Ncx, xb = xb-Ncx; end
                        if xb < 1, xb = xb+Ncx; end
                        yb = J(cc) + Qc(qi3(1),2);
                        if yb > Ncy, yb = yb-Ncy; end
                        if yb < 1, yb = yb+Ncy; end
                        X(xb,yb,z+1:Ncz) = X(xb,yb,z:Ncz-1);
                        U(xb,yb,z+1:Ncz) = U(xb,yb,z:Ncz-1);
                        X(xb,yb,z) = celltype;
                        U(xb,yb,z) = 0;
                        xd = xb;
                        yd = yb;
                    end
                    zt = Ncz;
                    while (zt>z)&&(X(xd,yd,zt)==0)
                        zt = zt-1;
                    end
                    Xe = [X(Ncx,:,zt-1); X(:,:,zt-1); X(1,:,zt-1)];
                    Xe = [Xe(:,Ncy), Xe, Xe(:,1)];
                    if sum(sum(Xe(xd:xd+2,yd:yd+2)>0.1))<=6
                        [xe,ye] = find(Xe(xd:xd+2,yd:yd+2)==0);
                        et = floor(rand(1)*length(xe))+1;
                        xet = xd + xe(et) - 2;
                        xet = xet + Ncx*((xet<1)-(xet>Ncx));
                        yet = yd + ye(et) - 2;
                        yet = yet + Ncy*((yet<1)-(yet>Ncy));
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
        X1m(z,ct) = X1lm(z,ct) + sum(sum(X(:,:,z)==0.5));
        X2m(z,ct) = X2lm(z,ct) + sum(sum(X(:,:,z)==1.5));
    end
    UAccL(ct) = sum(sum(sum(U.*(X==1))));
    UAccA(ct) = sum(sum(sum(U.*(X==2))));
    UW1(ct) = sum(sum(sum(U.*(X==0.5))));
    UW2(ct) = sum(sum(sum(U.*(X==1.5))));
    SLAccC(ct) = sum(sum(sum(SLd)));
    SLAccA(ct) = sum(sum(sum(SLa)));
    SAAccC(ct) = sum(sum(sum(SAd)));
    SAAccA(ct) = sum(sum(sum(SAa)));
    disp([tau  sum(X1lm(:,ct)) sum(X2lm(:,ct)) sum(X1m(:,ct)+X2m(:,ct))])

    % Update diffusion constant
    Dd = zeros(Ndx,Ndy,Ndz);
    for ii = 1:SC
        for jj = 1:SC
            for kk = 1:SC
                Dd = Dd + D1/SC^3*(X(ii:SC:Ncx,jj:SC:Ncy,kk:SC:Ncz)>0.1);
            end
        end
    end

    dt0 = min(0.05*d^2/D0,0.05*g^2/D1); % in seconds
    trng = linspace(0,3600*dtau-dt0,round(3600*dtau/dt0)); % in seconds
    dt = trng(2)-trng(1);

    %% Diffusion loop
    for t = trng
        % Update concentration at the boundary, SLb
        SLam = [SLa(Nax,:,Naz); SLa(:,:,Naz); SLa(1,:,Naz)];
        SLam = [SLam(:,Nay), SLam, SLam(:,1)];
        SLae = interp2(SD*(0:Nax+1),SD*(0:Nay+1),SLam',xxd,yyd,'linear')';
        SLb = SLb - dt/h * (D0/d*(SLb - SLae) + 1/g*Dd(:,:,1).*(SLb - SLd(:,:,1)));

        % Diffusion in the agar region
        SLbm = zeros(Nax,Nay);
        for ii = 1:SD
            for jj = 1:SD
                SLbm = SLbm + 1/SD^2*SLb(ii:SD:Ndx,jj:SD:Ndy);
            end
        end
        SLab = SLa(:,:,3:Naz);
        SLab(:,:,Naz-1) = SLbm;
        dSLa = dt/d^2 * D0*([SLa(Nax,:,2:Naz); SLa(1:Nax-1,:,2:Naz)] + [SLa(2:Nax,:,2:Naz); SLa(1,:,2:Naz)] + ...
            [SLa(:,Nay,2:Naz), SLa(:,1:Nay-1,2:Naz)] + [SLa(:,2:Nay,2:Naz), SLa(:,1,2:Naz)] + ...
            SLa(:,:,1:Naz-1) + SLab - 6*SLa(:,:,2:Naz));
        SLa(:,:,2:Naz) = SLa(:,:,2:Naz) + dSLa;
        SLa(:,:,1) = SLa(:,:,2);

        % Diffusion in the community region
        SLdb = zeros(Ndx,Ndy,Ndz-1);
        SLdb(:,:,1) = SLb;
        SLdb(:,:,2:Ndz-1) = SLd(:,:,1:Ndz-2);
        dSLd = dt/g^2 * Dd(:,:,1:Ndz-1).*([SLd(Ndx,:,1:Ndz-1); SLd(1:Ndx-1,:,1:Ndz-1)] + ...
            [SLd(2:Ndx,:,1:Ndz-1); SLd(1,:,1:Ndz-1)] + [SLd(:,Ndy,1:Ndz-1), SLd(:,1:Ndy-1,1:Ndz-1)] + ...
            [SLd(:,2:Ndy,1:Ndz-1), SLd(:,1,1:Ndz-1)] + SLdb + SLd(:,:,2:Ndz) - 6*SLd(:,:,1:Ndz-1)) + dt/g^2 * (...
            ([Dd(2:Ndx,:,1:Ndz-1); Dd(1,:,1:Ndz-1)]-Dd(:,:,1:Ndz-1)).*([SLd(2:Ndx,:,1:Ndz-1); SLd(1,:,1:Ndz-1)]-SLd(:,:,1:Ndz-1))+...
            ([Dd(:,2:Ndy,1:Ndz-1), Dd(:,1,1:Ndz-1)]-Dd(:,:,1:Ndz-1)).*([SLd(:,2:Ndy,1:Ndz-1), SLd(:,1,1:Ndz-1)]-SLd(:,:,1:Ndz-1))+...
            (Dd(:,:,2:Ndz)-Dd(:,:,1:Ndz-1)).*(SLd(:,:,2:Ndz)-SLd(:,:,1:Ndz-1)));
        SLd(:,:,1:Ndz-1) = SLd(:,:,1:Ndz-1) + dSLd;
        SLd(:,:,Ndz) = SLd(:,:,Ndz-1);

        % Update concentration at the boundary, SAb
        SAam = [SAa(Nax,:,Naz); SAa(:,:,Naz); SAa(1,:,Naz)];
        SAam = [SAam(:,Nay), SAam, SAam(:,1)];
        SAae = interp2(SD*(0:Nax+1),SD*(0:Nay+1),SAam',xxd,yyd,'linear')';
        SAb = SAb - dt/h * (D0/d*(SAb - SAae) + 1/g*Dd(:,:,1).*(SAb - SAd(:,:,1)));

        % Diffusion in the agar region
        SAbm = zeros(Nax,Nay);
        for ii = 1:SD
            for jj = 1:SD
                SAbm = SAbm + 1/SD^2*SAb(ii:SD:Ndx,jj:SD:Ndy);
            end
        end
        SAab = SAa(:,:,3:Naz);
        SAab(:,:,Naz-1) = SAbm;
        dSAa = dt/d^2 * D0*([SAa(Nax,:,2:Naz); SAa(1:Nax-1,:,2:Naz)] + [SAa(2:Nax,:,2:Naz); SAa(1,:,2:Naz)] + ...
            [SAa(:,Nay,2:Naz), SAa(:,1:Nay-1,2:Naz)] + [SAa(:,2:Nay,2:Naz), SAa(:,1,2:Naz)] + ...
            SAa(:,:,1:Naz-1) + SAab - 6*SAa(:,:,2:Naz));
        SAa(:,:,2:Naz) = SAa(:,:,2:Naz) + dSAa;
        SAa(:,:,1) = SAa(:,:,2);

        % Diffusion in the community region
        SAdb = zeros(Ndx,Ndy,Ndz-1);
        SAdb(:,:,1) = SAb;
        SAdb(:,:,2:Ndz-1) = SAd(:,:,1:Ndz-2);
        dSAd = dt/g^2 * Dd(:,:,1:Ndz-1).*([SAd(Ndx,:,1:Ndz-1); SAd(1:Ndx-1,:,1:Ndz-1)] + ...
            [SAd(2:Ndx,:,1:Ndz-1); SAd(1,:,1:Ndz-1)] + [SAd(:,Ndy,1:Ndz-1), SAd(:,1:Ndy-1,1:Ndz-1)] + ...
            [SAd(:,2:Ndy,1:Ndz-1), SAd(:,1,1:Ndz-1)] + SAdb + SAd(:,:,2:Ndz) - 6*SAd(:,:,1:Ndz-1)) + dt/g^2 * (...
            ([Dd(2:Ndx,:,1:Ndz-1); Dd(1,:,1:Ndz-1)]-Dd(:,:,1:Ndz-1)).*([SAd(2:Ndx,:,1:Ndz-1); SAd(1,:,1:Ndz-1)]-SAd(:,:,1:Ndz-1))+...
            ([Dd(:,2:Ndy,1:Ndz-1), Dd(:,1,1:Ndz-1)]-Dd(:,:,1:Ndz-1)).*([SAd(:,2:Ndy,1:Ndz-1), SAd(:,1,1:Ndz-1)]-SAd(:,:,1:Ndz-1))+...
            (Dd(:,:,2:Ndz)-Dd(:,:,1:Ndz-1)).*(SAd(:,:,2:Ndz)-SAd(:,:,1:Ndz-1)));
        SAd(:,:,1:Ndz-1) = SAd(:,:,1:Ndz-1) + dSAd;
        SAd(:,:,Ndz) = SAd(:,:,Ndz-1);

        %% Uptake
        dtu = dt/3600; % in hours
        % Nutrient uptake
        gammaL = interp1(HgL,gammaLH,SAd);
        for i1 = 1:SC
            for i2 = 1:SC
                for i3 = 1:SC
                    SLde(i1:SC:Ncx,i2:SC:Ncy,i3:SC:Ncz) = SLd + 1/dV*dtu*gammaL.*(X(i1:SC:Ncx,i2:SC:Ncy,i3:SC:Ncz)==2).*(1+1/alphaA*U(i1:SC:Ncx,i2:SC:Ncy,i3:SC:Ncz));
                    SAde(i1:SC:Ncx,i2:SC:Ncy,i3:SC:Ncz) = SAd + 1/dV*dtu*gammaA*(X(i1:SC:Ncx,i2:SC:Ncy,i3:SC:Ncz)==1).*(1+1/alphaL*U(i1:SC:Ncx,i2:SC:Ncy,i3:SC:Ncz));
                end
            end
        end

        UcL = (dtu*vmL*(1+1/alphaL*U).*(X==1).*(SLde.^nL)./(SLde.^nL+KmmL^nL)).*(SLde>0);
        UcA = (dtu*vmA*(1+1/alphaA*U).*(X==2).*(SAde.^nA)./(SAde.^nA+KmmA^nA)).*(SAde>0);
        U = U+UcL+UcA;
        SLde = SLde-1/dV*UcL;
        SAde = SAde-1/dV*UcA;

        SLd = zeros(Ndx,Ndy,Ndz);
        SAd = zeros(Ndx,Ndy,Ndz);
        for i1 = 1:SC
            for i2 = 1:SC
                for i3 = 1:SC
                    SLd = SLd + 1/SC^3*SLde(i1:SC:Ncx,i2:SC:Ncy,i3:SC:Ncz);
                    SAd = SAd + 1/SC^3*SAde(i1:SC:Ncx,i2:SC:Ncy,i3:SC:Ncz);
                end
            end
        end

    end
    SPLa(:,ct) = shiftdim(mean(mean(SLa,2),1));
    SPLd(:,ct) = shiftdim(mean(mean(SLd,2),1));
    SPAa(:,ct) = shiftdim(mean(mean(SAa,2),1));
    SPAd(:,ct) = shiftdim(mean(mean(SAd,2),1));

    if mod(tau+0.03,ts)<0.1
		cS = cS+1;
        save(strcat(savename,'_ts1_cs',num2str(cS)),'X','tau','SLd','SAd','SLa','SAa')
    end
    toc
    
end

figure
semilogy(linspace(0,tau-dtau,ct-1),sum(X1m(:,1:ct-1)),'r')
hold on
semilogy(linspace(0,tau-dtau,ct-1),sum(X2m(:,1:ct-1)),'g')
semilogy(linspace(0,tau-dtau,ct-1),sum(X1lm(:,1:ct-1)),'r:')
semilogy(linspace(0,tau-dtau,ct-1),sum(X2lm(:,1:ct-1)),'g:')
xlabel('Time (hours)')
ylabel('Population (# of cells)')
xlim([0 tau-dtau])
% legend('Live R','Live Y')

save(savename)
