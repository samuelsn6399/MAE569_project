addpath("J2000PlanetaryModel/")
load("J200.mat")

% TODO gauss TOF needs to do long way around

% TODO probe is NOT allowed into the center of the sun!

% TODO convert sams plotting routine to functions :)

% unknowns:
% time: launch, DSM, flyby, arrival
% positions(3): launch, DSM, flyby, arrival
% velocity(3): launch, DSM start/end, flyby start/end, arrival
% relations:
% tDSM-tlaunch, rlaunch, rDSM, vlaunch, vDSM1 (know 3 find 2)
% tflyby-tDSM, rDSM, rflyby, vDSM2, vflyby1 (know 3 find 2)
% tarrival-tflyby, rflyby, rarrival, vflyby2, varrival (know 3 find 2)
% rlaunch is a function of time
% rflyby is a function of time
% rarrival is a function of time
% vf1 and vf2 set single-var bounded options

% picks: tlaunch, tDSM
% rDSM, vlaunch, vDSM1 (know 3 find 2)
% tflyby-tDSM, rDSM, rflyby, vDSM2, vflyby1 (know 3 find 2)
% tarrival-tflyby, rflyby, rarrival, vflyby2, varrival (know 3 find 2)

% picks: tlaunch, tDSM, tflyby, tarrival, search vf1->vf2 (single,
% bounded), vlaunch (triple but assume single)
% vlaunch, vDSM1 (know 3 find 2)
% rDSM, vDSM2 (know 3 find 2) -> kepler
% vflyby2, varrival (know 3 find 2) -> gauss
% could do (rDSM, vDSM1) kepler, (vlaunch, vDSM1) gauss
% vlaunch(3) or rDSM(3) -- which one has the smaller search space?
% for the sake of getting this done, assume vlaunch tan to earth. Lin
% search. But we already fixed rDSM.

% define tlaunch, dtDSM, dtflyby, dtarrival (in series)
% total time = dtDSM + dtflyby + dtarrival
% dtDSM: 4 to 12 months seems reasonable?
% dtflyby: 9 months past launch to 24 months past launch
% dtflyby: 9 to 24 months + -dtDSM
% dtarrival hohmann transfer*(0.5 to 2)
% vf1->vf2 has set bounds (min and max alt). Might want to nonlin space rp.
% or linspace the delta? nah just linspace the vectors

iearth = 3;
ijptr = 5;
reartho = r_earth() + 500;
rjptro = r_jptr() + 10000;

sweepn = 5;
% tlaunch = 0; % TODO this should a function of indexlaunch
thohmannjptr = pi*sqrt(((J2000(iearth).a+J2000(ijptr).a)/2)^3/1);

mindvtotalkms = 1/0;
mindtDSM = 0;
mindtflyby = 0;
mindtarrival = 0;
mindelmix = 0;
minindexlaunch = 0;


% TODO copying stuff from Sam, package into functions?
figure(1);clf(1);hold on; grid on; legend('Location','bestoutside');
plot3(0,0,0,'Marker','.','MarkerSize',50,'Color',[1 0.8 0],'DisplayName','Sun')
xlabel('X [AU]')
ylabel('Y [AU]')
zlabel('Z [AU]')
title('Orbital Visualization')
subtitle('x-y Plane Defined By Earth Orbit')
XaxisCenter = 0e0;
XaxisOffset = 5e0;
YaxisCenter = 1e0;
YaxisOffset = 5e0;
ZaxisCenter = 1e0;
ZaxisOffset = 1e0;
view(45,45)
res = 300; % TODO copied val



for indexlaunch = round(linspace(1, 400, 20))
    % get known positions
    rearththrow = J2000(iearth).r(indexlaunch,:);
    vearththrow = J2000(iearth).v(indexlaunch,:);
    rjptrthrow = J2000(ijptr).r(indexlaunch,:);
    vjptrthrow = J2000(ijptr).v(indexlaunch,:);

    if exist('eartho', 'var')
        delete(eartho)
        delete(jptro)
    end
    earthperiod = 2*pi*sqrt(J2000(iearth).a^3/1);
    tplearth = linspace(0,earthperiod,res);
    rplearth = zeros(res,3);
    for index=1:res
        rplearth(index,:) = kepler_univar(rearththrow, vearththrow, tplearth(index), 1)';
    end
    eartho = plot3(rplearth(:,1), rplearth(:,2), rplearth(:,3),'LineStyle','-','MarkerSize',10,'Color',[0 0 1],'DisplayName','Earth Orbit');

    jptrperiod = 2*pi*sqrt(J2000(ijptr).a^3/1);
    tpljptr = linspace(0,jptrperiod,res);
    rpljptr = zeros(res,3);
    for index=1:res
        rpljptr(index,:) = kepler_univar(rjptrthrow, vjptrthrow, tpljptr(index), 1)';
    end
    jptro = plot3(rpljptr(:,1), rpljptr(:,2), rpljptr(:,3),'LineStyle','-','MarkerSize',10,'Color',[1 0 0],'DisplayName','Jupiter Orbit');

    for dtDSM=linspace(daystoTU(4*30), daystoTU(12*30), sweepn)
        for dtflyby=linspace(max(daystoTU(9*30),dtDSM), daystoTU(18*30), sweepn) - dtDSM
            for dtarrival=linspace(0.5,2,sweepn)*thohmannjptr
                [rearthflyby, vearthflyby] = kepler_univar(rearththrow, vearththrow, dtDSM+dtflyby, 1);
                [rjptrarrive, vjptrarrive] = kepler_univar(rjptrthrow, vjptrthrow, dtDSM+dtflyby+dtarrival, 1);

                % earth-jupiter transfer
                [vflyby2, varrival] = GaussToF(rearthflyby, rjptrarrive, dtarrival, true);
                if any(cross(vflyby2, rearthflyby).*[0 0 1]<0)
                    [vflyby2, varrival] = GaussToF(rearthflyby, rjptrarrive, dtarrival, false);
                end
                if any(isnan(vflyby2)); break; end % just skip bad values
                vinfjptr = norm(varrival-vjptrarrive);
                vinfjptrkms = AUtokm(vinfjptr)/TUtoDays(3600*24);
                vbojptrkms = sqrt(vinfjptrkms^2 + 2*mu_jptr()/rjptro);
                vo5kms = sqrt(mu_jptr()/rjptro);
                dv3kms = vbojptrkms - vo5kms;
                if dv3kms > mindvtotalkms
                    fprintf("insertion energy too high\n")
                    break;
                end
                if norm(vflyby2) >= sqrt(2*1/norm(rearthflyby))
                    fprintf("probe shouldn't escape system\n")
                    break;
                end
                if norm(vflyby2)^2/2 - 1/norm(rearthflyby) > -1/2/J2000(ijptr).a
                    fprintf("probe had higher energy than jupiter\n")
                    break
                end
                % discard orbits that go lower than venus orbit (ish)
                % h constant, at peri h=rp*vp
                % E constant
                % cross(r,v) = rp*vp
                % E = v^2/2 - mu/r = vp^2/2 - mu/rp
                % vp = cross(r,v)/rp
                % (v^2/2 - mu/r)*rp^2 + mu*rp - cross(r,v)^2/2 = 0
                % rp = (-mu +/- sqrt(mu^2 - 4*(v^2/2 -
                % mu/r)*cross(r,v)^2/2))/(v^2/2 - mu/r)/2
                quada = (norm(vflyby2)^2/2 - 1/norm(rearthflyby));
                quadb = 1;
                quadc = -norm(cross(rearthflyby, vflyby2))^2/2;
                rpsun = (-quadb + sqrt(quadb^2 - 4*quada*quadc))/2/quada;
                if rpsun < 0.5
                    fprintf("Transfer peri low: %f AU\n",rpsun);
                    break;
                end

                % transfer to jupiter defines flyby options
                rpmaxkm = reartho; % km
                vinfearthfly2 = vflyby2 - vearthflyby;
                vinfearthfly2km = AUtokm(vinfearthfly2)/TUtoDays(3600*24); % I think that's right?
                vpmaxkms = sqrt(norm(vinfearthfly2km)^2 + 2*mu_earth()/rpmaxkm); % km/s TODO units!
                eccmax = rpmaxkm*vpmaxkms^2/mu_earth()-1;
                delmax = 2*asin(1/eccmax);
                delmin = 0;
                for mixflybysweep=linspace(0,1,sweepn)
                    delmix = mixflybysweep*delmin + (1-mixflybysweep)*delmax;
                    normal = cross(vinfearthfly2, vearthflyby);
                    normal = normal/norm(normal);
                    W = [0, -normal(3), normal(2); normal(3), 0, -normal(1); -normal(2), normal(1), 0];
                    rot = eye(3) + sin(delmix)*W + 2*sin(delmix/2)^2*W^2;
                    vinfearthfly1 = (rot\vinfearthfly2')';
                    vflyby1 = vinfearthfly1 + vearthflyby;
                    [rDSM, vDSM2] = kepler_univar(rearthflyby, vflyby1, -dtflyby, 1);
                    if any(~isreal(rDSM))
                        break;
                    end
                    if norm(rDSM) > 0.75*J2000(ijptr).a
                        fprintf("DSM too high\n");
                    end
                    % TODO add check if going "wrong way" cross(rdsm,
                    % vdsm2) (3) < 0, but inly once both Gauss are working

                    % flyby determines DSM position
                    [vlaunch, vDSM1] = GaussToF(rearththrow, rDSM, dtDSM, true);
                    if any(cross(vlaunch, rearththrow).*[0 0 1] < 0)
                        [vlaunch, vDSM1] = GaussToF(rearththrow, rDSM, dtDSM, false);
                    end
                    dv2 = norm(vDSM2 - vDSM1);
                    dv2kms = AUtokm(dv2)/TUtoDays(3600*24);
                    vinflaunch = norm(vlaunch-vearththrow);
                    vinflaunchkms = AUtokm(vinflaunch)/TUtoDays(3600*24);
                    vbolaunchkms = sqrt(vinflaunchkms^2 + 2*mu_earth()/reartho);
                    vo0kms = sqrt(mu_earth()/reartho);
                    dv1kms = vbolaunchkms - vo0kms;

                    dvtotalkms = dv1kms + dv2kms + dv3kms;
                    if dvtotalkms < mindvtotalkms
                        mindvtotalkms = dvtotalkms;
                        mindtDSM = dtDSM;
                        mindtflyby = dtflyby;
                        mindtarrival = dtarrival;
                        minindexlaunch = indexlaunch;
                        mindelmix = delmix;
                        fprintf("new best!\n")
                    end
                end
            end
            if exist('markearth', 'var')
                delete(markearth)
                delete(markdsm)
                delete(markflyby)
                delete(markjptr)
            end
            markearth = plot3(rearththrow(1),rearththrow(2), rearththrow(3),'Marker','.','MarkerSize',20,'Color',[0 0 1],'DisplayName','Launch');
            markdsm = plot3(rDSM(1),rDSM(2), rDSM(3),'Marker','.','MarkerSize',20,'Color',[0 1 0],'DisplayName','DSM');
            markflyby = plot3(rearthflyby(1),rearthflyby(2), rearthflyby(3),'Marker','.','MarkerSize',20,'Color',[0 1 0],'DisplayName','Flyby');
            markjptr = plot3(rjptrarrive(1),rjptrarrive(2), rjptrarrive(3),'Marker','.','MarkerSize',20,'Color',[1 0 0],'DisplayName','Arrival');
            drawnow
        end
    end
end

delete(markearth)
delete(markdsm)
delete(markflyby)
delete(markjptr)
delete(eartho)
delete(jptro)


dtDSM = mindtDSM;
dtflyby = mindtflyby;
dtarrival = mindtarrival;
indexlaunch = minindexlaunch;
delmix = mindelmix;

rearththrow = J2000(iearth).r(indexlaunch,:);
vearththrow = J2000(iearth).v(indexlaunch,:);
rjptrthrow = J2000(ijptr).r(indexlaunch,:);
vjptrthrow = J2000(ijptr).v(indexlaunch,:);
[rearthflyby, vearthflyby] = kepler_univar(rearththrow, vearththrow, dtDSM+dtflyby, 1);
[rjptrarrive, vjptrarrive] = kepler_univar(rjptrthrow, vjptrthrow, dtDSM+dtflyby+dtarrival, 1);
[vflyby2, varrival] = GaussToF(rearthflyby, rjptrarrive, dtarrival, true);
if any(cross(vflyby2, rearthflyby).*[0 0 1]<0)
    fprintf("printing flip, transfer\n")
    [vflyby2, varrival] = GaussToF(rearthflyby, rjptrarrive, dtarrival, false);
end
[rDSM, vDSM2] = kepler_univar(rearthflyby, vflyby1, -dtflyby, 1);
[vlaunch, vDSM1] = GaussToF(rearththrow, rDSM, dtDSM);
if any(cross(vlaunch, vearththrow).*[0 0 1] < 0)
    fprintf("printing flip, launch->DSM\n")
    [vlaunch, vDSM1] = GaussToF(rearththrow, rDSM, dtDSM, false);
end


plot3(rearththrow(1),rearththrow(2), rearththrow(3),'Marker','.','MarkerSize',20,'Color',[0 0 1],'DisplayName','Launch')
plot3(rDSM(1),rDSM(2), rDSM(3),'Marker','.','MarkerSize',20,'Color',[0 1 0],'DisplayName','DSM')
plot3(rearthflyby(1),rearthflyby(2), rearthflyby(3),'Marker','.','MarkerSize',20,'Color',[0 1 0],'DisplayName','Flyby')
plot3(rjptrarrive(1),rjptrarrive(2), rjptrarrive(3),'Marker','.','MarkerSize',20,'Color',[1 0 0],'DisplayName','Arrival')



tplo1 = linspace(0,dtDSM,res);
rplo1 = zeros(res,3);
for index=1:res
    rplo1(index,:) = kepler_univar(rearththrow, vlaunch, tplo1(index), 1)';
    if norm(rplo1(index,:)) > 5
        fprintf("huh?\n")
        % TODO for some reason it starts just randomy throwing in massive
        % values???? It seems like it just really hates certain dt values.
        % I wonder why?
    end
end
plot3(rplo1(:,1), rplo1(:,2), rplo1(:,3),'LineStyle','-','MarkerSize',10,'Color',[0.5 0.8 0],'DisplayName','Sat1 Orbit')
% TODO seems to be an issue with the sat 1 orbit?

tplo1 = linspace(0,dtflyby,res);
rplo1 = zeros(res,3);
for index=1:res
    rplo1(index,:) = kepler_univar(rDSM, vDSM2, tplo1(index), 1)';
end
plot3(rplo1(:,1), rplo1(:,2), rplo1(:,3),'LineStyle','-','MarkerSize',10,'Color',[0 1 0],'DisplayName','Sat2 Orbit')

tplo1 = linspace(0,dtarrival,res);
rplo1 = zeros(res,3);
for index=1:res
    rplo1(index,:) = kepler_univar(rearthflyby, vflyby2, tplo1(index), 1)';
end
plot3(rplo1(:,1), rplo1(:,2), rplo1(:,3),'LineStyle','-','MarkerSize',10,'Color',[0 0.8 0.5],'DisplayName','Sat3 Orbit')

earthperiod = 2*pi*sqrt(J2000(iearth).a^3/1);
tplearth = linspace(0,earthperiod,res);
rplearth = zeros(res,3);
for index=1:res
    rplearth(index,:) = kepler_univar(rearththrow, vearththrow, tplearth(index), 1)';
end
eartho = plot3(rplearth(:,1), rplearth(:,2), rplearth(:,3),'LineStyle','-','MarkerSize',10,'Color',[0 0 1],'DisplayName','Earth Orbit');

jptrperiod = 2*pi*sqrt(J2000(ijptr).a^3/1);
tpljptr = linspace(0,jptrperiod,res);
rpljptr = zeros(res,3);
for index=1:res
    rpljptr(index,:) = kepler_univar(rjptrthrow, vjptrthrow, tpljptr(index), 1)';
end
jptro = plot3(rpljptr(:,1), rpljptr(:,2), rpljptr(:,3),'LineStyle','-','MarkerSize',10,'Color',[1 0 0],'DisplayName','Jupiter Orbit');

xlim([-6,6])
ylim([-6,6])
zlim([-1,1])