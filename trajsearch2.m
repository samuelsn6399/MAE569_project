reartho = 500;
rjptro = 10000;



sweepn = 5; % how many options do we check for each variabe?
sweepdatesn = 30; % how many dates do we check? This one behaves a little differently, which is why it's seperate

% structure: launch, DSM1, DSM2, flyby1, flyby2, arrival
paramslow.index = 9491;
paramshigh.index = paramslow.index + 400;
% launch is 0, unneeded, but keeps the indeces consistent. DSM and flyby are repeated bc they're the same for
% 1 and 2, but indeces consistent. arrival is given in multiples of hohmann
% transfer
paramslow.dt = [0, daystoTU(4*30), daystoTU(4*30), daystoTU(9*30), daystoTU(9*30), 0.5];
paramshigh.dt = [0, daystoTU(11*30), daystoTU(11*30), daystoTU(13*30), daystoTU(13*30), 2];

fid1 = figure(1);
fid2 = figure(2);

[mindv, params1] = optimizeE2J(paramslow, paramshigh, sweepn, sweepdatesn, reartho, rjptro, fid1, fid2)

paramslow.index = max(0, params1.index - 20);
paramshigh.index = params1.index + 20;
sweepn = 20;
sweepdatesn = 15;
fid1 = figure(3);
fid2 = figure(4);
[mindv, params2] = optimizeE2J(paramslow, paramshigh, sweepn, sweepdatesn, reartho, rjptro, fid1, fid2)


function [mindvtotalkms, params] = optimizeE2J(paramslow, paramshigh, sweepn, sweepdatesn, reartho, rjptro, fid1, fid2)

addpath("J2000PlanetaryModel/")
load("J200.mat", "J2000")

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

% picks: tlaunch, tDSM, tflyby, tarrival, search vf1->vf2 (single,
% bounded), vlaunch (triple but assume single)
% vlaunch, vDSM1 (know 3 find 2) -> gauss
% rDSM, vDSM2 (know 3 find 2) -> kepler
% vflyby2, varrival (know 3 find 2) -> gauss

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
reartho = r_earth() + reartho;
rjptro = r_jptr() + rjptro;


% find Hohmann transfer time as baseline
ahohmannjptr = (J2000(iearth).a+J2000(ijptr).a)/2;
thohmannjptr = pi*sqrt((ahohmannjptr)^3/1);

% define min vals to keep track of
mindvtotalkms = 1/0;
mindtDSM = 0;
mindtflyby = 0;
mindtarrival = 0;
mindelmix = 0;
minindexlaunch = 0;
minrDSM = [0 0 0];
minvlaunch = [0 0 0];
minvflyby1 = [0 0 0];
minvflyby2 = [0 0 0];
minvDSM2 = [0 0 0];
minvDSM1 = [0 0 0];
minvarrival = [0 0 0];


rDSM = [0, 0, 0]; % declare as a global var, i think something breaks without this lol
mindvlist = [];
dvlist = [];
iter = 0;

% plotting stuff...
figure(fid1);clf(fid1);
plot3([],[],[]);
ax1 = gca;
figure(fid2);clf(fid2);
plot([],[]);
ax2 = gca;

figure(fid1);hold on; grid on; legend('Location','bestoutside');
plot3(ax1, 0,0,0,'Marker','.','MarkerSize',50,'Color',[1 0.8 0],'DisplayName','Sun');
xlabel('X [AU]')
ylabel('Y [AU]')
zlabel('Z [AU]')
title('Orbital Visualization')
subtitle('x-y Plane Defined By Earth Orbit')
view(45,45)
res = 300;

% pre intializae some of this stuff for plotting
indexlaunch = 1;
rearththrow = J2000(iearth).r(indexlaunch,:);
vearththrow = J2000(iearth).v(indexlaunch,:);
rjptrthrow = J2000(ijptr).r(indexlaunch,:);
vjptrthrow = J2000(ijptr).v(indexlaunch,:);

% plot baseline orbits
earthperiod = 2*pi*sqrt(J2000(iearth).a^3/1);
eartho = plotorbit(ax1, rearththrow, vearththrow, earthperiod, res, 'Earth Orbit', [0 0 1]);
jptrperiod = 2*pi*sqrt(J2000(ijptr).a^3/1);
jptro = plotorbit(ax1, rjptrthrow, vjptrthrow, jptrperiod, res, 'Jupiter', [1 0 0]);

markearth = plot3(ax1, rearththrow(1),rearththrow(2), rearththrow(3),'Marker','.','MarkerSize',20,'Color',[0 0 1],'DisplayName','Launch');
markdsm = plot3(ax1, rearththrow(1),rearththrow(2), rearththrow(3),'Marker','.','MarkerSize',20,'Color',[0 1 0],'DisplayName','DSM');
markflyby = plot3(ax1, rearththrow(1),rearththrow(2), rearththrow(3),'Marker','.','MarkerSize',20,'Color',[0 1 0],'DisplayName','Flyby');
markjptr = plot3(ax1, rearththrow(1),rearththrow(2), rearththrow(3),'Marker','.','MarkerSize',20,'Color',[0 1 0],'DisplayName','Flyby');
drawnow

% main loop!
for indexlaunch = round(linspace(paramslow.index, paramshigh.index, sweepdatesn))
    % get known positions
    rearththrow = J2000(iearth).r(indexlaunch,:);
    vearththrow = J2000(iearth).v(indexlaunch,:);
    rjptrthrow = J2000(ijptr).r(indexlaunch,:);
    vjptrthrow = J2000(ijptr).v(indexlaunch,:);

    for dtDSM=linspace(paramslow.dt(2), paramshigh.dt(2), sweepn)
        for dtflyby=linspace(max(paramslow.dt(4),dtDSM), paramshigh.dt(4), sweepn) - dtDSM
            for dtarrival=linspace(paramslow.dt(6),paramshigh.dt(6),sweepn)*thohmannjptr
                [rearthflyby, vearthflyby] = kepler_univar(rearththrow, vearththrow, dtDSM+dtflyby, 1);
                [rjptrarrive, vjptrarrive] = kepler_univar(rjptrthrow, vjptrthrow, dtDSM+dtflyby+dtarrival, 1);

                % earth-jupiter transfer
                [vflyby2, varrival] = GaussToF(rearthflyby, rjptrarrive, dtarrival, true);
                if any(cross(rearthflyby, vflyby2).*[0 0 1]<0)
                    [vflyby2, varrival] = GaussToF(rearthflyby, rjptrarrive, dtarrival, false);
                    if any(cross(rearthflyby, vflyby2).*[0 0 1]<0)
                        % sometimes it's simply not possible to go the
                        % right way around in the given time
                        % TODO technically it could try a different dt
                        % here, but i'm not really sure if I have time for
                        % that now. Maybe try current dt + earth orbit
                        % time? basically re-target earth on a previous 
                        % go-around. Conditions that caused this bug:
                        % index: 300
                        % dtDSM: 4.6174
                        % dtflyby: 2.0371
                        % dtarrival: 34.3222
                        % delmix: 0.7485
                        continue;
                    end
                end
                [rtest, vtest] = kepler_univar(rearthflyby, vflyby2, dtarrival, 1);
                if max(abs(rtest-rjptrarrive)) > 1e-3 || max(abs(vtest-varrival)) > 1e-3
                    % doesn't need to be super accurate, just make sure we
                    % aren't trying to teleport/time travl...
                    % TODO i'm not sure why, but for some reason the
                    % optimizer keeps on exploiting some glitch in the
                    % Gauss solver that results in a teleporting
                    % spacecraft. The spacecraft travels on a
                    % reasonable trajectory, just ends up at the right
                    % place at the wrong time...
                    fprintf("arrive fail: GR[%4.2f %4.2f %4.2f] v KR[%4.2f %4.2f %4.2f]\n", rjptrarrive', rtest');
                    continue;
                end
                if any(isnan(vflyby2)); continue; end % just skip bad values
                vinfjptr = norm(varrival-vjptrarrive);
                vinfjptrkms = AUtokm(vinfjptr)/TUtoDays(3600*24);
                vbojptrkms = sqrt(vinfjptrkms^2 + 2*mu_jptr()/rjptro);
                vo5kms = sqrt(mu_jptr()/rjptro);
                dv3kms = vbojptrkms - vo5kms;
                if dv3kms > 0.9*mindvtotalkms % TODO maybe set to like 75% of min?
                    fprintf("insertion energy too high\n")
                    continue;
                end
                if norm(vflyby2)^2/2 - 1/norm(rearthflyby) > -1/2/J2000(ijptr).a
                    fprintf("probe had higher energy than jupiter\n")
                    continue
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
                    continue;
                end

                % transfer to jupiter defines flyby options
                rpmaxkm = reartho; % km
                vinfearthfly2 = vflyby2 - vearthflyby;
                vinfearthfly2km = AUtokm(vinfearthfly2)/TUtoDays(3600*24); % I think that's right?
                vpmaxkms = sqrt(norm(vinfearthfly2km)^2 + 2*mu_earth()/rpmaxkm); % km/s TODO units!
                eccmax = rpmaxkm*vpmaxkms^2/mu_earth()-1;
                delmax = 2*asin(1/eccmax);
                delmin = 0; % always an option to just skip the flyby

                for mixflybysweep=linspace(0,1,sweepn)
                    % TODO see if we can move this out one loop
                    delmix = mixflybysweep*delmin + (1-mixflybysweep)*delmax;
                    % three points define an orbital plane: earth, the sun,
                    % and jupiter. 
                    % plane must contain vflyby2, plane should contain
                    % jupiter r vector?
                    % 1. Probe may use all available del for slingshot (no
                    % plane change through del). This is the plane defined
                    % by vflyby2 and rjptrarrive.
                    % 2. Probe must use all available del (delmix =
                    % delmax).
                    % 3. Probe may use some of its turning del for plane
                    % change, but cannot "reverse" direction (i.e. probe
                    % should approach earth somewhere between the plane of
                    % earth and the final trajectory plane).
                    normal = cross(vflyby2, rjptrarrive);
                    % Rodriguez rotation formula
                    normal = normal/norm(normal);
                    W = [0, -normal(3), normal(2); normal(3), 0, -normal(1); -normal(2), normal(1), 0];
                    rot = eye(3) + sin(delmix)*W + 2*sin(delmix/2)^2*W^2;
                    vinfearthfly1 = (rot\vinfearthfly2')'; % inverse multiply since we want v1 from v2, could be faster to just use negative angle
                    vflyby1 = vinfearthfly1 + vearthflyby;
                    [rDSM, vDSM2] = kepler_univar(rearthflyby, vflyby1, -dtflyby, 1);
                    if any(~isreal(rDSM))
                        continue;
                    end
                    if norm(rDSM) > 0.75*J2000(ijptr).a
                        fprintf("DSM too high\n");
                    end
                    
                    % flyby determines DSM position
                    [vlaunch, vDSM1] = GaussToF(rearththrow, rDSM, dtDSM, true);
                    if any(cross(rearththrow, vlaunch).*[0 0 1] < 0)
                        [vlaunch, vDSM1] = GaussToF(rearththrow, rDSM, dtDSM, false);
                        if any(cross(rearththrow, vlaunch).*[0 0 1] < 0)
                            % see comment on previous
                            continue;
                        end
                    end
                    [rtest, vtest] = kepler_univar(rearththrow, vlaunch, dtDSM, 1);
                    if max(abs(rtest-rDSM)) > 1e-3 || max(abs(vtest-vDSM1)) > 1e-3
                        % see comment on previous
                        continue;
                    end
                    quada = (norm(vlaunch)^2/2 - 1/norm(rearththrow));
                    quadb = 1;
                    quadc = -norm(cross(rearththrow, vlaunch))^2/2;
                    rpsun = (-quadb + sqrt(quadb^2 - 4*quada*quadc))/2/quada;
                    if rpsun < 0.5
                        fprintf("Launch peri low: %f AU\n",rpsun);
                        continue;
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
                        minrDSM = rDSM;
                        minvlaunch = vlaunch;
                        minvflyby1 = vflyby1;
                        minvflyby2 = vflyby2;
                        minvDSM2 = vDSM2;
                        minvDSM1 = vDSM1;
                        minvarrival = varrival;
                        fprintf("new best!\n")
                    end
                    if dvtotalkms < 1.5*mindvtotalkms
                        dvlist = [dvlist; dvtotalkms];
                        mindvlist = [mindvlist; mindvtotalkms];
                    end
                    iter = iter + 1;
                end
            end
            % some plotting, hopefully not in a spot where it takes too
            % long... also only plot markers to show progress of optimizer,
            % since trajectories take longer to figure out
            delete(markearth)
            delete(markdsm)
            delete(markflyby)
            delete(markjptr)
            markearth = plot3(ax1, rearththrow(1),rearththrow(2), rearththrow(3),'Marker','.','MarkerSize',20,'Color',[0 0 1],'DisplayName','Launch');
            markdsm = plot3(ax1, rDSM(1),rDSM(2), rDSM(3),'Marker','.','MarkerSize',20,'Color',[0 1 0],'DisplayName','DSM');
            markflyby = plot3(ax1, rearthflyby(1),rearthflyby(2), rearthflyby(3),'Marker','.','MarkerSize',20,'Color',[0 1 0],'DisplayName','Flyby');
            markjptr = plot3(ax1, rjptrarrive(1),rjptrarrive(2), rjptrarrive(3),'Marker','.','MarkerSize',20,'Color',[1 0 0],'DisplayName','Arrival');
            drawnow
        end
        % no need to update the progress plot too often
        plot(ax2, [mindvlist, dvlist])
        drawnow
    end
end
% update optimizer plot to final version
loglog(ax2, [mindvlist, dvlist])
legend(ax2, "Minimum delta v", "Iteration delta v")
title(ax2, "Optimization Progress")
subtitle(ax2, "Only trajectories close to minimum are plotted")
xlabel(ax2, "Promising iteration number")
ylabel(ax2, "Delta V, km/s")

% clean chart
delete(markearth)
delete(markdsm)
delete(markflyby)
delete(markjptr)
delete(eartho)
delete(jptro)

% copy best vals
dtDSM = mindtDSM;
dtflyby = mindtflyby;
dtarrival = mindtarrival;
indexlaunch = minindexlaunch;
delmix = mindelmix;
vlaunch = minvlaunch;
vDSM2 = minvDSM2;
rDSM = minrDSM;
vflyby2 = minvflyby2;
vflyby1 = minvflyby1;
vDSM1 = minvDSM1;
varrival = minvarrival;

% find positions
rearththrow = J2000(iearth).r(indexlaunch,:);
vearththrow = J2000(iearth).v(indexlaunch,:);
rjptrthrow = J2000(ijptr).r(indexlaunch,:);
vjptrthrow = J2000(ijptr).v(indexlaunch,:);
[rearthflyby, ~] = kepler_univar(rearththrow, vearththrow, dtDSM+dtflyby, 1);
[rjptrarrive, ~] = kepler_univar(rjptrthrow, vjptrthrow, dtDSM+dtflyby+dtarrival, 1);

% plot orbits
plot3(ax1, rearththrow(1),rearththrow(2), rearththrow(3),'Marker','.','MarkerSize',20,'Color',[0 0 1],'DisplayName','Launch')
plot3(ax1, rDSM(1),rDSM(2), rDSM(3),'Marker','.','MarkerSize',20,'Color',[0 1 0],'DisplayName','DSM')
plot3(ax1, rearthflyby(1),rearthflyby(2), rearthflyby(3),'Marker','.','MarkerSize',20,'Color',[0 1 0],'DisplayName','Flyby')
plot3(ax1, rjptrarrive(1),rjptrarrive(2), rjptrarrive(3),'Marker','.','MarkerSize',20,'Color',[1 0 0],'DisplayName','Arrival')
plotorbit(ax1, rearththrow, vlaunch, dtDSM, res, "Launch->DSM",         [0 1  .0]);
plotorbit(ax1, rDSM, vDSM2, dtflyby, res, "DSM->Flyby",                 [0 .7 .3]);
plotorbit(ax1, rearthflyby, vflyby2, dtarrival, res, "Flyby->Arrival",  [0 .4 .6]);
earthperiod = 2*pi*sqrt(J2000(iearth).a^3/1);
eartho = plotorbit(ax1, rearththrow, vearththrow, earthperiod, res, "Earth", [0 0 1]);
jptrperiod = 2*pi*sqrt(J2000(ijptr).a^3/1);
jptro = plotorbit(ax1, rjptrthrow, vjptrthrow, jptrperiod, res, "Jupiter", [1 0 0]);

xlim(ax1, [-6,6])
axis(ax1, 'equal')

% return results
fprintf("Min delta v: %f\n", mindvtotalkms);

params.index = indexlaunch;
params.dt = [0, dtDSM, dtDSM, dtflyby, dtflyby, dtarrival*thohmannjptr];
params.r = [rearththrow; rDSM; rDSM; rearthflyby; rearthflyby; rjptrarrive];
params.v = [vlaunch; vDSM1; vDSM2; vflyby1; vflyby2; varrival];
params.delmix = delmix;
params.iter = iter;
params.maxiter = sweepdatesn*sweepn^4;

end