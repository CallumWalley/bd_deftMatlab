tic;
path_runs = '../Kaipara_Mx6Nx2_keps_new3/';
runid_keps = '20_12_Mx6Nx2_keps_new3';
folder_exp = 'Test_ExportPred/';
temp_array_store = 'temp_array_store';
mkdir(temp_array_store);
addpath('nz_mapping');
addpath('Functions')
addpath(genpath('matlab'));
expfig = 0;
exp = 1; % export output as MATfile
fid = 1;

%Distances of FD below water surface
delta_drift1day1_obs = 2.1;
delta_drift1day2_obs = 2.0;
delta_drift1day3_obs = 2.0;
delta_drift2day2_obs = 1.8;
delta_drift2day3_obs = 1.6;
delta_drift3day2_obs = 1.7;
delta_drift3day3_obs = 1.9;
%Distance between instrument (Aquadopp) and measurements
delta2_obs = 0.35;

%Limits of caxis for the colorplots showing the Lagrangian distribution of
%eps.
sc_l0 = -10.4;
sc_lend = 0;

xl_0 = 0;
xl_end = 18000;
yl_eps_0 = -6.8;
yl_eps_end = 0;

kpr_bnd = load('Kaipara/kaipara_bnd.txt');

%Load observations data
delta = 1; % 1;% 10;%delta of output (to optimize computation time)
obs = load('observations_vBD/AllData_days1and2and3.mat'); % obs = load('..\lagrangian_observations_vBD\AllData_days1and2and3.mat');
day1_obs = find(obs.alldata.daynum(1:delta:end) == 1);
day2_obs = find(obs.alldata.daynum(1:delta:end) == 2);
day3_obs = find(obs.alldata.daynum(1:delta:end) == 3);
drifter1 = find(obs.alldata.dn(1:delta:end) == 1);
drifter2 = find(obs.alldata.dn(1:delta:end) == 2);
drifter3 = find(obs.alldata.dn(1:delta:end) == 3);
drift1day1 = intersect(day1_obs, drifter1);
drift1day2 = intersect(day2_obs, drifter1);
drift1day3 = intersect(day3_obs, drifter1);
drift2day1 = intersect(day1_obs, drifter2);
drift2day2 = intersect(day2_obs, drifter2);
drift2day3 = intersect(day3_obs, drifter2);
drift3day1 = intersect(day1_obs, drifter3);
drift3day2 = intersect(day2_obs, drifter3);
drift3day3 = intersect(day3_obs, drifter3);

yd_obs = obs.alldata.yd; yd_obs = yd_obs(1:delta:end);
dist_obs = obs.alldata.dist; dist_obs = dist_obs(1:delta:end);
eps_obs = obs.alldata.eps; eps_obs = eps_obs(1:delta:end);
depth_obs = obs.alldata.depth; depth_obs = depth_obs(1:delta:end);
vel_obs = obs.alldata.vel; vel_obs = vel_obs(1:delta:end);
[x_obs, y_obs] = lltonztm(-obs.alldata.lat, obs.alldata.lon);
x_obs = x_obs(1:delta:end); y_obs = y_obs(1:delta:end);

%Distance below water surface at which measurements
offset_depth_obs = NaN(length(yd_obs), 1);
offset_depth_obs(drift1day1) = delta_drift1day1_obs + delta2_obs;
offset_depth_obs(drift1day2) = delta_drift1day2_obs + delta2_obs;
offset_depth_obs(drift1day3) = delta_drift1day3_obs + delta2_obs;
offset_depth_obs(drift2day2) = delta_drift2day2_obs + delta2_obs;
offset_depth_obs(drift2day3) = delta_drift2day3_obs + delta2_obs;
offset_depth_obs(drift3day2) = delta_drift3day2_obs + delta2_obs;
offset_depth_obs(drift3day3) = delta_drift3day3_obs + delta2_obs;

%Load model predictions
Out = qpfopen([path_runs, 'trim-kaipara_river_cal_3d_vBD_v', num2str(runid_keps), '.dat']); % Out = qpfopen([path_runs,'trim-kaipara_river_cal_3d_vBD_v', num2str(runid_keps), '.dat']);
% Out = qpfopen(['..\d3d_model_vBD_v3\Test09\trim-kaipara_river_cal_3d_vBD_v',num2str(runid_keps),'.dat']); % Out = qpfopen([path_runs,'trim-kaipara_river_cal_3d_vBD_v',num2str(runid_keps),'.dat']);
[DataFields, Dims, NVal] = qpread(Out);
s = qpread(Out, 'vertical eddy viscosity', 'size');
tmax = s(1);
mmax = s(3);
nmax = s(4);
kmax = s(5);

% DataFields FOR k-eps ONLY !!
%     01. 'morphologic grid'
%     02. 'hydrodynamic grid'
%     03. 'grid'
%     04. 'open boundaries'
%     05. 'closed boundaries'
%     06. 'thin dams'
%     07. 'temporarily inactive water level points'
%     08. 'temporarily inactive velocity points'
%     09. 'parallel partition numbers'
%     10. '-------'
%     11. 'water level (when dry: bed level)'
%     12. 'water level'
%     13. 'water depth'
%     14. 'depth averaged velocity'
%     15. 'staggered depth averaged velocities'
%     16. 'horizontal velocity'
%     17. 'staggered horizontal velocity'
%     18. 'velocity'
%     19. 'vertical velocity'
%     20. 'velocity in depth averaged flow direction'
%     21. 'velocity normal to depth averaged flow direction'
%     22. 'filtered depth averaged velocity'
%     23. 'froude number'
%     24. 'head'
%     25. '-------'
%     26. 'density'
%     27. 'hydrostatic pressure'
%     28. 'relative hydrostatic pressure'
%     29. 'salinity'
%     30. 'temperature'
%     31. 'turbulent energy'
%     32. 'energy dissipation'
%     33. 'vertical eddy viscosity'
%     34. 'vertical eddy diffusivity'
%     35. 'horizontal viscosity'
%     36. 'richardson number'
%     37. '-------'
%     38. 'bed shear stress'
%     39. 'staggered bed shear stress'
%     40. 'maximum bed shear stress'
%     41. 'initial bed level'
%     42. 'bed level in water level points'
%     43. 'bed level in velocity points'
%     44. '-------'
%     45. 'grid cell surface area'
%     46. '-------'
%     47. 'W-omega per layer in zeta point'
%     48. '-------'
%     49. 'X-coord. water level point in local system'
%     50. 'Y-coord. water level point in local system'
%     51. '1/-1 Active/Non-active bottom point ( w.r.t. coordinates )'
%     52. '1/-1 Active/Non-active water level point (w.r.t. coordinates )'

% t1 = 2000;
% eps = qpread(Out,'energy dissipation','griddata',1:tmax); % [tmax,mmax,nmax,kmax] matrix
eps = qpread(Out, 'energy dissipation', 'griddata', 1:tmax);
hv = qpread(Out, 'horizontal velocity', 'griddata', 1:tmax);
% wl = qpread(Out,'water level','griddata',1:tmax);
wl = qpread(Out, 'water level', 'griddata', 1:tmax);
wd = qpread(Out, 'water depth', 'griddata', 1:tmax);
sal = qpread(Out, 'salinity', 'griddata', 1:tmax);
tke = qpread(Out, 'turbulent energy', 'griddata', 1:tmax);
nuV = qpread(Out, 'vertical eddy viscosity', 'griddata', 1:tmax);

grd = wlgrid('read', '..\InputFiles\kpra_rvr_v01_Mx3.grd');
depth_i = wldep('read', '..\InputFiles\Test03_Glen_bank_v02_Mx3.dep', grd);

x_model = reshape(eps.X(1, :, :, 2), mmax, nmax); % x_model = reshape(eps1.X(1,2:mmax,2:nmax,2),mmax-1,nmax-1);
y_model = reshape(eps.Y(1, :, :, 2), mmax, nmax); % y_model = reshape(eps1.Y(1,2:mmax,2:nmax,2),mmax-1,nmax-1);

[Y, Mo, D, H, Mi, S] = datevec(eps.Time);
yd_model = yearday_DrJM(Y, Mo, D + H / 24 + Mi / (24 * 60) + S / (24 * 60 * 60)); % yd_model = yearday(Y,Mo,D + H/24 + Mi/(24*60) + S/(24*60*60));

% idx_obs = find(yd_obs <= yd_model1(end));
% yd_obs1 = yd_obs(idx_obs);
% x_obs1 = x_obs(idx_obs);
% y_obs1 = y_obs(idx_obs);


sc = SlurmManager();
%sc.debug = 1;
sc.sfor(@proc_single, 1:length(idx_obs));

for tt = 1:length(idx_obs)
    load(fullfile(temp_array_store,'wks_', tt, '.mat'));
    eps_pred_lag.idxt(tt) = idxt;
    eps_pred_lag.yd_pred_lag(tt) = yd_lag_eps;
    eps_pred_lag.x_pred_lag(tt) = x_pred_lag;
    eps_pred_lag.y_pred_lag(tt) = y_pred_lag;
    eps_pred_lag.z_pred_lag(tt) = z_lag_eps;
    eps_pred_lag.allz_pred_lag_old(tt) = allz_lag_eps;
    eps_pred_lag.allz_pred_lag(tt) = allz_lag_eps2;
    eps_pred_lag.m_lag(tt) = m_lag;
    eps_pred_lag.n_lag(tt) = n_lag;
    eps_pred_lag.k_lag(tt) = k_lag_eps;
    eps_pred_lag.yd_lag(tt) = yd_lag_eps;
    eps_pred_lag.wl_lag(tt) = wl_lag;
    eps_pred_lag.wd_lag(tt) = wd_lag;
    eps_pred_lag.valallz(tt) = eps_lag_allz;
    eps_pred_lag.val(tt) = eps_lag;
end

eps_pred_lag.idx_obs = idx_obs; 

if exp == 1
    save([folder_exp, 'Runid', runid_keps, 'eps_pred_lag_1to', num2str(t1), '.mat'], '-struct', 'eps_pred_lag')
end
toc;

profsave


function proc_single(tt)
    [idxt, ~] = near(yd_model, yd_obs(tt)); %find the index of the time step of the model corresponding to the time of observations yd_obs(tt)
    %     wd_model_lag = reshape(wd_model(idxt,:,:),s1-1,s2-1); %water depth at time step of interest
    %     z_model = reshape(eps.Z(idxt(tt),:,:,k),mmax,nmax); %     z_model(tt,:,:) = reshape(eps.Z(idxt(tt),:,:,k),mmax,nmax);
    %     eps_val = reshape(eps.Val(idxt(tt),:,:,k),mmax,nmax); %find corresponding dissipation distribution (i.e. whole grid at time yd_obs(tt))
    %     wl_val = reshape(wl.Val(idxt(tt),:,:),mmax,nmax);
    yd_lag_eps = yd_model(idxt(tt));

    [m_lag, n_lag] = xy2mn(x_model, y_model, x_obs1(tt), y_obs1(tt)); %     [idx_mnrange(tt,1),idx_mnrange(tt,2)] = xy2mn(x_model,y_model,x_obs(tt),y_obs(tt));%     [m,n] = xy2mn(x_model,y_model,x_obs(tt),y_obs(tt));%find grid cell corresponding to observations
    %     test(tt) = dsearchn([x_model,y_model],[x_obs(tt),y_obs(tt)]);
    d_i_lag = depth_i(m_lag, n_lag);
    k_lag_eps = near(eps.Z(idxt(tt), m_lag, n_lag, :), max(eps.Z(idxt(tt), m_lag, n_lag, :)) - offset_depth_obs(tt)); %     k_lag_eps(tt) = near(eps.Z(idxt(tt),m_lag,n_lag,:),max(eps.Z(idxt(tt),m_lag,n_lag,:))-1.5);
    z_lag_eps = eps.Z(idxt(tt), m_lag, n_lag, k_lag_eps);
    allz_lag_eps = max(eps.Z(idxt(tt), m_lag, n_lag, :)) - min(eps.Z(idxt(tt), m_lag, n_lag, :));
    allz_lag_eps2 = eps.Z(idxt(tt), m_lag, n_lag, :);
    k_lag_vel = near(hv.Z(idxt(tt), m_lag, n_lag, :), max(hv.Z(idxt(tt), m_lag, n_lag, :)) - offset_depth_obs(tt)); %     k_lag_vel(tt) = near(hv.Z(idxt(tt),m_lag,n_lag,:),max(hv.Z(idxt(tt),m_lag,n_lag,:))-1.5);
    z_lag_vel= hv.Z(idxt(tt), m_lag, n_lag, k_lag_vel);
    allz_lag_vel = max(hv.Z(idxt(tt), m_lag, n_lag, :)) - min(hv.Z(idxt(tt), m_lag, n_lag, :));
    allz_lag_vel2 = max(hv.Z(idxt(tt), m_lag, n_lag, :)) - min(hv.Z(idxt(tt), m_lag, n_lag, :)) + abs(mean(diff(hv.Z(idxt(tt), m_lag, n_lag, :))));
    x_pred_lag = x_model(m_lag, n_lag); 
    y_pred_lag = y_model(m_lag, n_lag);
    %     sin_lag = sinuosity_grid(m_lag1(tt),2);
    %     sin_deg_lag(tt) = sinuosity_grid(m_lag1(tt),3);
    wl_lag = wl.Val(idxt(tt), m_lag, n_lag);
    wd_lag = wd.Val(idxt(tt), m_lag, n_lag);
    sal_lag = sal.Val(idxt(tt), m_lag, n_lag, k_lag_vel);
    sal_lag_allz = sal.Val(idxt(tt), m_lag, n_lag, :);
    hv_lag_xcomp = hv.XComp(idxt(tt), m_lag, n_lag, k_lag_vel);
    hv_lag_ycomp = hv.YComp(idxt(tt), m_lag, n_lag, k_lag_vel);
    hv_lag_xcomp_allz = hv.XComp(idxt(tt), m_lag, n_lag, :);
    hv_lag_ycomp_allz = hv.YComp(idxt(tt), m_lag, n_lag, :);
    eps_lag = eps.Val(idxt(tt), m_lag, n_lag, k_lag_eps);
    eps_lag_allz = eps.Val(idxt(tt), m_lag, n_lag, :);
    tke_lag = tke.Val(idxt(tt), m_lag, n_lag, k_lag_eps);
    tke_lag_allz = tke.Val(idxt(tt), m_lag, n_lag, :);
    nuV_lag = nuV.Val(idxt(tt), m_lag, n_lag, k_lag_eps);
    nuV_lag_allz = nuV.Val(idxt(tt), m_lag, n_lag, :);

    save(fullfile(temp_array_store,'wks_', tt, '.mat'),'yd_lag_eps','m_lag','n_lag','k_lag_eps','z_lag_vel','allz_lag_vel','allz_lag_vel2','x_pred_lag','y_pred_lag','z_lag_eps','allz_lag_eps','allz_lag_eps2','eps_lag','eps_lag_allz','idxt','d_i_lag','wl_lag','wd_lag','hv_lag_xcomp','hv_lag_ycomp','hv_lag_xcomp_allz','hv_lag_ycomp_allz','sal_lag','sal_lag_allz','tke_lag','tke_lag_allz','nuV_lag','nuV_lag_allz');
end