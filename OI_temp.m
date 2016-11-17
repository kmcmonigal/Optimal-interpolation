% optimally interpolate T and S using similar method as OI_vel

load('ACTGEMData_skype.mat');
load('micro_b_d_struct_despike_filter.mat');
load('micro_b_noise.mat');
load('micro_d_noise.mat');

coast_lat=-33.2910;
coast_lon=27.4783;
c2_pos = [27.5167,-33.4232];
c3_pos = [27.5698,-33.5112];
a_pos = [27.595,-33.5583];
b_pos = [27.6428,-33.6674];
c_pos = [27.7152,-33.7996];
d_pos = [27.8603,-34.0435];
e_pos = [28.0316,-34.29];
f_pos = [28.17,-34.5380];
g_pos = [28.3453,-34.8213];

instr_pos = [c2_pos;c3_pos;a_pos;b_pos;c_pos;d_pos;e_pos;f_pos;g_pos];

for i=1:9
    dist(i)=1000*sw_dist([instr_pos(1,2) instr_pos(i,2)],[instr_pos(1,1) instr_pos(i,1)],'km');
end

B_int.temp=micro_b.temp;
D_int.temp=micro_d.temp;

B_int.pres=micro_b.pres;
D_int.pres=micro_d.pres;

B_int.sal=micro_b.sal;
D_int.sal=micro_d.sal;

% need to get theta
B_int.theta=gsw_pt_from_t(B_int.sal,B_int.temp,B_int.pres);
D_int.theta=gsw_pt_from_t(D_int.sal,D_int.temp,D_int.pres);

B_temp_var=nanvar(B_int.theta); %variance of the T measurements, to create a ratio with the noise level
D_temp_var=nanvar(D_int.theta);

dx_obs=[dist(4);dist(4);dist(4);dist(4);dist(4);dist(6);dist(6);dist(6);dist(6)];

for i=1:size(micro_b.temp,1)
    dp_obs(:,i)=[B_int.pres(i,:).';D_int.pres(i,:).'];
    temp_obs(:,i)=[B_int.theta(i,:).';D_int.theta(i,:).'];
end
        
var_obs=[B_temp_var.';D_temp_var.'];
noise_obs=[noise_micro_b(:,1);noise_micro_d(:,1)];

ratio_obs=noise_obs./var_obs;

ratio=zeros(9,9);
for i=1:9
    for j=1:9
        if i==j
            ratio(i,j)=ratio_obs(i);
        end
    end
end

clear x
clear xgrid
clear p
clear pgrid
x=0:500:1000*sw_dist([c2_pos(2) d_pos(2)],[c2_pos(1) d_pos(1)],'km'); %interpolate/extrapolate from the coast out to mooring B
p=10:10:4200; %just extrapolate down to 4200 dbar, then later will make anything under topography into NaN
for i=1:length(p)
    xgrid(i,:)=x(:);
end
for i=1:length(x)
    pgrid(:,i)=p(:);
end

xc_l=130*1000; %horizontal correlation length from cross correlations (m)
zc_l=790; %vertical correlation length from cross correlations (m)

x_corr_func=@(x) exp(-(x(:)/xc_l).^2).*cos(pi.*x(:)./(2.*xc_l)); %horiz correlation function
z_corr_func=@(z) exp(-(z(:)/zc_l).^2);

% find the closest grid point to subtract off 2013 data

for i=1:9
    for j=1:152
        grid_dist(i,j)=abs(x(1,j)-dx_obs(i));
    end
    for j=1:420
        for k=1:722
            dp_dist(i,j,k)=abs(p(1,j)-dp_obs(i,k));
        end
    end
end

for i=1:9
    [M(i),I(i)] = min(grid_dist(i,:));
    [Mp(i,:),Ip(i,:)]=min(dp_dist(i,:,:));
end

for time=1:size(micro_b.pres,1) 
    for k=1:9
        temp_obs_minus13(k,time)=temp_obs(k,time)-int_t_2013(Ip(k,time),I(k))-int_t2_2013(Ip(k,time),I(k)); %not totally sure this is right
    end
end

clear weight_corr
clear cross_corr
clear weights
for time=1:size(micro_b.pres,1) % this loop is where the optimal interpolation of temperature actually happens
    for i=1:length(ratio)
        for j=1:length(pgrid)
            for k=1:size(pgrid,2)
                weight_corr(i,j,k)=x_corr_func(abs(xgrid(j,k)-dx_obs(i)))*z_corr_func(abs(pgrid(j,k)-dp_obs(i)));
            end
        end
    end
    
    for i=1:length(ratio)
        for j=1:length(ratio)
            cross_corr(i,j)=x_corr_func(abs(dx_obs(i)-dx_obs(j)))*z_corr_func(abs(dp_obs(i)-dp_obs(j)));
        end
    end

    for j=1:length(pgrid)
        for k=1:size(pgrid,2)
            weights(:,j,k)=(ratio+cross_corr)\weight_corr(:,j,k);
        end
    end
    
    for j=1:length(pgrid)
        for k=1:size(pgrid,2)
            int_t(j,k,time)=weights(:,j,k).'*temp_obs_minus13(:,time); %check to see if dimensions work?
        end
    end
end

% subtract and do it again on small scale

for time=1:size(micro_b.pres,1)
    for k=1:9
        temp_obs_anom(k,time)=temp_obs_minus13(k,time)-int_t(Ip(k,time),I(k),time); %not totally sure this is right
    end
end

xc=130*1000*.2; %horizontal correlation length from cross correlations (m)
zc=790*.2; %vertical correlation length from cross correlations (m)

x_corr_func=@(x) exp(-(x(:)/xc).^2).*cos(pi.*x(:)./(2.*xc)); %horiz correlation function
z_corr_func=@(z) exp(-(z(:)/zc).^2);
ratio=ratio*.2;

clear weight_corr
clear cross_corr
clear weights
for time=1:size(micro_b.pres,1) % this loop is where the optimal interpolation of temperature actually happens
    for i=1:length(ratio)
        for j=1:length(pgrid)
            for k=1:size(pgrid,2)
                weight_corr(i,j,k)=x_corr_func(abs(xgrid(j,k)-dx_obs(i)))*z_corr_func(abs(pgrid(j,k)-dp_obs(i)));
            end
        end
    end
    
    for i=1:length(ratio)
        for j=1:length(ratio)
            cross_corr(i,j)=x_corr_func(abs(dx_obs(i)-dx_obs(j)))*z_corr_func(abs(dp_obs(i)-dp_obs(j)));
        end
    end

    for j=1:length(pgrid)
        for k=1:size(pgrid,2)
            weights(:,j,k)=(ratio+cross_corr)\weight_corr(:,j,k);
        end
    end
    
    for j=1:length(pgrid)
        for k=1:size(pgrid,2)
            int_t2(j,k,time)=weights(:,j,k).'*temp_obs_anom(:,time); %check to see if dimensions work?
        end
    end
end

% calculate the expected error
background_error=nan(420,152);
for i=1:size(pgrid,1)
    for j=1:size(pgrid,2)
        background_error(i,j)=1; %let's say it's .01 C everywhere, see what happens
    end
end

clear analysis_error
for i=1:size(pgrid,1)
    for j=1:size(pgrid,2)
        %term2(i,j)=(weight_corr(:,i,j).'*inv(ratio+cross_corr)*weight_corr(:,i,j));
        analysis_error(i,j)=background_error(i,j)-(weight_corr(:,i,j).'*inv(ratio+cross_corr)*weight_corr(:,i,j));
    end
end

% graph up the time mean
for i=1:420
    for j=1:152
        mean_t(i,j)=nanmean(int_t(i,j,:));
        mean_t2(i,j)=nanmean(int_t2(i,j,:));
    end
end

% delete parts under topography

figure
hold on
[C,h]=contourf(xgrid,pgrid,int_t_2013(:,1:152)+int_t2_2013(:,1:152),0:1:20)
axis 'ij'
xlabel('Distance from coast (m)','FontSize',24)
ylabel('Pressure (dbar)','FontSize',24)
title('Background temperature','FontSize',24)
caxis([0 20])
clabel(C,h,0:5:20)
cmocean('thermal')
colorbar
brighten(.3)
scatter(dx_obs,dp_obs(:,1),'k','filled')

figure
hold on
[C,h]=contourf(xgrid,pgrid,mean_t,-5:.1:5)
axis 'ij'
xlabel('Distance from coast (m)','FontSize',24)
ylabel('Pressure (dbar)','FontSize',24)
title('Time mean interpolated temperature: large scale','FontSize',24)
%caxis
clabel(C,h,-5:.2:5)
cmocean('balance','zero')
colorbar

figure
hold on
[C,h]=contourf(xgrid,pgrid,mean_t2,-5:.1:5)
axis 'ij'
xlabel('Distance from coast (m)','FontSize',24)
ylabel('Pressure (dbar)','FontSize',24)
title('Time mean interpolated temperature: small scale','FontSize',24)
clabel(C,h,-5:.2:5)
cmocean('balance','zero')
colorbar


figure
hold on
[C,h]=contourf(xgrid,pgrid,mean_t+mean_t2+int_t_2013(:,1:152)+int_t2_2013(:,1:152),0:1:20)
axis 'ij'
xlabel('Distance from coast (m)','FontSize',24)
ylabel('Pressure (dbar)','FontSize',24)
title('Time mean interpolated temperature: large+small scale+background','FontSize',24)
caxis([0 20])
clabel(C,h,0:5:20)
cmocean('thermal')
colorbar
brighten(.3)

figure
hold on
[C,h]=contourf(xgrid,pgrid,analysis_error,0:.125:1)
axis 'ij'
xlabel('Distance from coast (m)','FontSize',24)
ylabel('Pressure (dbar)','FontSize',24)
title('Expected error ratio','FontSize',24)
clabel(C,h,0:.25:1)
cmocean('curl','zero')
colorbar
brighten(.3)

%% interpolate salinity in a similar fashion
xc=77*1000; 
zc=171; 

x_corr_func=@(x) exp(-(x(:)/xc).^2).*cos(pi.*x(:)./(2.*xc)); %horiz correlation function
z_corr_func=@(z) exp(-(z(:)/zc).^2);

for i=1:9
    for j=1:152
        grid_dist(i,j)=abs(x(1,j)-dx_obs(i));
    end
    for j=1:420
        for k=1:722
            dp_dist(i,j,k)=abs(p(1,j)-dp_obs(i,k));
        end
    end
end

for i=1:9
    [M(i),I(i)] = min(grid_dist(i,:));
    [Mp(i,:),Ip(i,:)]=min(dp_dist(i,:,:));
end

B_sal_var=nanvar(B_int.sal); %variance of the T measurements, to create a ratio with the noise level
D_sal_var=nanvar(D_int.sal);

%dx_obs=[B_dx;B_dx;B_dx;B_dx;B_dx;D_dx;D_dx;D_dx;D_dx];
for i=1:size(micro_b.temp,1)
    dp_obs(:,i)=[B_int.pres(i,:).';D_int.pres(i,:).'];
    sal_obs(:,i)=[B_int.sal(i,:).';D_int.sal(i,:).'];
end

for time=1:size(micro_b.pres,1) 
    for k=1:9
        sal_obs_minus13(k,time)=sal_obs(k,time)-int_s_2013(Ip(k,time),I(k))-int_s2_2013(Ip(k,time),I(k)); %not totally sure this is right
    end
end

var_obs=[B_sal_var.';D_sal_var.'];
noise_obs=[noise_micro_b(:,4);noise_micro_d(:,4)];

ratio_obs=noise_obs./var_obs;

ratio=zeros(9,9);
for i=1:9
    for j=1:9
        if i==j
            ratio(i,j)=ratio_obs(i);
        end
    end
end

clear weight_corr
clear cross_corr
clear weights
for time=1:size(micro_b.pres,1) % this loop is where the optimal interpolation of salinity actually happens
    for i=1:length(ratio)
        for j=1:length(pgrid)
            for k=1:size(pgrid,2)
                weight_corr(i,j,k)=x_corr_func(abs(xgrid(j,k)-dx_obs(i)))*z_corr_func(abs(pgrid(j,k)-dp_obs(i)));
            end
        end
    end
    
    for i=1:length(ratio)
        for j=1:length(ratio)
            cross_corr(i,j)=x_corr_func(abs(dx_obs(i)-dx_obs(j)))*z_corr_func(abs(dp_obs(i)-dp_obs(j)));
        end
    end

    for j=1:length(pgrid)
        for k=1:size(pgrid,2)
            weights(:,j,k)=(ratio+cross_corr)\weight_corr(:,j,k);
        end
    end
    
    for j=1:length(pgrid)
        for k=1:size(pgrid,2)
            int_s(j,k,time)=weights(:,j,k).'*sal_obs_minus13(:,time); %check to see if dimensions work?
        end
    end
end

% subtract, do smaller scale
xc=xc*.2;
zc=zc*.2
ratio=ratio*.2;

for time=1:size(micro_b.pres,1)
    for k=1:9
        sal_obs_anom(i,time)=sal_obs_minus13(k,time)-int_s(Ip(k),I(k),time); %not totally sure this is right
    end
end

clear weight_corr
clear cross_corr
clear weights
for time=1:size(micro_b.pres,1) % this loop is where the optimal interpolation of salinity actually happens
    for i=1:length(ratio)
        for j=1:length(pgrid)
            for k=1:size(pgrid,2)
                weight_corr(i,j,k)=x_corr_func(abs(xgrid(j,k)-dx_obs(i)))*z_corr_func(abs(pgrid(j,k)-dp_obs(i)));
            end
        end
    end
    
    for i=1:length(ratio)
        for j=1:length(ratio)
            cross_corr(i,j)=x_corr_func(abs(dx_obs(i)-dx_obs(j)))*z_corr_func(abs(dp_obs(i)-dp_obs(j)));
        end
    end

    for j=1:length(pgrid)
        for k=1:size(pgrid,2)
            weights(:,j,k)=(ratio+cross_corr)\weight_corr(:,j,k);
        end
    end
    
    for j=1:length(pgrid)
        for k=1:size(pgrid,2)
            int_s2(j,k,time)=weights(:,j,k).'*sal_obs_anom(:,time); %check to see if dimensions work?
        end
    end
end


for i=1:420
    for j=1:152
        mean_s(i,j)=nanmean(int_s(i,j,:));
        mean_s2(i,j)=nanmean(int_s2(i,j,:));
    end
end

% expected error
background_error=nan(420,152);
for i=1:size(pgrid,1)
    for j=1:size(pgrid,2)
        background_error(i,j)=1; %let's say it's .01 C everywhere, see what happens
    end
end

clear analysis_error
for i=1:size(pgrid,1)
    for j=1:size(pgrid,2)
        %term2(i,j)=(weight_corr(:,i,j).'*inv(ratio+cross_corr)*weight_corr(:,i,j));
        analysis_error(i,j)=background_error(i,j)-(weight_corr(:,i,j).'*inv(ratio+cross_corr)*weight_corr(:,i,j));
    end
end

% figures
figure
hold on
[C,h]=contourf(xgrid,pgrid,int_s_2013(:,1:152)+int_s2_2013(:,1:152),0:.25:40)
axis 'ij'
xlabel('Distance from coast','FontSize',24)
ylabel('Pressure (dbar)','FontSize',24)
title('Background salinity','FontSize',24)
caxis([34 36])
clabel(C,h,0:.5:40)
cmocean('haline')
colorbar

figure
hold on
[C,h]=contourf(xgrid,pgrid,mean_s,-3:.01:3)
axis 'ij'
xlabel('Distance from coast (m)','FontSize',24)
ylabel('Pressure (dbar)','FontSize',24)
title('Large scale salinity','FontSize',24)
clabel(C,h,-3:.02:3)
%caxis
cmocean('balance','zero')
colorbar

% figure
% hold on
% contourf(xgrid,pgrid,mean_s2)
% axis 'ij'
% xlabel('Distance from coast (m)','FontSize',24)
% ylabel('Pressure (dbar)','FontSize',24)
% title('Small scale salinity')
% %clabel
% %caxis
% cmocean('balance','zero')
% colorbar


figure
hold on
[C,h]=contourf(xgrid,pgrid,mean_s2+mean_s+int_s_2013(:,1:152)+int_s2_2013(:,1:152),0:.25:40)
axis 'ij'
xlabel('Distance from coast (m)','FontSize',24)
ylabel('Pressure (dbar)','FontSize',24)
title('Background+large+small scale salinity','FontSize',24)
clabel(C,h,0:.5:40)
caxis([34 36])
cmocean('haline')
colorbar

figure
hold on
[C,h]=contourf(xgrid,pgrid,analysis_error,0:.125:1)
axis 'ij'
xlabel('Distance from coast (m)','FontSize',24)
ylabel('Pressure (dbar)','FontSize',24)
title('Expected error ratio','FontSize',24)
clabel(C,h,0:.25:1)
cmocean('curl','zero')
colorbar
brighten(.3)


%% make parts under topography NaN

% FIGURE THIS OUT.
overall_lat = interp1([0,x(152)],[instr_pos(1,2),instr_pos(6,2)],[0:500:x(152)]); 
overall_lon = interp1([0,x(152)],[instr_pos(1,1),instr_pos(6,1)],[0:500:x(152)]);

% now compare our interpolated points to ETOPO2 and find the closest match
% to each point
% will want to use high res bottom topography for real thing but this will
% give a general idea - issue - this is depth not pressure
[elev,long,lat]=m_etopo2([27 28 -34 -33]);

dist = nan(size(overall_lat,2),size(lat,1),size(long,1));
for i=1:size(overall_lat,2)
    for j=1:size(lat,1)
        for k=1:size(long,1)
            dist(i,j,k) = sw_dist([overall_lat(1,i) lat(j,k)],[overall_lon(1,i) long(j,k)],'km');
        end 
    end
end

clear M
clear I
closest = nan(182,1);
for i=1:size(overall_lat,2)
    min_temp = dist(i,:,:);
    temp = min_temp(:);
    [M,I] = min(temp);
    closest(i,1) = I;
end

for i=1:size(overall_lat,2)
    [I_row(i),I_col(i)] = ind2sub([61,61],closest(i));
end

topo = nan(182,1);
for i=1:size(overall_lat,2)
    topo(i,1) = elev(I_row(i),I_col(i));
    topo_p(i,1)=sw_pres(topo(i,1),b_pos(2)); % change into pressure coordinates using seawater package - again not the most precise way, will need to change way of defining topography later
end
% turn this into pressure - seawater package???

pres_int = nan(420,size(micro_b.pres,1),152);
for i=1:size(micro_b.pres,1)
    for j=1:182 %indices
        for k=1:420
            pres_int(k,i,j) = 10*k;
        end
    end
end

for i=1:722
    for j=1:152
        for k=1:420
            if pres_int(k,i,j) > abs(topo_p(j,1))
                int_t(k,j,i) = NaN;
                int_t2(k,j,i) = NaN;
                int_s(k,j,i) = NaN;
                int_s2(k,j,i) = NaN;
            end
        end
    end
end


