% optimally interpolate T and S using similar method as OI_vel

load('ACTGEMData_skype.mat');
load('micro_b_d_struct_despike_filter.mat');
load('micro_b_noise.mat');
load('micro_d_noise.mat');

coast_lat=-33.2910;
coast_lon=27.4783;

% create a climatology of T, S to subtract as the "large scale" term (ex:
% Roemmich 1983)
ctdg.months = month(datetime(ctdg.datenum,'ConvertFrom','datenum'));

for k=1:size(ctdg.temp,1)
    for j=1:size(ctdg.temp,2)
        for i=1:12
            if ctdg.months(1,j) == i
                monthly_clim_t(k,j,i) = ctdg.temp(k,j);
                monthly_clim_s(k,j,i) = ctdg.sal(k,j);
            end
        end
    end
end

monthly_clim_t(monthly_clim_t==0)=NaN;
monthly_clim_s(monthly_clim_s==0)=NaN;

for i=1:12
    for j=1:size(ctdg.temp,1)
        ave_temp(j,i) = nanmean(monthly_clim_t(j,:,i));
        number_measures(j,i) = sum(~isnan(monthly_clim_t(j,:,i)));
        ave_sal(j,i) = nanmean(monthly_clim_s(j,:,i));
    end
end

% Below 2000m, make ave_temp(:,i) be the same for each month since there
% aren't many observations
for i=1:12
    for j=1:200 %surface to 1990 m
        ave_temp2(j,i) = ave_temp(j,i);
        ave_sal2(j,i) = ave_sal(j,i);
    end
    for j=201:501 %2000 m to bottom
        ave_temp2(j,i) = nanmean(ave_temp(j,:));
        ave_sal2(j,i) = nanmean(ave_sal(j,:));
    end
end

% climatology in pressure coordinates - use microcats in p coords

% now subtract the climatology and do optimal interpolation on the
% anomalies

b_pos = [27.6428,-33.6674]; % position of mooring B
d_pos = [27.8603,-34.0435]; % position of mooring D

B_dx=1000*sw_dist([coast_lat b_pos(2)],[coast_lon b_pos(1)],'km'); %distance from coast to B
D_dx=1000*sw_dist([coast_lat d_pos(2)],[coast_lon d_pos(1)],'km'); %distance from coast to D

xc=130*1000; %horizontal correlation length from cross correlations (m)
zc=790; %vertical correlation length from cross correlations (m)

x_corr_func=@(x) exp(-(x(:)/xc).^2).*cos(pi.*x(:)./(2.*xc)); %horiz correlation function
z_corr_func=@(z) exp(-(z(:)/zc).^2); %vertical correlation function

B_int.temp=micro_b.temp(:,:);
D_int.temp=micro_d.temp(6:727,:); % to make the microcats cover the same time period 

B_int.pres=micro_b.pres(:,:);
D_int.pres=micro_d.pres(6:727,:);

%this section is to subtract a weighted clim mean of the 2 closest 10 dbar bin
%find the 2 closest 10 dbar levels and how far away they are
for i=1:size(micro_b.pres,1)
    for j=1:size(micro_b.pres,2)
        closest1(i,j)=round(B_int.pres(i,j),-1);
        dp1(i,j)=abs(B_int.pres(i,j)-closest1(i,j));
        weight1(i,j)=((dp1(i,j)*dp2(i,j)))/(((dp1(i,j)+dp2(i,j)))*dp1(i,j));
        if B_int.pres(i,j)<closest1(i,j) %rounded value deeper than original
            closest2(i,j)=closest1(i,j)-10;
        else %rounded value shallower than original
            closest2(i,j)=closest1(i,j)+10;
        end
        dp2(i,j)=abs(B_int.pres(i,j)-closest2(i,j));
        weight2(i,j)=((dp1(i,j)*dp2(i,j)))/(((dp1(i,j)+dp2(i,j)))*dp2(i,j));
    end
end

B_int.tempanom=nan(size(micro_b.pres,1),size(micro_b.pres,2));
B_int.months=month(micro_b.date);
for i=1:size(micro_b.pres,1)
    for j=1:size(micro_b.pres,2)
        for k=1:12
            if B_int.months(i,j)==k
                B_int.tempanom(i,j)=B_int.temp(i,j)-(weight1(i,j)*ave_temp2((closest1(i,j)/10)+1,k)+weight2(i,j)*ave_temp2((closest2(i,j)/10)+1,k)); %subtract weighted mean of closest pres bins
            end
        end
    end
end

% same for D - subtract clim temp by finding the 2 closest 10 dbar bins and
% weighting them by distance from instrument
for i=1:size(micro_b.pres,1)
    for j=1:size(micro_d.pres,2)
        closest1(i,j)=round(D_int.pres(i,j),-1);
        dp1(i,j)=abs(D_int.pres(i,j)-closest1(i,j));
        weight1(i,j)=((dp1(i,j)*dp2(i,j)))/(((dp1(i,j)+dp2(i,j)))*dp1(i,j));
        if D_int.pres(i,j)<closest1(i,j) %rounded value deeper than original
            closest2(i,j)=closest1(i,j)-10;
        else %rounded value shallower than original
            closest2(i,j)=closest1(i,j)+10;
        end
        dp2(i,j)=abs(D_int.pres(i,j)-closest2(i,j));
        weight2(i,j)=((dp1(i,j)*dp2(i,j)))/(((dp1(i,j)+dp2(i,j)))*dp2(i,j));
    end
end

D_int.tempanom=nan(size(micro_b.pres,1),size(micro_d.pres,2));
D_int.months=month(micro_b.date);
for i=1:size(micro_b.pres,1)
    for j=1:size(micro_d.pres,2)
        for k=1:12
            if D_int.months(i,j)==k
                D_int.tempanom(i,j)=D_int.temp(i,j)-(weight1(i,j)*ave_temp2((closest1(i,j)/10)+1,k)+weight2(i,j)*ave_temp2((closest2(i,j)/10)+1,k)); %subtract weighted mean of closest pres bins
            end
        end
    end
end

B_temp_var=nanvar(B_int.temp); %variance of the T measurements, to create a ratio with the noise level
D_temp_var=nanvar(D_int.temp);

dx_obs=[B_dx;B_dx;B_dx;B_dx;B_dx;D_dx;D_dx;D_dx;D_dx];
for i=1:size(micro_b.temp,1)
    dp_obs(:,i)=[B_int.pres(i,:).';D_int.pres(i,:).'];
    temp_obs(:,i)=[B_int.tempanom(i,:).';D_int.tempanom(i,:).'];
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


x=0:500:1000*sw_dist([coast_lat d_pos(2)],[coast_lon d_pos(1)],'km'); %interpolate/extrapolate from the coast out to mooring B
p=0:20:5000; %just extrapolate down to 5000 dbar, then later will make anything under topography into NaN
for i=1:length(p)
    xgrid(i,:)=x(:);
end
for i=1:length(x)
    pgrid(:,i)=p(:);
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
            anom_value_t(j,k,time)=weights(:,j,k).'*temp_obs(:,time); %check to see if dimensions work?
        end
    end
end

% interpolate salinity in a similar fashion
xc=77*1000; 
zc=171; 

B_int.sal=micro_b.sal(:,:);
D_int.sal=micro_d.sal(6:727,:);

for i=1:size(micro_b.pres,1)
    for j=1:size(micro_b.pres,2)
        closest1(i,j)=round(B_int.pres(i,j),-1);
        dp1(i,j)=abs(B_int.pres(i,j)-closest1(i,j));
        weight1(i,j)=((dp1(i,j)*dp2(i,j)))/(((dp1(i,j)+dp2(i,j)))*dp1(i,j));
        if B_int.pres(i,j)<closest1(i,j) %rounded value deeper than original
            closest2(i,j)=closest1(i,j)-10;
        else %rounded value shallower than original
            closest2(i,j)=closest1(i,j)+10;
        end
        dp2(i,j)=abs(B_int.pres(i,j)-closest2(i,j));
        weight2(i,j)=((dp1(i,j)*dp2(i,j)))/(((dp1(i,j)+dp2(i,j)))*dp2(i,j));
    end
end

B_int.tempanom=nan(size(micro_b.pres,1),size(micro_b.pres,2));
B_int.months=month(micro_b.date);
for i=1:size(micro_b.pres,1)
    for j=1:size(micro_b.pres,2)
        for k=1:12
            if B_int.months(i,j)==k
                B_int.salanom(i,j)=B_int.sal(i,j)-(weight1(i,j)*ave_sal2((closest1(i,j)/10)+1,k)+weight2(i,j)*ave_sal2((closest2(i,j)/10)+1,k)); %subtract weighted mean of closest pres bins
            end
        end
    end
end

% same for D - subtract clim temp by finding the 2 closest 10 dbar bins and
% weighting them by distance from instrument
for i=1:size(micro_b.pres,1)
    for j=1:size(micro_d.pres,2)
        closest1(i,j)=round(D_int.pres(i,j),-1);
        dp1(i,j)=abs(D_int.pres(i,j)-closest1(i,j));
        weight1(i,j)=((dp1(i,j)*dp2(i,j)))/(((dp1(i,j)+dp2(i,j)))*dp1(i,j));
        if D_int.pres(i,j)<closest1(i,j) %rounded value deeper than original
            closest2(i,j)=closest1(i,j)-10;
        else %rounded value shallower than original
            closest2(i,j)=closest1(i,j)+10;
        end
        dp2(i,j)=abs(D_int.pres(i,j)-closest2(i,j));
        weight2(i,j)=((dp1(i,j)*dp2(i,j)))/(((dp1(i,j)+dp2(i,j)))*dp2(i,j));
    end
end

D_int.tempanom=nan(size(micro_b.pres,1),size(micro_d.pres,2));
D_int.months=month(micro_b.date);
for i=1:size(micro_b.pres,1)
    for j=1:size(micro_d.pres,2)
        for k=1:12
            if D_int.months(i,j)==k
                D_int.salanom(i,j)=D_int.sal(i,j)-(weight1(i,j)*ave_sal2((closest1(i,j)/10)+1,k)+weight2(i,j)*ave_sal2((closest2(i,j)/10)+1,k)); %subtract weighted mean of closest pres bins
            end
        end
    end
end

B_sal_var=nanvar(B_int.sal); %variance of the T measurements, to create a ratio with the noise level
D_sal_var=nanvar(D_int.sal);

dx_obs=[B_dx;B_dx;B_dx;B_dx;B_dx;D_dx;D_dx;D_dx;D_dx];
for i=1:size(micro_b.temp,1)
    dp_obs(:,i)=[B_int.pres(i,:).';D_int.pres(i,:).'];
    sal_obs(:,i)=[B_int.salanom(i,:).';D_int.salanom(i,:).'];
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
            anom_value_s(j,k,time)=weights(:,j,k).'*sal_obs(:,time); %check to see if dimensions work?
        end
    end
end

% add the anomaly values to a climatology

micro_months=month(micro_b.date);

for i=1:size(anom_value_t,1) %make sure indices work out
    for j=1:size(anom_value_t,2)
        for k=1:size(micro_b.pres,1)
            for l=1:12
                if micro_months(k,1)==l
                    int_t(i,j,k)=anom_value_t(i,j,k)+ave_temp2(i,l);
                    int_s(i,j,k)=anom_value_s(i,j,k)+ave_sal2(i,l);
                end
            end
        end
    end
end

% make parts under topography NaN

overall_lat = interp1([0,D_dx],[coast_lat,d_pos(2)],[0:500:D_dx]); 
overall_lon = interp1([0,D_dx],[coast_lon,d_pos(1)],[0:500:D_dx]);

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

pres_int = nan(251,size(micro_b.pres,1),182);
for i=1:size(micro_b.pres,1)
    for j=1:182 %indices
        for k=1:251
            pres_int(k,i,j) = 20*(k-1);
        end
    end
end

for i=1:722
    for j=1:182
        for k=1:251
            if pres_int(k,i,j) > abs(topo_p(j,1))
                int_t(k,j,i) = NaN;
                int_s(k,j,i) = NaN;
            end
        end
    end
end


