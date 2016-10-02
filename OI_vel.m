% optimal interpolation of u and v data
% should I be subtracting off a "large scale field" or no? seems like
% different authors do it differently. see how each works out...

% vertical correlation can't be isotropic - speed needs to increase as you
% get more shallow
load('ASCA0416currents_filter_sub_rotate.mat');
load('cm_noise.mat');

a_pos = [27.595,-33.5583];
b_pos = [27.6428,-33.6674];
c_pos = [27.7152,-33.7996];
d_pos = [27.8603,-34.0435];

A_int.v=A.vtmp3_filt(:,1:725);
B_int.v=B.vtmp4_filt(:,:);
C_int.v=C.vtmp4_filt(:,2:726);
D_int.v=D.vtmp4_filt(:,3:727);

all_obs=[A_int.v;B_int.v;C_int.v;D_int.v];

A_int.z=A.ztmp3(:,1:725);
B_int.z=B.ztmp4(:,:);
C_int.z=C.ztmp4(:,2:726);
D_int.z=D.ztmp4(:,3:727);

x=0:500:1000*sw_dist([coast_lat d_pos(2)],[coast_lon d_pos(1)],'km');
z=0:20:5000; % then later will make anything under topography into NaN
for i=1:length(z)
    xgrid(i,:)=x(:);
end
for i=1:length(x)
    zgrid(:,i)=z(:);
end

A_dx=1000*sw_dist([coast_lat a_pos(2)],[coast_lon a_pos(1)],'km');
B_dx=1000*sw_dist([coast_lat b_pos(2)],[coast_lon b_pos(1)],'km');
C_dx=1000*sw_dist([coast_lat c_pos(2)],[coast_lon c_pos(1)],'km');
D_dx=1000*sw_dist([coast_lat d_pos(2)],[coast_lon d_pos(1)],'km');

xc=50*1000; %horizontal decorrelation length
zc=2200; %vertical decorrelation length

x_corr_func=@(x) exp(-(x(:)/xc).^2).*cos(pi.*x(:)./(2.*xc));
z_corr_func=@(z) exp(-(z(:)/zc).^2);

% optimal interpolation of v
for time=1%:size(A_int.v,2) % change back after running
    clear Noise2
    clear u_obs_time
    clear dx
    clear weights
    clear cross_corr
    clear ratio
    clear weight_corr
    v_obs_time=[A_int.v(~isnan(A_int.v(:,time)),:);B_int.v(~isnan(B_int.v(:,time)),:);C_int.v(~isnan(C_int.v(:,time)),:);D_int.v(~isnan(D_int.v(:,time)),:)];
    dz_obs=[A_int.z(~isnan(A_int.v(:,time)),time);B_int.z(~isnan(B_int.v(:,time)),time);C_int.z(~isnan(C_int.v(:,time)),time);D_int.z(~isnan(D_int.v(:,time)),time)];
    v_obs=[A_int.v(~isnan(A_int.v(:,time)),time);B_int.v(~isnan(B_int.v(:,time)),time);C_int.v(~isnan(C_int.v(:,time)),time);D_int.v(~isnan(D_int.v(:,time)),time)];
    % setting up Noise2 matrix
    for i=1:sum(~isnan(A_int.v(:,time)))
        Noise2(i)=Noise(1);
    end
    for i=sum(~isnan(A_int.v(:,time)))+1:sum(~isnan(A_int.v(:,time)))+sum(~isnan(B_int.v(:,time)))-1
        Noise2(i)=Noise(2);
    end
    Noise2(sum(~isnan(A_int.v(:,time)))+sum(~isnan(B_int.v(:,time))))=Noise(3);
    for i=sum(~isnan(A_int.v(:,time)))+sum(~isnan(B_int.v(:,time)))+1:sum(~isnan(A_int.v(:,time)))+sum(~isnan(B_int.v(:,time)))+sum(~isnan(C_int.v(:,time)))-3
        Noise2(i)=Noise(4);
    end
    Noise2(sum(~isnan(A_int.v(:,time)))+sum(~isnan(B_int.v(:,time)))+sum(~isnan(C_int.v(:,time)))-2)=Noise(5);
    Noise2(sum(~isnan(A_int.v(:,time)))+sum(~isnan(B_int.v(:,time)))+sum(~isnan(C_int.v(:,time)))-1)=Noise(6);
    Noise2(sum(~isnan(A_int.v(:,time)))+sum(~isnan(B_int.v(:,time)))+sum(~isnan(C_int.v(:,time))))=Noise(7);
    for i=sum(~isnan(A_int.v(:,time)))+sum(~isnan(B_int.v(:,time)))+sum(~isnan(C_int.v(:,time)))+1:sum(~isnan(A_int.v(:,time)))+sum(~isnan(B_int.v(:,time)))+sum(~isnan(C_int.v(:,time)))+sum(~isnan(D_int.v(:,time)))-5
        Noise2(i)=Noise(8);
    end
    Noise2(sum(~isnan(A_int.v(:,time)))+sum(~isnan(B_int.v(:,time)))+sum(~isnan(C_int.v(:,time)))+sum(~isnan(D_int.v(:,time)))-4)=Noise(9);
    Noise2(sum(~isnan(A_int.v(:,time)))+sum(~isnan(B_int.v(:,time)))+sum(~isnan(C_int.v(:,time)))+sum(~isnan(D_int.v(:,time)))-3)=Noise(10);
    Noise2(sum(~isnan(A_int.v(:,time)))+sum(~isnan(B_int.v(:,time)))+sum(~isnan(C_int.v(:,time)))+sum(~isnan(D_int.v(:,time)))-2)=Noise(11);
    Noise2(sum(~isnan(A_int.v(:,time)))+sum(~isnan(B_int.v(:,time)))+sum(~isnan(C_int.v(:,time)))+sum(~isnan(D_int.v(:,time)))-1)=Noise(12);
    Noise2(sum(~isnan(A_int.v(:,time)))+sum(~isnan(B_int.v(:,time)))+sum(~isnan(C_int.v(:,time)))+sum(~isnan(D_int.v(:,time))))=Noise(13);
    % set up ratio matrix
    ratio=zeros(length(Noise2),length(Noise2));
    Noise2=Noise2/100;
    for i=1:length(Noise2)
        for j=1:length(Noise2)
            if i==j
                ratio(i,j)=Noise2(i)/nanvar(v_obs_time(i,:)); % seems bad that these values are so big. is this right? what was it on the one obs case?
            end
        end
    end
    % set up dx matrix
    for i=1:sum(~isnan(A_int.v(:,time)))
        dx(i)=A_dx;
    end
    for i=sum(~isnan(A_int.v(:,time)))+1:sum(~isnan(A_int.v(:,time)))+sum(~isnan(B_int.v(:,time)))
        dx(i)=B_dx;
    end
    for i=sum(~isnan(A_int.v(:,time)))+sum(~isnan(B_int.v(:,time)))+1:sum(~isnan(A_int.v(:,time)))+sum(~isnan(B_int.v(:,time)))+sum(~isnan(C_int.v(:,time)))
        dx(i)=C_dx;
    end
    for i=sum(~isnan(A_int.v(:,time)))+sum(~isnan(B_int.v(:,time)))+sum(~isnan(C_int.v(:,time)))+1:sum(~isnan(A_int.v(:,time)))+sum(~isnan(B_int.v(:,time)))+sum(~isnan(C_int.v(:,time)))+sum(~isnan(D_int.v(:,time)))
        dx(i)=D_dx;
    end
    % now get cross correlations between instruments, grid points
    for i=1:length(Noise2)
        for j=1:length(zgrid)
            for k=1:size(zgrid,2)
                weight_corr(i,j,k)=x_corr_func(abs(xgrid(j,k)-dx(i)))*z_corr_func(abs(zgrid(j,k)-dz_obs(i)));
            end
        end
    end
    % get cross corr between instruments - should this be theoretical or
    % "real" though
    for i=1:length(Noise2)
        for j=1:length(Noise2)
            cross_corr(i,j)=x_corr_func(abs(dx(i)-dx(j)))*z_corr_func(abs(dz_obs(i)-dz_obs(j)));
        end
    end
    % solve for weights
    for j=1:length(zgrid)
        for k=1:size(zgrid,2)
            weights(:,j,k)=(ratio+cross_corr)\weight_corr(:,j,k); % may need to do transform of weight_corr
        end
    end

    for j=1:length(zgrid)
        for k=1:size(zgrid,2)
            anom_value_v(j,k,time)=weights(:,j,k).'*v_obs; % mean already subtracted
        end
    end
    
end

% same process for u

A_int.u=A.utmp3_filt(:,1:725);
B_int.u=B.utmp4_filt(:,:);
C_int.u=C.utmp4_filt(:,2:726);
D_int.u=D.utmp4_filt(:,3:727);

all_obs=[A_int.u;B_int.u;C_int.u;D_int.u];

for time=1:size(A_int.v,2) % change back after running
    clear Noise2
    clear u_obs_time
    clear dx
    clear weights
    clear cross_corr
    clear ratio
    clear weight_corr
    u_obs_time=[A_int.u(~isnan(A_int.u(:,time)),:);B_int.u(~isnan(B_int.u(:,time)),:);C_int.u(~isnan(C_int.u(:,time)),:);D_int.u(~isnan(D_int.u(:,time)),:)];
    dz_obs=[A_int.z(~isnan(A_int.u(:,time)),time);B_int.z(~isnan(B_int.u(:,time)),time);C_int.z(~isnan(C_int.u(:,time)),time);D_int.z(~isnan(D_int.u(:,time)),time)];
    u_obs=[A_int.u(~isnan(A_int.u(:,time)),time);B_int.u(~isnan(B_int.u(:,time)),time);C_int.u(~isnan(C_int.u(:,time)),time);D_int.u(~isnan(D_int.u(:,time)),time)];
    % setting up Noise2 matrix
    for i=1:sum(~isnan(A_int.u(:,time)))
        Noise2(i)=Noise(1);
    end
    for i=sum(~isnan(A_int.u(:,time)))+1:sum(~isnan(A_int.u(:,time)))+sum(~isnan(B_int.u(:,time)))-1
        Noise2(i)=Noise(2);
    end
    Noise2(sum(~isnan(A_int.u(:,time)))+sum(~isnan(B_int.u(:,time))))=Noise(3)
    for i=sum(~isnan(A_int.u(:,time)))+sum(~isnan(B_int.u(:,time)))+1:sum(~isnan(A_int.u(:,time)))+sum(~isnan(B_int.u(:,time)))+sum(~isnan(C_int.u(:,time)))-3
        Noise2(i)=Noise(4);
    end
    Noise2(sum(~isnan(A_int.u(:,time)))+sum(~isnan(B_int.u(:,time)))+sum(~isnan(C_int.u(:,time)))-2)=Noise(5);
    Noise2(sum(~isnan(A_int.u(:,time)))+sum(~isnan(B_int.u(:,time)))+sum(~isnan(C_int.u(:,time)))-1)=Noise(6);
    Noise2(sum(~isnan(A_int.u(:,time)))+sum(~isnan(B_int.u(:,time)))+sum(~isnan(C_int.u(:,time))))=Noise(7);
    for i=sum(~isnan(A_int.u(:,time)))+sum(~isnan(B_int.u(:,time)))+sum(~isnan(C_int.u(:,time)))+1:sum(~isnan(A_int.u(:,time)))+sum(~isnan(B_int.u(:,time)))+sum(~isnan(C_int.u(:,time)))+sum(~isnan(D_int.u(:,time)))-5
        Noise2(i)=Noise(8);
    end
    Noise2(sum(~isnan(A_int.u(:,time)))+sum(~isnan(B_int.u(:,time)))+sum(~isnan(C_int.u(:,time)))+sum(~isnan(D_int.u(:,time)))-4)=Noise(9);
    Noise2(sum(~isnan(A_int.u(:,time)))+sum(~isnan(B_int.u(:,time)))+sum(~isnan(C_int.u(:,time)))+sum(~isnan(D_int.u(:,time)))-3)=Noise(10);
    Noise2(sum(~isnan(A_int.u(:,time)))+sum(~isnan(B_int.u(:,time)))+sum(~isnan(C_int.u(:,time)))+sum(~isnan(D_int.u(:,time)))-2)=Noise(11);
    Noise2(sum(~isnan(A_int.u(:,time)))+sum(~isnan(B_int.u(:,time)))+sum(~isnan(C_int.u(:,time)))+sum(~isnan(D_int.u(:,time)))-1)=Noise(12);
    Noise2(sum(~isnan(A_int.u(:,time)))+sum(~isnan(B_int.u(:,time)))+sum(~isnan(C_int.u(:,time)))+sum(~isnan(D_int.u(:,time))))=Noise(13);
    % set up ratio matrix
    ratio=zeros(length(Noise2),length(Noise2));
    Noise2=Noise2/100;
    for i=1:length(Noise2)
        for j=1:length(Noise2)
            if i==j
                ratio(i,j)=Noise2(i)/nanvar(u_obs_time(i,:)); % seems bad that these values are so big. is this right? what was it on the one obs case?
            end
        end
    end
    % set up dx matrix
    for i=1:sum(~isnan(A_int.u(:,time)))
        dx(i)=A_dx;
    end
    for i=sum(~isnan(A_int.u(:,time)))+1:sum(~isnan(A_int.u(:,time)))+sum(~isnan(B_int.u(:,time)))
        dx(i)=B_dx;
    end
    for i=sum(~isnan(A_int.u(:,time)))+sum(~isnan(B_int.u(:,time)))+1:sum(~isnan(A_int.u(:,time)))+sum(~isnan(B_int.u(:,time)))+sum(~isnan(C_int.u(:,time)))
        dx(i)=C_dx;
    end
    for i=sum(~isnan(A_int.u(:,time)))+sum(~isnan(B_int.u(:,time)))+sum(~isnan(C_int.u(:,time)))+1:sum(~isnan(A_int.u(:,time)))+sum(~isnan(B_int.u(:,time)))+sum(~isnan(C_int.u(:,time)))+sum(~isnan(D_int.u(:,time)))
        dx(i)=D_dx;
    end
    % now get cross correlations between instruments, grid points
    for i=1:length(Noise2)
        for j=1:length(zgrid)
            for k=1:size(zgrid,2)
                weight_corr(i,j,k)=x_corr_func(abs(xgrid(j,k)-dx(i)))*z_corr_func(abs(zgrid(j,k)-dz_obs(i)));
            end
        end
    end
    % get cross corr between instruments - should this be theoretical or
    % "real" though
    for i=1:length(Noise2)
        for j=1:length(Noise2)
            cross_corr(i,j)=x_corr_func(abs(dx(i)-dx(j)))*z_corr_func(abs(dz_obs(i)-dz_obs(j)));
        end
    end
    % solve for weights
    for j=1:length(zgrid)
        for k=1:size(zgrid,2)
            weights(:,j,k)=(ratio+cross_corr)\weight_corr(:,j,k); % may need to do transform of weight_corr
        end
    end

    for j=1:length(zgrid)
        for k=1:size(zgrid,2)
            anom_value_u(j,k,time)=weights(:,j,k).'*u_obs; % mean already subtracted
        end
    end
    
end

% make NaN if below topography

B_dx=1000*sw_dist([coast_lat b_pos(2)],[coast_lon b_pos(1)],'km');
D_dx=1000*sw_dist([coast_lat d_pos(2)],[coast_lon d_pos(1)],'km');

overall_lat = interp1([0,D_dx],[coast_lat,d_pos(2)],[0:500:D_dx]); 
overall_lon = interp1([0,D_dx],[coast_lon,d_pos(1)],[0:500:D_dx]);

% now compare our interpolated points to ETOPO2 and find the closest match
% to each point
[elev,long,lat]=m_etopo2([27 28 -34 -33]);

dist = nan(182,61,61);
for i=1:182
    for j=1:61
        for k=1:61
            dist(i,j,k) = sw_dist([overall_lat(1,i) lat(j,k)],[overall_lon(1,i) long(j,k)],'km');
        end % why overall_lat, overall_lon shorter than my interpolated values
    end
end

clear M
clear I
closest = nan(182,1);
for i=1:182
    min_temp = dist(i,:,:);
    temp = min_temp(:);
    [M,I] = min(temp);
    closest(i,1) = I;
end

for i=1:182
    [I_row(i),I_col(i)] = ind2sub([61,61],closest(i));
end

topo = nan(182,1);
for i=1:182
    topo(i,1) = elev(I_row(i),I_col(i));
end

depth_int = nan(251,725,182);
for i=1:725
    for j=1:182
        for k=1:251
            depth_int(k,i,j) = 20*(k-1);
        end
    end
end

for i=1:725
    for j=1:182
        for k=1:251
            if depth_int(k,i,j) > abs(topo(j,1))
                u_vel(k,j,i) = NaN;
                v_vel(k,j,i) = NaN;
            end
        end
    end
end

% rotate into along, cross stream
[across,along,dist,angle]=rotate_hydro([b_pos(2) d_pos(2)],[b_pos(1) d_pos(1)],u_vel,v_vel);