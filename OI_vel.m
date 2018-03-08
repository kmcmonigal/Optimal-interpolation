% interpolate the ACT velocities at every time step using OI
% first and second deployments done seperately

test_act=1;
timesteps=1;

if test_act
    load('act1-A_cm_velocities.mat')
    load('act1-B_cm_velocities.mat')
    load('act1-C_cm_velocities.mat')
    load('act1-D_cm_velocities.mat')
    load('act1-E_cm_velocities.mat')
    load('act1-F_cm_velocities.mat')
    load('act1-G_cm_velocities.mat')
    load('ACT_CPIES_velocities_20h.mat')
    load('P2_241_currents.mat')
    
    % make all the records start and same time, get rid of gaps
    A.u2=nan(50,13645-14-10);
    A.v2=nan(50,13645-14-10);
    A.z2=nan(50,13645-14-10);
    A.u2(:,1:9480)=A.u_adcp(:,186:186+9480-1);
    A.u2(:,9489:10956)=A.u_adcp(:,186+9480:186+10956-9);
    A.u2(:,10965:12439)=A.u_adcp(:,186+10957-9:186+12439-17);
    A.u2(:,12448:13645-14-10)=A.u_adcp(:,186+12440-17:13792-10);
    A.v2(:,1:9480)=A.v_adcp(:,186:186+9480-1);
    A.v2(:,9489:10956)=A.v_adcp(:,186+9480:186+10956-9);
    A.v2(:,10965:12439)=A.v_adcp(:,186+10957-9:186+12439-17);
    A.v2(:,12448:13645-14-10)=A.v_adcp(:,186+12440-17:13792-10);
    A.z2(:,1:9480)=A.z_adcp(:,186:186+9480-1);
    A.z2(:,9489:10956)=A.z_adcp(:,186+9480:186+10956-9);
    A.z2(:,10965:12439)=A.z_adcp(:,186+10957-9:186+12439-17);
    A.z2(:,12448:13645-14-10)=A.z_adcp(:,186+12440-17:13792-10);
    
    B.u2=B.u(:,151+14:13785);
    B.v2=B.v(:,151+14:13795-10);
    B.z2=B.z(:,151+14:13795-10);
    C.u2=C.u(:,129+14:13773-10);
    C.v2=C.v(:,129+14:13773-10);
    C.z2=C.z(:,129+14:13773-10);
    D.u2=D.u(:,78+14:13722-10);
    D.v2=D.v(:,78+14:13722-10);
    D.z2=D.z(:,78+14:13722-10);
    E.u2=E.u(:,58+14:13702-10);
    E.v2=E.v(:,58+14:13702-10);
    E.z2=E.z(:,58+14:13702-10);
    F.u2=F.u(:,30+14:13674-10);
    F.v2=F.v(:,30+14:13674-10);
    F.z2=F.z(:,30+14:13674-10);
    G.u2=G.u(:,1+14:13645-10);
    G.v2=G.v(:,1+14:13645-10);
    G.z2=G.z(:,1+14:13645-10);
    
    hr=20; % do 20 hour subsampling
    A.u=A.u2(:,1:hr:end);
    B.u=B.u2(:,1:hr:end);
    C.u=C.u2(:,1:hr:end);
    D.u=D.u2(:,1:hr:end);
    E.u=E.u2(:,1:hr:end);
    F.u=F.u2(:,1:hr:end);
    G.u=G.u2(:,1:hr:end);
    A.v=A.v2(:,1:hr:end);
    B.v=B.v2(:,1:hr:end);
    C.v=C.v2(:,1:hr:end);
    D.v=D.v2(:,1:hr:end);
    E.v=E.v2(:,1:hr:end);
    F.v=F.v2(:,1:hr:end);
    G.v=G.v2(:,1:hr:end);
    A.z=A.z2(:,1:hr:end);
    B.z=B.z2(:,1:hr:end);
    C.z=C.z2(:,1:hr:end);
    D.z=D.z2(:,1:hr:end);
    E.z=E.z2(:,1:hr:end);
    F.z=F.z2(:,1:hr:end);
    G.z=G.z2(:,1:hr:end);
    
    timesteps=10;
end

% make some choices
use_background=0; % subtract a "background field" and map anomalies vs mapping field itself
Lx=63*1000; % horizontal decorrelation length scale
Lz=1913; % vertical decorrelation length scale
lx=.3*63*1000; % small scale decorrelation length scales
lz=.3*1913;
N=.07; % noise ratio

% correlation functions
Xc=@(x) exp(-(x(:)/Lx).^2).*cos(pi.*x(:)./(2.*Lx))+exp(-(x(:)/lx).^2).*cos(pi.*x(:)./(2.*lx));
Zc=@(z) exp(-(z(:)/Lz).^2)+exp(-(z(:)/lz).^2);

% instrument locations - will need to edit to exact locations
a_pos = [27.595,-33.5583];
b_pos = [27.6428,-33.6674];
c_pos = [27.7152,-33.7996];
d_pos = [27.8603,-34.0435];
e_pos = [28.0316,-34.29];
f_pos = [28.17,-34.5380];
g_pos = [28.3453,-34.8213];

% map with position of mooring a as left bound, mooring g as right bound
B_dx=1000*sw_dist([a_pos(2) b_pos(2)],[a_pos(1) b_pos(1)],'km');
C_dx=1000*sw_dist([a_pos(2) c_pos(2)],[a_pos(1) c_pos(1)],'km');
D_dx=1000*sw_dist([a_pos(2) d_pos(2)],[a_pos(1) d_pos(1)],'km');
E_dx=1000*sw_dist([a_pos(2) e_pos(2)],[a_pos(1) e_pos(1)],'km');
F_dx=1000*sw_dist([a_pos(2) f_pos(2)],[a_pos(1) f_pos(1)],'km');
G_dx=1000*sw_dist([a_pos(2) g_pos(2)],[a_pos(1) g_pos(1)],'km');

% make a matrix of the distance of each observation from A, will need to
% fill in
dist=zeros(1,size(A.u,1)+size(B.u,1)+size(C.u,1)+size(D.u,1)+size(E.u,1)+size(F.u,1)+size(G.u,1));
for i=1:size(A.u,1)
    dist(i)=0;
end
for i=size(A.u,1)+1:size(A.u,1)+size(B.u,1)
    dist(i)=B_dx;
end
for i=size(A.u,1)+size(B.u,1)+1:size(A.u,1)+size(B.u,1)+size(C.u,1)
    dist(i)=C_dx;
end
for i=size(A.u,1)+size(B.u,1)+size(C.u,1)+1:size(A.u,1)+size(B.u,1)+size(C.u,1)+size(D.u,1)
    dist(i)=D_dx;
end
for i=size(A.u,1)+size(B.u,1)+size(C.u,1)+size(D.u,1)+1:size(A.u,1)+size(B.u,1)+size(C.u,1)+size(D.u,1)+size(E.u,1)
    dist(i)=E_dx;
end
for i=size(A.u,1)+size(B.u,1)+size(C.u,1)+size(D.u,1)+size(E.u,1):size(A.u,1)+size(B.u,1)+size(C.u,1)+size(D.u,1)+size(E.u,1)+size(F.u,1)
    dist(i)=F_dx;
end
for i=size(A.u,1)+size(B.u,1)+size(C.u,1)+size(D.u,1)+size(E.u,1)+size(F.u,1):size(A.u,1)+size(B.u,1)+size(C.u,1)+size(D.u,1)+size(E.u,1)+size(F.u,1)+size(G.u,1)
    dist(i)=G_dx;
end

% make a grid to map data onto
x=0:500:G_dx;
z=0:20:4500; % then later will make anything under topography nan

[xgrid,zgrid]=meshgrid(x,z);

% set up a noise matrix by finding variace of each measurement
% need to go back and review why exactly i was doing it this way, what is
% best way to quantify noise

noisev=N./nanvar([A.v;B.v;C.v;D.v;E.v;F.v;G.v].'); % right the problem with this is that number of measurements changes at each time point
noiseu=N./nanvar([A.u;B.u;C.u;D.u;E.u;F.v;G.u].');

% map time mean, do two iterations
noise=complex(diag(noiseu),diag(noisev));
noise(isnan(noise))=complex(.01,.01);
z_mean=[nanmean(A.z,2);nanmean(B.z,2);nanmean(C.z,2);nanmean(D.z,2);nanmean(E.z,2);nanmean(F.z,2);nanmean(G.z,2)];

weight_corr=nan(length(noise),size(zgrid,1),size(zgrid,2));
for j=1:length(noise)
    for k=1:size(zgrid,1) % was length(zgrid)
        for l=1:size(zgrid,2)
            weight_corr(j,k,l)=Xc(abs(xgrid(k,l)-dist(j)))*Zc(abs(zgrid(k,l)-z_mean(j)));
        end
    end
end

% get cross corr between instruments
cross_corr=nan(length(noise),length(noise));
for j=1:length(noise)
    for k=1:length(noise)
        cross_corr(j,k)=Xc(abs(dist(j)-dist(k)))*Zc(abs(z_mean(j)-z_mean(k)));
    end
end

% solve for weights
weights=nan(size(weight_corr,1),size(zgrid,1),size(zgrid,2));
for j=1:size(zgrid,1)
    for k=1:size(zgrid,2)
        weights(:,j,k)=((noise+cross_corr)\weight_corr(:,j,k));
    end
end

u=[nanmean(A.u,2);nanmean(B.u,2);nanmean(C.u,2);nanmean(D.u,2);nanmean(E.u,2);nanmean(F.u,2);nanmean(G.u,2)]; % don't get along track velocity from cpies. let's assume not using cpies
v=[nanmean(A.v,2);nanmean(B.v,2);nanmean(C.v,2);nanmean(D.v,2);nanmean(E.v,2);nanmean(F.v,2);nanmean(G.v,2)];% cpies34(:),i);cpies45(:,i)];

mean_u=nanmean(nanmean(u));
mean_v=nanmean(nanmean(v));

w=complex(u,v)-complex(mean_u,mean_v);
weights(isnan(w))=0;
w(isnan(w))=0;
for j=1:size(zgrid,1)
    for k=1:size(zgrid,2)
        vel_mean(j,k)=weights(:,j,k).'*w+complex(mean_u,mean_v); % if w has nan's it gives us all nan. need to fix that. do weights have to add to something to get a right answer???
    end
end

% now map anomalies relative to time mean
tic
for i=1:timesteps
    % put all the obervations from that point in time into a matrix
    u=[A.u(:,i);B.u(:,i);C.u(:,i);D.u(:,i);E.u(:,i);F.u(:,i);G.u(:,i)];
    v=[A.v(:,i);B.v(:,i);C.v(:,i);D.v(:,i);E.v(:,i);F.v(:,i);G.v(:,i)];
    
    % depth of each of those observations
    z=[A.z(:,i);B.z(:,i);C.z(:,i);D.z(:,i);E.z(:,i);F.z(:,i);G.z(:,i)];
    
    % see if any of the observations are nan, get rid of those
    gaps=[isnan(A.u(:,i));isnan(B.u(:,i));isnan(C.u(:,i));isnan(D.u(:,i));isnan(E.u(:,i));isnan(F.u(:,i));isnan(G.u(:,i))];
    
    u2=u(gaps==0);
    v2=v(gaps==0);
    z2=z(gaps==0);
    dist2=dist(gaps==0);
    noisev2=noisev(gaps==0);
    noiseu2=noiseu(gaps==0);

    % make diagonal noise matrix
    noise=complex(diag(noiseu2),diag(noisev2));
    
    % let's map as complex data: w=u+iv
    w=complex(u2,v2);
    
    for j=1:size(z2,1)
        for k=1:size(zgrid,1)
            dz_dist(j,k)=abs(zgrid(k,1)-z2(j));
        end
        [Mz,Iz] = min(dz_dist(j,:));
        for k=1:size(xgrid,2)
            grid_dist(j,k)=abs(xgrid(1,k)-dist2(j));
        end
        [M,I] = min(grid_dist(j,:));
        vel_mean_anom(j)=w(j)-vel_mean(Iz,I);
    end
    
    vel_mean_anom2=vel_mean_anom(gaps==0);
    
    % now get cross correlations between instruments, grid points
    weight_corr=nan(length(noise),length(zgrid),size(zgrid,2));
    for j=1:length(noise)
        for k=1:size(zgrid,1) % was length(zgrid)
            for l=1:size(zgrid,2)
                weight_corr(j,k,l)=Xc(abs(xgrid(k,l)-dist2(j)))*Zc(abs(zgrid(k,l)-z2(j)));
            end
        end
    end
    
    % get cross corr between instruments
    cross_corr=nan(length(noise),length(noise));
    for j=1:length(noise)
        for k=1:length(noise)
            cross_corr(j,k)=Xc(abs(dist2(j)-dist2(k)))*Zc(abs(z2(j)-z2(k)));
        end
    end
    
    % solve for weights
    weights=nan(size(weight_corr,1),size(zgrid,1),size(zgrid,2));
    for j=1:size(zgrid,1)
        for k=1:size(zgrid,2)
            weights(:,j,k)=((noise+cross_corr)\weight_corr(:,j,k)); % what are ratio_ac, ac_obs
            temp(j,k)=weights(:,j,k).'*vel_mean_anom2.';
        end
    end
    vel_total(:,:,i)=temp+vel_mean;

    disp([i toc])
end
