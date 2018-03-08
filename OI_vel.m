% interpolate the ACT velocities at every time step using OI
% first and second deployments done seperately

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

time1=G.sd(1,1+14:20:13645-10); % all records on this time step

cpies34=v34(:,1:682);
cpies45=v45(:,1:682);

% downscale cpies measurements to every 100 m to speed up computation
cpies45_=cpies45(1:5:end,:);
cpies34_=cpies34(1:5:end,:);
dpth45_=dpth45(1:5:end);
dpth34_=dpth34(1:5:end);

% deep current meter from p2
p2.u=p2.u(1:20:682*20);
p2.v=p2.v(1:20:682*20);

% instrument locations, use to get distances between instruments

a_pos = [27.595,-33.5583];
b_pos = [27.6428,-33.6674];
c_pos = [27.7152,-33.7996];
d_pos = [27.8603,-34.0435];
e_pos = [28.0316,-34.29];
f_pos = [28.17,-34.5380];
g_pos = [28.3453,-34.8213];
cpies34_pos = [(28.3453+28.6248)/2,-(34.8213+35.2453)/2]; % half way between the 2 cpies
cpies45_pos = [(28.6248+28.9)/2,-(35.2453+35.7338)/2];
p2_pos = [p2.lon(1),p2.lat(1)];

x=0:500:1000*sw_dist([a_pos(2) cpies45_pos(2)],[a_pos(1) cpies45_pos(1)],'km');
z=0:20:4500; % then later will make anything under topography into NaN 10:10:4200
for i=1:length(z)
    xgrid(i,:)=x(:);
end
for i=1:length(x)
    zgrid(:,i)=z(:);
end

[A.ac,A.al]=rotate_hydro([a_pos(2) cpies45_pos(2)],[a_pos(1) cpies45_pos(1)],A.u,A.v);
[B.ac,B.al]=rotate_hydro([a_pos(2) cpies45_pos(2)],[a_pos(1) cpies45_pos(1)],B.u,B.v);
[C.ac,C.al]=rotate_hydro([a_pos(2) cpies45_pos(2)],[a_pos(1) cpies45_pos(1)],C.u,C.v);
[D.ac,D.al]=rotate_hydro([a_pos(2) cpies45_pos(2)],[a_pos(1) cpies45_pos(1)],D.u,D.v);
[E.ac,E.al]=rotate_hydro([a_pos(2) cpies45_pos(2)],[a_pos(1) cpies45_pos(1)],E.u,E.v);
[F.ac,F.al]=rotate_hydro([a_pos(2) cpies45_pos(2)],[a_pos(1) cpies45_pos(1)],F.u,F.v);
[G.ac,G.al]=rotate_hydro([a_pos(2) cpies45_pos(2)],[a_pos(1) cpies45_pos(1)],G.u,G.v);
[p2.ac,p2.al]=rotate_hydro([p2_pos(2) cpies45_pos(2)],[a_pos(1) cpies45_pos(1)],p2.u,p2.v);

B_dx=1000*sw_dist([a_pos(2) b_pos(2)],[a_pos(1) b_pos(1)],'km');
C_dx=1000*sw_dist([a_pos(2) c_pos(2)],[a_pos(1) c_pos(1)],'km');
D_dx=1000*sw_dist([a_pos(2) d_pos(2)],[a_pos(1) d_pos(1)],'km');
E_dx=1000*sw_dist([a_pos(2) e_pos(2)],[a_pos(1) e_pos(1)],'km');
F_dx=1000*sw_dist([a_pos(2) f_pos(2)],[a_pos(1) f_pos(1)],'km');
G_dx=1000*sw_dist([a_pos(2) g_pos(2)],[a_pos(1) g_pos(1)],'km');
cpies34_dx=1000*sw_dist([a_pos(2) cpies34_pos(2)],[a_pos(1) cpies34_pos(1)],'km');
cpies45_dx=1000*sw_dist([a_pos(2) cpies45_pos(2)],[a_pos(1) cpies45_pos(1)],'km');
p2_dx=1000*sw_dist([a_pos(2) p2_pos(2)],[a_pos(1) p2_pos(1)],'km');


% assign decorrelation length scales, noise levels, and difference in scale
% between passes 1 and 2

xc=63*1000; %horizontal decorrelation length
zc=1913; %vertical decorrelation length

x_corr_func=@(x) exp(-(x(:)/xc).^2).*cos(pi.*x(:)./(2.*xc));
z_corr_func=@(z) exp(-(z(:)/zc).^2);

int_ac2=nan(size(xgrid,1),size(xgrid,2),size(A.ac,2));
int_ac=nan(size(xgrid,1),size(xgrid,2),size(A.ac,2));
parfor time=1:size(A.u,2)
    Noise2=[];
    ac_obs_time=[];
    dx=[];
    weights_ac=[];
    cross_corr=[];
    ratio_ac=[];
    weight_corr=[];
    ac_obs=[];
    temp1=[];
    temp2=[];
    grid_dist=[];
    Mz=[];
    Iz=[];
    M=[];
    I=[];
    dz_dist=[];
    ac_obs_anom=[];
    temp3=[];
    xc=63*1000; %horizontal decorrelation length
    zc=1913; %vertical decorrelation length
    ac_obs_time=[A.ac(~isnan(A.ac(:,time)),:);B.ac(~isnan(B.ac(:,time)),:);C.ac(~isnan(C.ac(:,time)),:);D.ac(~isnan(D.ac(:,time)),:);E.ac(~isnan(E.ac(:,time)),:);F.ac(~isnan(F.ac(:,time)),:);p2.ac.';G.ac(~isnan(G.ac(:,time)),:);cpies34_(~isnan(cpies34_(:,time)),:);cpies45_(~isnan(cpies45_(:,time)),:)];
    dz_obs=[A.z(~isnan(A.v(:,time)),time);B.z(~isnan(B.v(:,time)),time);C.z(~isnan(C.v(:,time)),time);D.z(~isnan(D.v(:,time)),time);E.z(~isnan(E.v(:,time)),time);F.z(~isnan(F.v(:,time)),time);4183;G.z(~isnan(G.v(:,time)),time);dpth34_(~isnan(cpies34_(:,time)));dpth45_(~isnan(cpies45_(:,time)))];
    ac_obs=[A.ac(~isnan(A.ac(:,time)),time);B.ac(~isnan(B.ac(:,time)),time);C.ac(~isnan(C.ac(:,time)),time);D.ac(~isnan(D.ac(:,time)),time);E.ac(~isnan(E.ac(:,time)),time);F.ac(~isnan(F.ac(:,time)),time);p2.ac(time);G.ac(~isnan(G.ac(:,time)),time);cpies34_(~isnan(cpies34_(:,time)),time);cpies45_(~isnan(cpies45_(:,time)),time)];
    % set up noise matrix
    Noise2=zeros(sum(~isnan(A.u(:,time)))+sum(~isnan(B.u(:,time)))+sum(~isnan(C.u(:,time)))+sum(~isnan(D.u(:,time)))+sum(~isnan(E.u(:,time)))+sum(~isnan(F.u(:,time)))+sum(~isnan(G.u(:,time)))+sum(~isnan(cpies34_(:,time)))+sum(~isnan(cpies45_(:,time)))+1);
    for i=1:length(Noise2)
        Noise2(i)=.07;
    end
    % set up ratio matrix
    ratio_ac=zeros(length(Noise2),length(Noise2));
    for i=1:length(Noise2)
        for j=1:length(Noise2)
            if i==j
                ratio_ac(i,j)=Noise2(i)/nanvar(ac_obs_time(i,:)); % seems bad that these values are so big. is this right? what was it on the one obs case?
            end
        end
    end
    % set up dx matrix
    for i=1:sum(~isnan(A.v(:,time)))
        dx(i)=0; % I am only going to interpolate starting at A - will treat inshore portion differently
    end
    for i=sum(~isnan(A.v(:,time)))+1:sum(~isnan(A.v(:,time)))+sum(~isnan(B.v(:,time)))
        dx(i)=B_dx;
    end
    for i=sum(~isnan(A.v(:,time)))+sum(~isnan(B.v(:,time)))+1:sum(~isnan(A.v(:,time)))+sum(~isnan(B.v(:,time)))+sum(~isnan(C.v(:,time)))
        dx(i)=C_dx;
    end
    for i=sum(~isnan(A.v(:,time)))+sum(~isnan(B.v(:,time)))+sum(~isnan(C.v(:,time)))+1:sum(~isnan(A.v(:,time)))+sum(~isnan(B.v(:,time)))+sum(~isnan(C.v(:,time)))+sum(~isnan(D.v(:,time)))
        dx(i)=D_dx;
    end
    for i=sum(~isnan(A.v(:,time)))+sum(~isnan(B.v(:,time)))+sum(~isnan(C.v(:,time)))+sum(~isnan(D.v(:,time)))+1:sum(~isnan(A.v(:,time)))+sum(~isnan(B.v(:,time)))+sum(~isnan(C.v(:,time)))+sum(~isnan(D.v(:,time)))+sum(~isnan(E.v(:,time)))
        dx(i)=E_dx;
    end
    for i=sum(~isnan(A.v(:,time)))+sum(~isnan(B.v(:,time)))+sum(~isnan(C.v(:,time)))+sum(~isnan(D.v(:,time)))+sum(~isnan(E.v(:,time)))+1:sum(~isnan(A.v(:,time)))+sum(~isnan(B.v(:,time)))+sum(~isnan(C.v(:,time)))+sum(~isnan(D.v(:,time)))+sum(~isnan(E.v(:,time)))+sum(~isnan(F.v(:,time)))
        dx(i)=F_dx;
    end
    for i=sum(~isnan(A.v(:,time)))+sum(~isnan(B.v(:,time)))+sum(~isnan(C.v(:,time)))+sum(~isnan(D.v(:,time)))+sum(~isnan(E.v(:,time)))+sum(~isnan(F.v(:,time)))+1
        dx(i)=p2_dx;
    end
    for i=sum(~isnan(A.v(:,time)))+sum(~isnan(B.v(:,time)))+sum(~isnan(C.v(:,time)))+sum(~isnan(D.v(:,time)))+sum(~isnan(E.v(:,time)))+sum(~isnan(F.v(:,time)))+2:sum(~isnan(A.v(:,time)))+sum(~isnan(B.v(:,time)))+sum(~isnan(C.v(:,time)))+sum(~isnan(D.v(:,time)))+sum(~isnan(E.v(:,time)))+sum(~isnan(F.v(:,time)))+sum(~isnan(G.v(:,time)))+1
        dx(i)=G_dx;
    end
    for i=sum(~isnan(A.v(:,time)))+sum(~isnan(B.v(:,time)))+sum(~isnan(C.v(:,time)))+sum(~isnan(D.v(:,time)))+sum(~isnan(E.v(:,time)))+sum(~isnan(F.v(:,time)))+sum(~isnan(G.v(:,time)))+2:sum(~isnan(A.v(:,time)))+sum(~isnan(B.v(:,time)))+sum(~isnan(C.v(:,time)))+sum(~isnan(D.v(:,time)))+sum(~isnan(E.v(:,time)))+sum(~isnan(F.v(:,time)))+sum(~isnan(G.v(:,time)))+sum(~isnan(cpies34_(:,time)))+1
        dx(i)=cpies34_dx;
    end
    for i=sum(~isnan(A.v(:,time)))+sum(~isnan(B.v(:,time)))+sum(~isnan(C.v(:,time)))+sum(~isnan(D.v(:,time)))+sum(~isnan(E.v(:,time)))+sum(~isnan(F.v(:,time)))+sum(~isnan(G.v(:,time)))+sum(~isnan(cpies34_(:,time)))+2:sum(~isnan(A.v(:,time)))+sum(~isnan(B.v(:,time)))+sum(~isnan(C.v(:,time)))+sum(~isnan(D.v(:,time)))+sum(~isnan(E.v(:,time)))+sum(~isnan(F.v(:,time)))+sum(~isnan(G.v(:,time)))+sum(~isnan(cpies34_(:,time)))+sum(~isnan(cpies45_(:,time)))+1
        dx(i)=cpies45_dx;
    end
    
    % now get cross correlations between instruments, grid points
    weight_corr=nan(length(Noise2),length(zgrid),size(zgrid,2));
    for i=1:length(Noise2)
        for j=1:size(zgrid,1) % was length(zgrid)
            for k=1:size(zgrid,2)
                weight_corr(i,j,k)=x_corr_func(abs(xgrid(j,k)-dx(i)))*z_corr_func(abs(zgrid(j,k)-dz_obs(i)));
            end
        end
    end
    
    % get cross corr between instruments
    cross_corr=nan(length(Noise2),length(Noise2));
    for i=1:length(Noise2)
        for j=1:length(Noise2)
            cross_corr(i,j)=x_corr_func(abs(dx(i)-dx(j)))*z_corr_func(abs(dz_obs(i)-dz_obs(j)));
        end
    end
    
    % solve for weights
    weights_ac=nan(size(weight_corr,1),size(zgrid,1),size(zgrid,2));
    for j=1:size(zgrid,1)
        for k=1:size(zgrid,2)
            weights_ac(:,j,k)=(ratio_ac+cross_corr)\weight_corr(:,j,k);
        end
    end

    for j=1:size(zgrid,1)
        for k=1:size(zgrid,2)
            temp1(j,k)=weights_ac(:,j,k).'*ac_obs;
        end
    end
    int_ac(:,:,time)=temp1(:,:);
    % this will give you the error predicted from the OI but will slow down
    % the total computation time significantly
%     background_error=ones(size(xgrid,1),size(xgrid,2));
%     for i=1:size(xgrid,1)
%         for j=1:size(xgrid,2)
%             temp3(i,j)=background_error(i,j)-(weight_corr(:,i,j).'*inv(ratio+cross_corr)*weight_corr(:,i,j));
%         end
%     end
%     analysis_error(:,:,time)=temp3(:,:);
    
    % repeat on small scale
    small_scale_ratio=.2; % may want to test different values
    
    xc=63*1000*small_scale_ratio; %horizontal decorrelation length
    zc=1913*small_scale_ratio; %vertical decorrelation length
    for i=1:size(dz_obs,1)
        for j=1:size(zgrid,1)
            dz_dist(i,j)=abs(zgrid(j,1)-dz_obs(i));
        end
    end
    
    for i=1:size(dz_dist,1)
        [Mz(i),Iz(i)] = min(dz_dist(i,:));
    end
    
    for i=1:size(dx,2)
        for j=1:size(xgrid,2)
            grid_dist(i,j)=abs(xgrid(1,j)-dx(i));
        end
    end
    
    for i=1:size(dx,2)
        [M(i),I(i)] = min(grid_dist(i,:));
    end
    
    for i=1:size(ac_obs) % this will come in for the small scale part
        ac_obs_anom(i)=ac_obs(i)-temp1(Iz(i),I(i));
    end
    
    cross_corr=nan(length(Noise2),length(Noise2));
    for i=1:length(Noise2)
        for j=1:length(Noise2)
            cross_corr(i,j)=x_corr_func(abs(dx(i)-dx(j)))*z_corr_func(abs(dz_obs(i)-dz_obs(j)));
        end
    end
    
    % weights are different than before
    
    weights_ac=nan(size(weight_corr,1),size(zgrid,1),size(zgrid,2));
    for j=1:size(zgrid,1)
        for k=1:size(zgrid,2)
            weights_ac(:,j,k)=(ratio_ac+cross_corr)\weight_corr(:,j,k);
        end
    end
    
    for j=1:size(zgrid,1)
        for k=1:size(zgrid,2)
            temp2(j,k)=weights_ac(:,j,k).'*ac_obs_anom.';
        end
    end
    int_ac2(:,:,time)=temp2(:,:);
end
vel_d1=int_ac+int_ac2;



% deployment 2

load('act2-A_cm_velocities.mat')
load('act2-B_cm_velocities.mat')
load('act2-C_cm_velocities.mat')
load('act2-D_cm_velocities.mat')
load('act2-E_cm_velocities.mat')
load('act2-F_cm_velocities.mat')
load('act2-G_cm_velocities.mat')
load('ACT_CPIES_velocities_20h.mat')
load('P2_241_currents.mat')

hr=20; % interpolate once every 20 hours

A.u=A.u_adcp(:,182:hr:11162);
B.u=B.u(:,175:hr:11155);
C.u=C.u(:,152:hr:11132);
D.u=D.u(:,131:hr:11111);
E.u=E.u(:,104:hr:11084);
F.u=F.u(:,58:hr:11038);
G.u=G.u(:,10:hr:10990);

A.v=A.v_adcp(:,182:hr:11162);
B.v=B.v(:,175:hr:11155);
C.v=C.v(:,152:hr:11132);
D.v=D.v(:,131:hr:11111);
E.v=E.v(:,104:hr:11084);
F.v=F.v(:,58:hr:11038);
G.v=G.v(:,10:hr:10990);

A.z=A.z_adcp(:,182:hr:11162);
B.z=B.z(:,175:hr:11155);
C.z=C.z(:,152:hr:11132);
D.z=D.z(:,131:hr:11111);
E.z=E.z(:,104:hr:11084);
F.z=F.z(:,58:hr:11038);
G.z=G.z(:,10:hr:10990);

time2=G.sd(1,10:hr:10990);

cpies34=v34(:,700:1249);
cpies45=v45(:,700:1249);

% instrument locations, use to get distances between instruments

a_pos = [27.595,-33.5583];
b_pos = [27.6428,-33.6674];
c_pos = [27.7152,-33.7996];
d_pos = [27.8603,-34.0435];
e_pos = [28.0316,-34.29];
f_pos = [28.17,-34.5380];
g_pos = [28.3453,-34.8213];
cpies34_pos = [(28.3453+28.6248)/2,-(34.8213+35.2453)/2];
cpies45_pos = [(28.6248+28.9)/2,-(35.2453+35.7338)/2];
p2_pos = [p2.lon(1),p2.lat(1)];

x=0:500:1000*sw_dist([a_pos(2) cpies45_pos(2)],[a_pos(1) cpies45_pos(1)],'km');
z=0:20:4500; 
for i=1:length(z)
    xgrid(i,:)=x(:);
end
for i=1:length(x)
    zgrid(:,i)=z(:);
end

p2.u=p2.u(700*20:20:1249*20);
p2.v=p2.v(700*20:20:1249*20);

[A.ac,A.al]=rotate_hydro([a_pos(2) cpies45_pos(2)],[a_pos(1) cpies45_pos(1)],A.u,A.v);
[B.ac,B.al]=rotate_hydro([a_pos(2) cpies45_pos(2)],[a_pos(1) cpies45_pos(1)],B.u,B.v);
[C.ac,C.al]=rotate_hydro([a_pos(2) cpies45_pos(2)],[a_pos(1) cpies45_pos(1)],C.u,C.v);
[D.ac,D.al]=rotate_hydro([a_pos(2) cpies45_pos(2)],[a_pos(1) cpies45_pos(1)],D.u,D.v);
[E.ac,E.al]=rotate_hydro([a_pos(2) cpies45_pos(2)],[a_pos(1) cpies45_pos(1)],E.u,E.v);
[F.ac,F.al]=rotate_hydro([a_pos(2) cpies45_pos(2)],[a_pos(1) cpies45_pos(1)],F.u,F.v);
[G.ac,G.al]=rotate_hydro([a_pos(2) cpies45_pos(2)],[a_pos(1) cpies45_pos(1)],G.u,G.v);
[p2.ac,p2.al]=rotate_hydro([p2_pos(2) cpies45_pos(2)],[a_pos(1) cpies45_pos(1)],p2.u,p2.v);

B_dx=1000*sw_dist([a_pos(2) b_pos(2)],[a_pos(1) b_pos(1)],'km');
C_dx=1000*sw_dist([a_pos(2) c_pos(2)],[a_pos(1) c_pos(1)],'km');
D_dx=1000*sw_dist([a_pos(2) d_pos(2)],[a_pos(1) d_pos(1)],'km');
E_dx=1000*sw_dist([a_pos(2) e_pos(2)],[a_pos(1) e_pos(1)],'km');
F_dx=1000*sw_dist([a_pos(2) f_pos(2)],[a_pos(1) f_pos(1)],'km');
G_dx=1000*sw_dist([a_pos(2) g_pos(2)],[a_pos(1) g_pos(1)],'km');
cpies34_dx=1000*sw_dist([a_pos(2) cpies34_pos(2)],[a_pos(1) cpies34_pos(1)],'km');
cpies45_dx=1000*sw_dist([a_pos(2) cpies45_pos(2)],[a_pos(1) cpies45_pos(1)],'km');
p2_dx=1000*sw_dist([a_pos(2) p2_pos(2)],[a_pos(1) p2_pos(1)],'km');

%downscale cpies to one measurement every 100 m to make computation faster
cpies45_=cpies45(1:5:end,:);
cpies34_=cpies34(1:5:end,:);
dpth45_=dpth45(1:5:end);
dpth34_=dpth34(1:5:end);

% assign decorrelation length scales, noise levels, and difference in scale
% between passes 1 and 2
xc=48*1000; %horizontal decorrelation length
zc=1009;
x_corr_func=@(x) exp(-(x(:)/xc).^2).*cos(pi.*x(:)./(2.*xc));
z_corr_func=@(z) exp(-(z(:)/zc).^2);

int_act=nan(size(xgrid,1),size(xgrid,2),size(A.ac,2));
int_act2=nan(size(xgrid,1),size(xgrid,2),size(A.ac,2));
parfor time=1:size(A.u,2)
    xc=48*1000; %horizontal decorrelation length
    zc=1009; %vertical decorrelation length
    Noise2=[];
    ac_obs_time=[];
    ac_obs_anom=[];
    dx=[];
    weights_ac=[];
    cross_corr=[];
    ratio_ac=[];
    weight_corr=[];
    ac_obs=[];
    M=[];
    I=[];
    Mz=[];
    Iz=[];
    temp1=[];
    temp2=[];
    grid_dist=[];
    dz_dist=[];
    temp3=[];
    ac_obs_time=[A.ac(~isnan(A.ac(:,time)),:);B.ac(~isnan(B.ac(:,time)),:);C.ac(~isnan(C.ac(:,time)),:);D.ac(~isnan(D.ac(:,time)),:);E.ac(~isnan(E.ac(:,time)),:);F.ac(~isnan(F.ac(:,time)),:);p2.ac.';G.ac(~isnan(G.ac(:,time)),:);cpies34_(~isnan(cpies34_(:,time)),:);cpies45_(~isnan(cpies45_(:,time)),:)];
    dz_obs=[A.z(~isnan(A.v(:,time)),time);B.z(~isnan(B.v(:,time)),time);C.z(~isnan(C.v(:,time)),time);D.z(~isnan(D.v(:,time)),time);E.z(~isnan(E.v(:,time)),time);F.z(~isnan(F.v(:,time)),time);4183;G.z(~isnan(G.v(:,time)),time);dpth34_(~isnan(cpies34_(:,time)));dpth45_(~isnan(cpies45_(:,time)))];
    ac_obs=[A.ac(~isnan(A.ac(:,time)),time);B.ac(~isnan(B.ac(:,time)),time);C.ac(~isnan(C.ac(:,time)),time);D.ac(~isnan(D.ac(:,time)),time);E.ac(~isnan(E.ac(:,time)),time);F.ac(~isnan(F.ac(:,time)),time);p2.ac(time);G.ac(~isnan(G.ac(:,time)),time);cpies34_(~isnan(cpies34_(:,time)),time);cpies45_(~isnan(cpies45_(:,time)),time)];
    % set up noise matrix
    Noise2=zeros(sum(~isnan(A.u(:,time)))+sum(~isnan(B.u(:,time)))+sum(~isnan(C.u(:,time)))+sum(~isnan(D.u(:,time)))+sum(~isnan(E.u(:,time)))+sum(~isnan(F.u(:,time)))+sum(~isnan(G.u(:,time)))+sum(~isnan(cpies34_(:,time)))+sum(~isnan(cpies45_(:,time)))+1);
    for i=1:length(Noise2)
        Noise2(i)=.07;
    end
    % set up ratio matrix
    ratio_ac=zeros(length(Noise2),length(Noise2));
    for i=1:length(Noise2)
        for j=1:length(Noise2)
            if i==j
                ratio_ac(i,j)=Noise2(i)/nanvar(ac_obs_time(i,:));
            end
        end
    end
    % set up dx matrix
    for i=1:sum(~isnan(A.v(:,time)))
        dx(i)=0; % I am only going to interpolate starting at A - will treat inshore portion differently
    end
    for i=sum(~isnan(A.v(:,time)))+1:sum(~isnan(A.v(:,time)))+sum(~isnan(B.v(:,time)))
        dx(i)=B_dx;
    end
    for i=sum(~isnan(A.v(:,time)))+sum(~isnan(B.v(:,time)))+1:sum(~isnan(A.v(:,time)))+sum(~isnan(B.v(:,time)))+sum(~isnan(C.v(:,time)))
        dx(i)=C_dx;
    end
    for i=sum(~isnan(A.v(:,time)))+sum(~isnan(B.v(:,time)))+sum(~isnan(C.v(:,time)))+1:sum(~isnan(A.v(:,time)))+sum(~isnan(B.v(:,time)))+sum(~isnan(C.v(:,time)))+sum(~isnan(D.v(:,time)))
        dx(i)=D_dx;
    end
    for i=sum(~isnan(A.v(:,time)))+sum(~isnan(B.v(:,time)))+sum(~isnan(C.v(:,time)))+sum(~isnan(D.v(:,time)))+1:sum(~isnan(A.v(:,time)))+sum(~isnan(B.v(:,time)))+sum(~isnan(C.v(:,time)))+sum(~isnan(D.v(:,time)))+sum(~isnan(E.v(:,time)))
        dx(i)=E_dx;
    end
    for i=sum(~isnan(A.v(:,time)))+sum(~isnan(B.v(:,time)))+sum(~isnan(C.v(:,time)))+sum(~isnan(D.v(:,time)))+sum(~isnan(E.v(:,time)))+1:sum(~isnan(A.v(:,time)))+sum(~isnan(B.v(:,time)))+sum(~isnan(C.v(:,time)))+sum(~isnan(D.v(:,time)))+sum(~isnan(E.v(:,time)))+sum(~isnan(F.v(:,time)))
        dx(i)=F_dx;
    end
    for i=sum(~isnan(A.v(:,time)))+sum(~isnan(B.v(:,time)))+sum(~isnan(C.v(:,time)))+sum(~isnan(D.v(:,time)))+sum(~isnan(E.v(:,time)))+sum(~isnan(F.v(:,time)))+1
        dx(i)=p2_dx;
    end
    for i=sum(~isnan(A.v(:,time)))+sum(~isnan(B.v(:,time)))+sum(~isnan(C.v(:,time)))+sum(~isnan(D.v(:,time)))+sum(~isnan(E.v(:,time)))+sum(~isnan(F.v(:,time)))+2:sum(~isnan(A.v(:,time)))+sum(~isnan(B.v(:,time)))+sum(~isnan(C.v(:,time)))+sum(~isnan(D.v(:,time)))+sum(~isnan(E.v(:,time)))+sum(~isnan(F.v(:,time)))+sum(~isnan(G.v(:,time)))+1
        dx(i)=G_dx;
    end
    for i=sum(~isnan(A.v(:,time)))+sum(~isnan(B.v(:,time)))+sum(~isnan(C.v(:,time)))+sum(~isnan(D.v(:,time)))+sum(~isnan(E.v(:,time)))+sum(~isnan(F.v(:,time)))+sum(~isnan(G.v(:,time)))+2:sum(~isnan(A.v(:,time)))+sum(~isnan(B.v(:,time)))+sum(~isnan(C.v(:,time)))+sum(~isnan(D.v(:,time)))+sum(~isnan(E.v(:,time)))+sum(~isnan(F.v(:,time)))+sum(~isnan(G.v(:,time)))+sum(~isnan(cpies34_(:,time)))+1
        dx(i)=cpies34_dx;
    end
    for i=sum(~isnan(A.v(:,time)))+sum(~isnan(B.v(:,time)))+sum(~isnan(C.v(:,time)))+sum(~isnan(D.v(:,time)))+sum(~isnan(E.v(:,time)))+sum(~isnan(F.v(:,time)))+sum(~isnan(G.v(:,time)))+sum(~isnan(cpies34_(:,time)))+2:sum(~isnan(A.v(:,time)))+sum(~isnan(B.v(:,time)))+sum(~isnan(C.v(:,time)))+sum(~isnan(D.v(:,time)))+sum(~isnan(E.v(:,time)))+sum(~isnan(F.v(:,time)))+sum(~isnan(G.v(:,time)))+sum(~isnan(cpies34_(:,time)))+sum(~isnan(cpies45_(:,time)))+1
        dx(i)=cpies45_dx;
    end
    
    % now get cross correlations between instruments, grid points
    weight_corr=nan(length(Noise2),size(zgrid,1),size(zgrid,2));
    for i=1:length(Noise2)
        for j=1:size(zgrid,1)
            for k=1:size(zgrid,2)
                weight_corr(i,j,k)=x_corr_func(abs(xgrid(j,k)-dx(i)))*z_corr_func(abs(zgrid(j,k)-dz_obs(i)));
            end
        end
    end
    % get cross corr between instruments
    cross_corr=nan(length(Noise2),length(Noise2));
    for i=1:length(Noise2)
        for j=1:length(Noise2)
            cross_corr(i,j)=x_corr_func(abs(dx(i)-dx(j)))*z_corr_func(abs(dz_obs(i)-dz_obs(j)));
        end
    end
    
    % solve for weights
    weights_ac=nan(size(weight_corr,1),size(zgrid,1),size(zgrid,2));
    for j=1:size(zgrid,1)
        for k=1:size(zgrid,2)
            weights_ac(:,j,k)=(ratio_ac+cross_corr)\weight_corr(:,j,k);
        end
    end

    for j=1:size(zgrid,1)
        for k=1:size(zgrid,2)
            temp1(j,k)=weights_ac(:,j,k).'*ac_obs;
        end
    end
    int_act(:,:,time)=temp1(:,:);
    % this will give you the error predicted from the OI but will slow down
    % the total computation time significantly
    %     background_error=ones(size(xgrid,1),size(xgrid,2));
    %     for i=1:size(xgrid,1)
    %         for j=1:size(xgrid,2)
    %             temp3(i,j)=background_error(i,j)-(weight_corr(:,i,j).'*inv(ratio+cross_corr)*weight_corr(:,i,j));
    %         end
    %     end
    %     analysis_error2(:,:,time)=temp3(:,:);
    
    
    %small scale
    small_scale_ratio=.2;
    
    xc=48*1000*small_scale_ratio; %horizontal decorrelation length
    zc=1009*small_scale_ratio; %vertical decorrelation length
    for i=1:size(dz_obs,1)
        for j=1:size(zgrid,1)
            dz_dist(i,j)=abs(zgrid(j,1)-dz_obs(i));
        end
    end
    
    for i=1:size(dz_dist,1)
        [Mz(i),Iz(i)] = min(dz_dist(i,:));
    end
    
    for i=1:size(dx,2)
        for j=1:size(xgrid,2)
            grid_dist(i,j)=abs(xgrid(1,j)-dx(i));
        end
    end
    
    for i=1:size(dx,2)
        [M(i),I(i)] = min(grid_dist(i,:));
    end
    
    for i=1:size(ac_obs) % this will come in for the small scale part
        ac_obs_anom(i)=ac_obs(i)-temp1(Iz(i),I(i));
    end
    
    cross_corr=nan(length(Noise2),length(Noise2));
    for i=1:length(Noise2)
        for j=1:length(Noise2)
            cross_corr(i,j)=x_corr_func(abs(dx(i)-dx(j)))*z_corr_func(abs(dz_obs(i)-dz_obs(j)));
        end
    end
    
    % weights are different than before
    weights_ac=nan(size(weight_corr,1),size(zgrid,1),size(zgrid,2));
    for j=1:size(zgrid,1)
        for k=1:size(zgrid,2)
            weights_ac(:,j,k)=(ratio_ac+cross_corr)\weight_corr(:,j,k);
        end
    end
    
    for j=1:size(zgrid,1)
        for k=1:size(zgrid,2)
            temp2(j,k)=weights_ac(:,j,k).'*ac_obs_anom.';
        end
    end
    int_act2(:,:,time)=temp2(:,:);
end

vel_d2=int_act_int_act2;
