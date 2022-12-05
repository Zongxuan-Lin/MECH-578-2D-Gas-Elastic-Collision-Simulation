clc;clear all;close all;
%% given parameters and table setup
t=0; %time start at 0
ancol=0;
boxsize=21; %box size
num=(boxsize-1)^2;%A total of 400 molecules in the simulation
bconst=1.0; %fictional boltzman constant in the problem
temp=1/2;%initial temperature
ncol=1;%number of collisions
moleculed=0.1; %diameter and mass of molecule
molecuelm=1;
m2=molecuelm*2;
rangepair=[-5 -3 0;5 7 10];
bolfac=4;%size of sample scaling factor for boltz
intervals1=trapezoid([-2 -1 1 2], num);
intervals2=trapezoid([-2 -1 1 2], num);
intervals3=triangle([-2 0 1], num*bolfac);%for boltzmann distribution function
intervals4=triangle([-1 0 2], num*bolfac);%for boltzmann distribution function
%vel=radinrange(intervals1(1:num,:),intervals2(1:num,:),num);
vel=randomdis([2 4],num);
velbol=radinrange(intervals3,intervals4,num*bolfac);%boltzman velocity 
%vel=settemp(vel,temp,num);
velbol=settemp(velbol,temp,num*bolfac);%boltzman velocity 
pos=zeros(num,2);%position table of the particles
volrec=zeros(ncol+1,num,2);
posrec=zeros(ncol+1,num,2);
spdrec=zeros(ncol+1,num,1);
speed=zeros(num,1);
ancols=zeros(ncol+1,1);
time=zeros(ncol+1,1);
dts=zeros(ncol,1);
Cnum=0;
afb=(num+1)*ones(num,1);%affected ball id in collision
%% initial Position
for i=1:num
    x=fix(i/(boxsize-1))+1;
    y=mod(i,20);
    if y==0
        y=20; x=x-1;
    end
    pos(i,1)=x;pos(i,2)=y;%initialize position
end
speed=sqrt(sum(vel.^2,2));
volrec(1,:,:)=vel; posrec(1,:,:)=pos; spdrec(1,:,:)=speed;
%% First collision time calculation
comb=nchoosek(num,2);%total number of pairs to consider
colliBall=zeros(comb,3);%potential collision between particles table, update between every collision
colliWall=zeros(num,3);%potential collision between particles and wall table, update between every collision
temppos=1; %navigator for collision table
%first run of collision between ball simulation
for i=1:(num-1)
    for j=(i+1):num
        colliBall(temppos,1:2)=[i j];%initializing pair id
        rel_pos=pos(i,:)-pos(j,:);
        rel_vel=vel(i,:)-vel(j,:);
        mark1=dot(rel_pos, rel_vel);%Dot Product of rel position and rel velocity
        velsq=norm(rel_vel)^2;
        possq=norm(rel_pos)^2;
        mark2=mark1^2-velsq*(possq-moleculed^2);
        if mark1<0 && mark2>=0
            th1=(-mark1+sqrt(mark2))/velsq; th2=(-mark1-sqrt(mark2))/velsq;
            colliBall(temppos,3)=compare(th1,th2);
        else
            colliBall(temppos,3)=inf;
        end
        temppos=temppos+1;
    end
end
%wall 1 left  0
%wall 2 right     21
%wall 3 bottom    0
%wall 4 top   21
for i=1:num
    colliWall(i,1)=i;
    posx=pos(i,1);  posy=pos(i,2);  velx=vel(i,1);  vely=vel(i,2);
    tx=0;ty=0;wx=0;wy=0;
    if (velx>0)
        tx=(boxsize-posx-moleculed/2)/velx;wx=2;
    elseif(velx<0)
        tx=(0-posx+moleculed/2)/velx;wx=1;
    else
        tx=inf;wx=5;
    end
    if (vely>0)
        ty=(boxsize-posy-moleculed/2)/vely;wy=4;
    elseif(vely<0)
        ty=(0-posy+moleculed/2)/vely;wy=3;
    else
        ty=inf;wy=5;
    end
    if(tx<0)
        tx=0;
    end
    if(ty<0)
        ty=0;
    end
    if (tx<ty)
        colliWall(i,2)=wx;colliWall(i,3)=tx;
    elseif(tx>ty)
        colliWall(i,2)=wy; colliWall(i,3)=ty;
    elseif(tx~=inf && tx==ty)
        colliWall(i,2)=6;colliWall(i,3)=ty;%code 6 collision with both wall together
    else
        colliWall(i,2)=5;colliWall(i,3)=inf;
    end
end

%% first post-collision Velocity Analysis 
%find the first collision happend
%2d elastic collision between balls calc credit to fiflulu
   temporarystorage1=colliBall(:,3);temporarystorage2=colliWall(:,3);
   [mint_ball, minBallInd] = min(temporarystorage1);
   [mint_wall, minwallInd] = min(temporarystorage2);
   ancol=ancol+1;
%% Collision with Wall Velocity Analysis 
if mint_wall<mint_ball%hit wall first
    wallid=colliWall(minwallInd,2);deltat=mint_wall;
    t=t+deltat;
    dts(1)=deltat;
    time(2)=t;
    pos=pos+vel*deltat;
    if wallid <=2%x direction
        vel(minwallInd,1)=-vel(minwallInd,1);
    elseif wallid>2 && wallid <=4%y direction
        vel(minwallInd,2)=-vel(minwallInd,2);
    elseif wallid==6%in the corner
        vel(minwallInd,1)=-vel(minwallInd,1);
        vel(minwallInd,2)=-vel(minwallInd,2);
    end
    afb(1)=minwallInd;Cnum=1;
    for i=1:comb
        xx=colliBall(i,1);yx=colliBall(i,2);
        dist=norm(pos(xx,:)-pos(yx,:));
        if dist<moleculed
            Cnum=Cnum+2;
            afb(Cnum-1)=xx;afb(Cnum)=yx;
            r12unit=(pos(yx,:)-pos(xx,:))/(norm(pos(yx,:)-pos(xx,:)));
            cmnterm=dot((vel(xx,:)-vel(yx,:)),r12unit)*r12unit;
            vb1=vel(xx,:)-cmnterm;vel(xx,:)=vb1;
            vb2=vel(yx,:)+cmnterm;vel(yx,:)=vb2;
            speed(xx)=norm(vel(xx,:));speed(yx)=norm(vel(yx,:));
            posxx=pos(xx,:)+vb1*(moleculed-dist)/speed(xx);pos(xx,:)=posxx;
            posyx=pos(yx,:)+vb2*(moleculed-dist)/speed(yx);pos(yx,:)=posyx;
            ancol=ancol+1;
        end
    end
    volrec(2,:,:)=vel;posrec(2,:,:)=pos;spdrec(2,:,:)=speed;
    colliWall(:,3)=colliWall(:,3)-deltat;colliWall(:,3)=colliWall(:,3)-deltat;

%% Collision with Ball Velocity Analysis 
elseif mint_ball<mint_wall  %hit ball first
    ball1=colliBall(minBallInd,1); ball2=colliBall(minBallInd,2);
    deltat=mint_ball; t=t+deltat;dts(1)=deltat;
    time(2)=t;
    % Simple procedure found in Halie
    r12unit=(pos(ball2,:)-pos(ball1,:))/(norm(pos(ball2,:)-pos(ball1,:)));
    cmnterm=dot((vel(ball1,:)-vel(ball2,:)),r12unit)*r12unit;
    vb1=vel(ball1,:)-cmnterm;vel(ball1,:)=vb1;
    vb2=vel(ball2,:)+cmnterm;vel(ball2,:)=vb2;
    speed(ball1)=norm(vel(ball1,:));speed(ball2)=norm(vel(ball2,:));
    pos=pos+vel*deltat;
    afb(1)=ball1;afb(2)=ball2;Cnum=2;
    for i=1:comb
        xx=colliBall(i,1);yx=colliBall(i,2);
        dist=norm(pos(xx,:)-pos(yx,:));
        if dist<moleculed
            Cnum=Cnum+2;
            afb(Cnum-1)=xx;afb(Cnum)=yx;
            r12unit=(pos(yx,:)-pos(xx,:))/(norm(pos(yx,:)-pos(xx,:)));
            cmnterm=dot((vel(xx,:)-vel(yx,:)),r12unit)*r12unit;
            vb1=vel(xx,:)-cmnterm;vel(xx,:)=vb1;
            vb2=vel(yx,:)+cmnterm;vel(yx,:)=vb2;
            speed(xx)=norm(vel(xx,:));speed(yx)=norm(vel(yx,:));
            posxx=pos(xx,:)+vb1*(moleculed-dist)/speed(xx);pos(xx,:)=posxx;
            posyx=pos(yx,:)+vb2*(moleculed-dist)/speed(yx);pos(yx,:)=posyx;
            ancol=ancol+1;
        end
    end
    volrec(2,:,:)=vel;posrec(2,:,:)=pos;spdrec(2,:,:)=speed;
    colliWall(:,3)=colliWall(:,3)-deltat;colliWall(:,3)=colliWall(:,3)-deltat;

%% Collision with Wall and Ball Simultaneous Velocity Analysis 
else%hit wall and ball simultaneously
    wallid=colliWall(minwallInd,2);
    deltat=colliWall(minwallInd,3);
    t=t+deltat;time(2)=t;dts(1)=deltat;
    if wallid <=2
        vel(minwallInd,1)=-vel(minwallInd,1);
    elseif wallid>2 && wallid <=4
        vel(minwallInd,2)=-vel(minwallInd,2);
    elseif wallid==6%in the corner
        vel(minwallInd,1)=-vel(minwallInd,1);
        vel(minwallInd,2)=-vel(minwallInd,2);
    end

    ball1=colliBall(minBallInd,1); ball2=colliBall(minBallInd,2);
    % Simple procedure found in Halie
    r12unit=(pos(ball2,:)-pos(ball1,:))/(norm(pos(ball2,:)-pos(ball1,:)));
    cmnterm=dot((vel(ball1,:)-vel(ball2,:)),r12unit)*r12unit;
    vb1=vel(ball1,:)-cmnterm;vel(ball1,:)=vb1;
    vb2=vel(ball2,:)+cmnterm;vel(ball2,:)=vb2;
    speed(ball1)=norm(vel(ball1,:));speed(ball2)=norm(vel(ball2,:));
    pos=pos+vel*deltat;
    afb(1)=ball1;afb(2)=ball2;afb(3)=minwallInd;Cnum=3;
    for i=1:comb
        xx=colliBall(i,1);yx=colliBall(i,2);
        dist=norm(pos(xx,:)-pos(yx,:));
        if dist<moleculed
            Cnum=Cnum+2;
            afb(Cnum-1)=xx;afb(Cnum)=yx;
            r12unit=(pos(yx,:)-pos(xx,:))/(norm(pos(yx,:)-pos(xx,:)));
            cmnterm=dot((vel(xx,:)-vel(yx,:)),r12unit)*r12unit;
            vb1=vel(xx,:)-cmnterm;vel(xx,:)=vb1;
            vb2=vel(yx,:)+cmnterm;vel(yx,:)=vb2;
            speed(xx)=norm(vel(xx,:));speed(yx)=norm(vel(yx,:));
            posxx=pos(xx,:)+vb1*(moleculed-dist)/speed(xx);pos(xx,:)=posxx;
            posyx=pos(yx,:)+vb2*(moleculed-dist)/speed(yx);pos(yx,:)=posyx;
            ancol=ancol+1;
        end
    end
    volrec(2,:,:)=vel;posrec(2,:,:)=pos;spdrec(2,:,:)=speed;
    colliWall(:,3)=colliWall(:,3)-deltat;colliWall(:,3)=colliWall(:,3)-deltat;

end
ancols(2)=ancol;
for Simnum=3:(ncol+1)
    for i=1:comb
        x=colliBall(i,1);y=colliBall(i,2);
        c1= concern(x,afb,Cnum);c2=concern(y,afb,Cnum);
        if (c1 || c2)
            rel_pos=pos(x,:)-pos(y,:);
            rel_vel=vel(x,:)-vel(y,:);
            mark1=dot(rel_pos, rel_vel);%Dot Product of rel position and rel velocity
            velsq=dot(rel_vel, rel_vel);
            possq=dot(rel_pos, rel_pos);
            mark2=mark1^2-velsq*(possq-moleculed^2);
            if mark1<0 && mark2>=0
                th1=(-mark1+sqrt(mark2))/velsq; th2=(-mark1-sqrt(mark2))/velsq;
                colliBall(i,3)=compare(th1,th2);
            else
                colliBall(i,3)=inf;
            end
        end
    end
    Cnum=0;
    afb=(num+1)*ones(num,1);
    for j=1:num
            posx=pos(j,1);  posy=pos(j,2);  velx=vel(j,1);  vely=vel(j,2);
            tx=0;ty=0;wx=0;wy=0;
            if (velx>0)
                tx=(boxsize-posx-moleculed/2)/velx;wx=2;
            elseif(velx<0)
                tx=(0-posx+moleculed/2)/velx;wx=1;
            elseif(velx==0)
                tx=inf;wx=5;
            end
            if (vely>0)
                ty=(boxsize-posy-moleculed/2)/vely;wy=4;
            elseif(vely<0)
                ty=(0-posy+moleculed/2)/vely;wy=3;
            elseif(vely==0)
                ty=inf;wy=5;
            end
            if(tx<0)
                tx=0;
            end
            if(ty<0)
                ty=0;
            end
            if (tx<ty)
                colliWall(j,2)=wx;colliWall(j,3)=tx;
            elseif(tx>ty)
                colliWall(j,2)=wy; colliWall(j,3)=ty;
            elseif(tx~=inf && tx==ty)
                colliWall(j,2)=6;colliWall(j,3)=ty;%code 6 collision with both wall together
            else
                colliWall(j,2)=5;colliWall(j,3)=inf;
            end
   end
   temporarystorage1=colliBall(:,3);temporarystorage2=colliWall(:,3);
   [mint_ball, minBallInd] = min(temporarystorage1);
   [mint_wall, minwallInd] = min(temporarystorage2);
   ancol=ancol+1;
   if mint_wall<mint_ball%hit wall first
       wallid=colliWall(minwallInd,2);deltat=mint_wall;t=t+deltat;
       time(Simnum)=t;dts(Simnum-1)=deltat;
       pos=pos+vel*deltat;
   if wallid <=2%x direction
       vel(minwallInd,1)=-vel(minwallInd,1);
   elseif wallid>2 && wallid <=4%y direction
       vel(minwallInd,2)=-vel(minwallInd,2);
   elseif wallid==6%in the corner
       vel(minwallInd,1)=-vel(minwallInd,1);
       vel(minwallInd,2)=-vel(minwallInd,2);
   end
   afb(1)=minwallInd;Cnum=1;
   for i=1:comb
       xx=colliBall(i,1);yx=colliBall(i,2);
       dist=norm(pos(xx,:)-pos(yx,:));
       if dist<moleculed
           Cnum=Cnum+2;
           afb(Cnum-1)=xx;afb(Cnum)=yx;
           r12unit=(pos(yx,:)-pos(xx,:))/(norm(pos(yx,:)-pos(xx,:)));
           cmnterm=dot((vel(xx,:)-vel(yx,:)),r12unit)*r12unit;
           vb1=vel(xx,:)-cmnterm;vel(xx,:)=vb1;
           vb2=vel(yx,:)+cmnterm;vel(yx,:)=vb2;
           speed(xx)=norm(vel(xx,:));speed(yx)=norm(vel(yx,:));
           posxx=pos(xx,:)+vb1*(moleculed-dist)/speed(xx);pos(xx,:)=posxx;
           posyx=pos(yx,:)+vb2*(moleculed-dist)/speed(yx);pos(yx,:)=posyx;
           ancol=ancol+1;
       end
   end
   volrec(Simnum,:,:)=vel;posrec(Simnum,:,:)=pos;spdrec(Simnum,:,:)=speed;
   colliWall(:,3)=colliWall(:,3)-deltat;colliWall(:,3)=colliWall(:,3)-deltat;
%% Collision with Ball Velocity Analysis 
   elseif mint_ball<mint_wall  %hit ball first
       ball1=colliBall(minBallInd,1); ball2=colliBall(minBallInd,2);
       deltat=mint_ball; t=t+deltat;dts(Simnum-1)=deltat;
       time(Simnum)=t;
       % Simple procedure found in Halie
       r12unit=(pos(ball2,:)-pos(ball1,:))/(norm(pos(ball2,:)-pos(ball1,:)));
       cmnterm=dot((vel(ball1,:)-vel(ball2,:)),r12unit)*r12unit;
       vb1=vel(ball1,:)-cmnterm;vel(ball1,:)=vb1;
       vb2=vel(ball2,:)+cmnterm;vel(ball2,:)=vb2;
       speed(ball1)=norm(vel(ball1,:));speed(ball2)=norm(vel(ball2,:));
       pos=pos+vel*deltat;
        afb(1)=ball1;afb(2)=ball2;Cnum=2;
        for i=1:comb
            xx=colliBall(i,1);yx=colliBall(i,2);
            dist=norm(pos(xx,:)-pos(yx,:));
            if dist<moleculed
                Cnum=Cnum+2;
                afb(Cnum-1)=xx;afb(Cnum)=yx;
                r12unit=(pos(yx,:)-pos(xx,:))/(norm(pos(yx,:)-pos(xx,:)));
                cmnterm=dot((vel(xx,:)-vel(yx,:)),r12unit)*r12unit;
                vb1=vel(xx,:)-cmnterm;vel(xx,:)=vb1;
                vb2=vel(yx,:)+cmnterm;vel(yx,:)=vb2;
                speed(xx)=norm(vel(xx,:));speed(yx)=norm(vel(yx,:));
                posxx=pos(xx,:)+vb1*(moleculed-dist)/speed(xx);pos(xx,:)=posxx;
                posyx=pos(yx,:)+vb2*(moleculed-dist)/speed(yx);pos(yx,:)=posyx;
                ancol=ancol+1;
            end
        end
       volrec(Simnum,:,:)=vel;posrec(Simnum,:,:)=pos;spdrec(Simnum,:,:)=speed;
       colliWall(:,3)=colliWall(:,3)-deltat;colliWall(:,3)=colliWall(:,3)-deltat;

%% Collision with Wall and Ball Simultaneous Velocity Analysis 
    else%hit wall and ball simultaneously
        wallid=colliWall(minwallInd,2);
        deltat=colliWall(minwallInd,3);
        t=t+deltat;time(Simnum)=t;dts(Simnum-1)=deltat;
        if wallid <=2
            vel(minwallInd,1)=-vel(minwallInd,1);
        elseif wallid>2 && wallid <=4
            vel(minwallInd,2)=-vel(minwallInd,2);
        elseif wallid==6%in the corner
            vel(minwallInd,1)=-vel(minwallInd,1);
            vel(minwallInd,2)=-vel(minwallInd,2);
        end
        ball1=colliBall(minBallInd,1); ball2=colliBall(minBallInd,2);
        % Simple procedure found in Halie
        r12unit=(pos(ball2,:)-pos(ball1,:))/(norm(pos(ball2,:)-pos(ball1,:)));
        cmnterm=dot((vel(ball1,:)-vel(ball2,:)),r12unit)*r12unit;
        vb1=vel(ball1,:)-cmnterm;vel(ball1,:)=vb1;
        vb2=vel(ball2,:)+cmnterm;vel(ball2,:)=vb2;
        speed(ball1)=norm(vel(ball1,:));speed(ball2)=norm(vel(ball2,:));
        pos=pos+vel*deltat;
        afb(1)=ball1;afb(2)=ball2;afb(3)=minwallInd;Cnum=3;
        for i=1:comb
            xx=colliBall(i,1);yx=colliBall(i,2);
            dist=norm(pos(xx,:)-pos(yx,:));
            if dist<moleculed
                Cnum=Cnum+2;
                afb(Cnum-1)=xx;afb(Cnum)=yx;
                r12unit=(pos(yx,:)-pos(xx,:))/(norm(pos(yx,:)-pos(xx,:)));
                cmnterm=dot((vel(xx,:)-vel(yx,:)),r12unit)*r12unit;
                vb1=vel(xx,:)-cmnterm;vel(xx,:)=vb1;
                vb2=vel(yx,:)+cmnterm;vel(yx,:)=vb2;
                speed(xx)=norm(vel(xx,:));speed(yx)=norm(vel(yx,:));
                posxx=pos(xx,:)+vb1*(moleculed-dist)/speed(xx);pos(xx,:)=posxx;
                posyx=pos(yx,:)+vb2*(moleculed-dist)/speed(yx);pos(yx,:)=posyx;
                ancol=ancol+1;
            end
        end
        volrec(Simnum,:,:)=vel;posrec(Simnum,:,:)=pos;spdrec(Simnum,:,:)=speed;
        colliWall(:,3)=colliWall(:,3)-deltat;colliWall(:,3)=colliWall(:,3)-deltat;
   end
   ancols(Simnum)=ancol;
end
        



%% Final State Plot
min(dts)
temp=1/2;
figure('Color', 'white','units','normalized','outerposition',[0 0 1 1])
p1=subplot(2,2,1);
h1 = plot(posrec(ncol+1,:,1),posrec(ncol+1,:,2),...
	'o','Color','red','MarkerFaceColor','red','MarkerSize',2);
xlim([0 21]);ylim([0 21]);interval=0.05;
edges=(-5:interval:5);edgeb=(0:interval:5);
coe1=(sqrt(molecuelm/(2*pi*temp*bconst)));%probability histogram vx vy
coe3=-(molecuelm/(2*temp*bconst));%exp stuff
coe4=(molecuelm/(temp*bconst));%probability histogram speed
curve1=zeros(length(edges),1);
for i=1:length(edges)
    curve1(i)=coe1*exp(coe3*edges(i)^2);
end
curve3=zeros(length(edgeb),1);
for i=1:length(edgeb)
    curve3(i)=coe4*edgeb(i)*exp(coe3*edgeb(i)^2);
end
simtitle=strcat('Simulation at time   ::   ', string(time(ncol+1)),'    Number of collisions   ::', string(ancol));
p1.XLabel.String = simtitle;
q2=subplot(2,2,2);
g =histogram(volrec(ncol+1,:,1),edges,'Normalization','pdf'); hold on;
cp1=plot(edges, curve1); hold on;
q2.XLabel.String = 'Vx (u/t)';q2.YLabel.String = 'Distribution';
p3=subplot(2,2,3);
k1 =histogram(volrec(ncol+1,:,2),edges,'Normalization','pdf');hold on;
cp2=plot(edges, curve1); hold on;
p3.XLabel.String = 'Vy (u/t)';p3.YLabel.String = 'Distribution';
p4=subplot(2,2,4);
l1=histogram(spdrec(ncol+1,:,1),edgeb,'Normalization','pdf');hold on;
cp3=plot(edgeb, curve3); hold on;
p4.XLabel.String = 'Speed (u/t)';p4.YLabel.String = 'Distribution';
%% Animated Distribution and particle movement
figure('Color', 'white','units','normalized','outerposition',[0 0 1 1])
p1=subplot(2,2,1);
h = plot(posrec(1,:,1),posrec(1,:,2),...
	'o','Color','red','MarkerFaceColor','red','MarkerSize',2);
xlim([0 21]);ylim([0 21]);interval=0.05;
simtitle=strcat('Simulation at time   ::   ', string(time(1)),'    Number of collisions   ::', string(ancols(1)));
p1.XLabel.String = simtitle;
q2=subplot(2,2,2);
g =histogram(volrec(1,:,1),edges,'Normalization','pdf');hold on;
cp1=plot(edges, curve1); hold on;
q2.XLabel.String = 'Vx (u/t)';q2.YLabel.String = 'Distribution';
p3=subplot(2,2,3);
k =histogram(volrec(1,:,2),edges,'Normalization','pdf');hold on;
cp2=plot(edges, curve1); hold on;
p3.XLabel.String = 'Vy (u/t)';p3.YLabel.String = 'Distribution';
p4=subplot(2,2,4);
l=histogram(spdrec(1,:,1),edgeb,'Normalization','pdf');hold on;
cp3=plot(edgeb, curve3); hold on;
p4.XLabel.String = 'Speed (u/t)';p4.YLabel.String = 'Distribution';
figtitle='Initial Condition';
saveas(gcf,figtitle,'png')
for i=2:ncol+1
    pause(0.002);
    h.XData=posrec(i,:,1);
    h.YData=posrec(i,:,2);
    simtitle=strcat('Simulation at time  ::   ', string(time(i)),'    Number of collisions   ::', string(ancols(i)));
    p1.XLabel.String = simtitle;
    subplot(2,2,2)
    hold off;
    g =histogram(volrec(i,:,1),edges,'Normalization','pdf');hold on;
    cp1=plot(edges, curve1); hold on;
    q2.XLabel.String = 'Vx (u/t)';
    q2.YLabel.String = 'Distribution';
    subplot(2,2,3)
    hold off;
    k =histogram(volrec(i,:,2),edges,'Normalization','pdf');hold on;
    cp2=plot(edges, curve1); hold on;
    p3.XLabel.String = 'Vy (u/t)';
    p3.YLabel.String = 'Distribution';
    subplot(2,2,4)
    hold off;
    l=histogram(spdrec(i,:,1),edgeb,'Normalization','pdf');hold on;
    cp3=plot(edgeb, curve3); hold on;
    p4.XLabel.String = 'Speed (u/t)';
    p4.YLabel.String = 'Distribution';
    if mod((i-1),200)==0
        figtitle=strcat('collision number',string(i-1));
        saveas(gcf,figtitle,'png')
    end
end
figtitle=strcat('collision number',string(ancols(i)));
saveas(gcf,figtitle,'png')

function toupdate = concern(i,afb,Cnum)
    toupdate=false;
    for j=1:Cnum
        if(i==afb(j))
            toupdate=true;
            break
        end
    end
end
function approtime=compare(time1,time2)
    if (time1>=0 && time2>=0)
        approtime=min(time1,time2);
    elseif(time1<0 && time2>=0)
        approtime=time2;
    elseif(time2<0 && time1>=0)
        approtime=time1;
    else
        approtime=inf;
    end
end

%% Different Distribution Initialization
%determine range of velocity randomization for 
function intervals=tophat(rangepair, num)
%range pair n*m m=2, then it is tophats, each interval is ranked from
%if m=3, the third parameter should be population in the subrange
%smallest to largest
intervals=zeros(num,4);
if size(rangepair,2)==2 && size(rangepair,1)>1
    tempr=linspace(0,sum(rangepair(:,2)-rangepair(:,1)),num+1);
    intervals(:,1)=tempr(1:num)+rangepair(1,1);intervals(:,2)=tempr(2:num+1)+rangepair(1,1);
    for i=2:size(rangepair,1)
        decision=(intervals(:,2)>rangepair(i-1,2))+(intervals(:,1)>rangepair(i-1,2));
        [~, idx] = max((decision==1));
        temp=intervals(idx,2);intervals(idx,2)=rangepair(i-1,2);
        intervals(idx,3)=rangepair(i,1);
        intervals(idx,4)=rangepair(i,1)+temp-rangepair(i-1,2);
        intervals(:,1)=intervals(:,1)+(decision==2)*(rangepair(i,1)-rangepair(i-1,2));
        intervals(:,2)=intervals(:,2)+(decision==2)*(rangepair(i,1)-rangepair(i-1,2));
    end
elseif size(rangepair,2)==2 && size(rangepair,1)==1
    tempr=linspace(0,(rangepair(2)-rangepair(1)),num+1);
    intervals(:,1)=tempr(1:num)+rangepair(1);intervals(:,2)=tempr(2:num+1)+rangepair(1);
elseif size(rangepair,2)==3
    step=1;
    for i=1:size(rangepair,1)
        tempr=linspace(rangepair(i,1),rangepair(i,2),rangepair(i,3)+1);
        intervals(step:rangepair(i,3)+step-1,1)=tempr(1:rangepair(i,3));
        intervals(step:rangepair(i,3)+step-1,2)=tempr(2:rangepair(i,3)+1);
        step=rangepair(i,3)+step;
    end
end
end
function intervals=triaglint(rangepair,tempr1,height)
%range pair[a,b,c] not for right angle triangles
    base1=(rangepair(2)-rangepair(1));base2=(rangepair(3)-rangepair(2));
    tempr=sqrt(tempr1*2*base1/height)+rangepair(1);
    num=length(tempr1)-1;
    intervals=zeros(num,2);
    intervals(:,1)=tempr(1:num);intervals(:,2)=tempr(2:num+1);
    decision=(intervals(:,2)>rangepair(2))+(intervals(:,1)>rangepair(2));
    [~, idx] = max((decision==1));
    intervals(idx,2)=rangepair(2);
    leftarea=base1*height/2;
    decision2=(tempr1-leftarea)>0;
    tempr2=decision2.*base2*height/2-(tempr1-leftarea).*decision2;
    tempr2=-sqrt(tempr2*2*base2/height)+rangepair(3)*decision2;
    intervals2=zeros(num,2);
    intervals2(:,1)=tempr2(1:num);intervals2(:,2)=tempr2(2:num+1);
    intervals(idx,2)=intervals2(idx,2);intervals2(idx,2)=0;intervals2(idx,1)=0;
    intervals(idx+1:num,:)=0;
    intervals=intervals+intervals2;
    intervals(num,2)=rangepair(3);
end
function intervals=thint(rangepair,tempr,height)
    num=length(tempr)-1;
    tempr=tempr/height+rangepair(1);
    intervals=zeros(num,2);
    intervals(:,1)=tempr(1:num);intervals(:,2)=tempr(2:num+1);
end
function intervals=triangle(rangepair, num) 
%range pair n*m m=3, then it is tophats, each interval is ranked from
%smallest to largest like y coordinate of the 3 points from smallest to
%largest
%if m=4 the forth coordinatewould be the population in that range
intervals=zeros(num,4);
if size(rangepair,2)==3 && size(rangepair,1)==1
    height=num*2/(rangepair(3)-rangepair(1));
    base1=(rangepair(2)-rangepair(1));base2=(rangepair(3)-rangepair(2));
    tempr=linspace(0,400,num+1);
    tempr=sqrt(tempr*2*base1/height)+rangepair(1);
    intervals(:,1)=tempr(1:num);intervals(:,2)=tempr(2:num+1);
    decision=(intervals(:,2)>rangepair(2))+(intervals(:,1)>rangepair(2));
    [~, idx] = max((decision==1));
    intervals(idx,2)=rangepair(2);
    stepsleft=num-idx+1;
    tempr2=linspace(stepsleft,0,stepsleft+1);
    tempr2=-sqrt(tempr2*2*base2/height)+rangepair(3);
    intervals2=zeros(num,4);
    intervals2(idx:num,1)=tempr2(1:length(tempr2)-1);intervals2(idx:num,2)=tempr2(1:length(tempr2)-1);
    intervals(idx,2)=intervals2(idx,2);intervals2(idx,2)=0;intervals2(idx,1)=0;
    intervals(idx+1:num,:)=0;intervals=intervals+intervals2;
elseif size(rangepair,2)==3 && size(rangepair,1)>1
    height=num*2/sum(rangepair(:,3)-rangepair(:,1));
    areas=height*(rangepair(:,3)-rangepair(:,1))/2;
    loc=1;
    st=1;
    for i=1:length(areas)
        atc=linspace(0,floor(areas(i)),floor(areas(i))+1);
        if st~=1
            atc=atc+st;
            atc=[0 atc];
        end
        if i~=length(areas)  && (floor(areas(i)))~=areas(i)
            atc=[atc areas(i)];
        end
        st=1-(areas(i)-floor(areas(i))-st);
        if st==2
            st=1;
        end
        subset=triaglint(rangepair(i,:),atc,height);
        if loc==1
            intervals(loc:length(subset),1)=subset(:,1);
            intervals(loc:length(subset),2)=subset(:,2);
            loc=length(subset);
        else
            if st~=1
                intervals(loc,3)=subset(1,1);
                intervals(loc,4)=subset(1,2);
                intervals(loc+1:length(subset)+loc-1,1)=subset(2:length(subset),1);
                intervals(loc+1:length(subset)+loc-1,2)=subset(2:length(subset),2);
            else
                intervals(loc+1:length(subset)+loc,1)=subset(1:length(subset),1);
                intervals(loc+1:length(subset)+loc,2)=subset(1:length(subset),2);
            end
        end
    end
    st=0;
elseif size(rangepair,2)==4 % the forth parameter is the number of particle in the range
        step=1;
    for i=1:size(rangepair,1)
        tempr=0:1:rangepair(i,4);
        height=rangepair(i,4)*2/(rangepair(i,3)-rangepair(i,1));
        subset=triaglint(rangepair(i,1:3),tempr,height);
        intervals(step:rangepair(i,4)+step-1,1)=subset(:,1);
        intervals(step:rangepair(i,4)+step-1,2)=subset(:,2);
        step=rangepair(i,4)+step;
    end
end
end
function intervals=trapezoid(rangepair, num)
%range pair n*m m=4, then it is tophats, each interval is ranked from
%smallest to largest like y coordinate of the 4 points from smallest to
%largest
%if m=5 the fifth coordinatewould be the population in that range
intervals=zeros(num,4);    
if size(rangepair,2)==4 && size(rangepair,1)==1
    height=num/((rangepair(4)-rangepair(1)+rangepair(3)-rangepair(2))/2);
    a1=height*(rangepair(2)-rangepair(1))/2;
    a2=height*(rangepair(3)-rangepair(2));
    a3=height*(rangepair(4)-rangepair(3))/2;
    step=1;
    st=0;
    if a1~=0
        atc=linspace(0,floor(a1),floor(a1)+1);
        if (floor(a1))~=a1
            atc=[atc a1];
        end
        st=1-(a1-floor(a1)-st);
        subset1=rtriaglint(rangepair(1:2),atc,height);
        intervals(step:length(subset1)+step-1,1)=subset1(:,1);
        intervals(step:length(subset1)+step-1,2)=subset1(:,2);
        step=length(subset1);
    end
    atc=linspace(0,floor(a2),floor(a2)+1);
    if st==2
        st=1;
    end
    if st~=1
        atc=atc+st;
        atc=[0 atc];
    end
    if (floor(a2))~=a2
        atc=[atc a2];
    end
    subset2=thint(rangepair(2:3),atc,height);
    if st~=1
        intervals(step,2)=subset2(1,2);
        intervals(step+1:length(subset2)+step-1,1)=subset2(2:length(subset2),1);
        intervals(step+1:length(subset2)+step-1,2)=subset2(2:length(subset2),2);
        
    else
        st=1;
        intervals(step+1:length(subset2)+step,1)=subset2(1:length(subset2),1);
        intervals(step+1:length(subset2)+step,2)=subset2(1:length(subset2),2);
    end
    step=length(subset2)+step;
    st=1-(a2-floor(a2)-st);
    if st==2
        st=1;
    end
    if a3~=0
        atc=linspace(0,floor(a3),floor(a3)+1);
        if st~=1
            atc=atc+st;
            atc=[0 atc];
        end
        if (floor(a3))~=a3
            atc=[atc a3];
        end
        subset3=ltrianglint(rangepair(3:4),atc,height);
        if st~=1
            intervals(step,2)=subset3(1,2);
            intervals(step+1:length(subset3)+step-1,1)=subset3(2:length(subset3),1);
            intervals(step+1:length(subset3)+step-1,2)=subset3(2:length(subset3),2);
        else
            st=1;
            intervals(step+1:length(subset3)+step,1)=subset3(1:length(subset3),1);
            intervals(step+1:length(subset3)+step,2)=subset3(1:length(subset3),2);
        end
    end
    step=length(subset3)+step;
    st=1-(a3-floor(a3)-st);
elseif size(rangepair,2)==4 && size(rangepair,1)>1
        height=num/sum((rangepair(:,4)-rangepair(:,1)+rangepair(:,3)-rangepair(:,2))/2);
        st=1;step=1;
        for i=1:size(rangepair,1) 
            a1=height*(rangepair(i,2)-rangepair(i,1))/2;
            a2=height*(rangepair(i,3)-rangepair(i,2));
            a3=height*(rangepair(i,4)-rangepair(i,3))/2;
            if a1~=0
                atc=linspace(0,floor(a1),floor(a1)+1);
                if (floor(a1))~=a1
                    atc=[atc a1];
                end
                subset1=rtriaglint(rangepair(i,1:2),atc,height);
                if i==1
                    intervals(step:length(subset1)+step-1,1)=subset1(:,1);
                    intervals(step:length(subset1)+step-1,2)=subset1(:,2);
                    step=length(subset1);
                else
                    if st~=1
                        intervals(step,3)=subset1(1,1);intervals(step,4)=subset1(1,2);
                        intervals(step+1:length(subset1)+step-1,1)=subset1(2:length(subset1),1);
                        intervals(step+1:length(subset1)+step-1,2)=subset1(2:length(subset1),2);
                    else
                        intervals(step+1:length(subset1)+step,1)=subset1(1:length(subset1),1);
                        intervals(step+1:length(subset1)+step,2)=subset1(1:length(subset1),2);
                    end
                    step=length(subset1)+step;
                end
            end
            st=1-(a1-floor(a1)-st);
            atc=linspace(0,floor(a2),floor(a2)+1);
            if st==2
                st=1;
            end
            if st~=1
                atc=atc+st;
                atc=[0 atc];
            end
            if (floor(a2))~=a2
                atc=[atc a2];
            end
            subset2=thint(rangepair(i,2:3),atc,height);
            if st~=1
                intervals(step,2)=subset2(1,2);
                intervals(step+1:length(subset2)+step-1,1)=subset2(2:length(subset2),1);
                intervals(step+1:length(subset2)+step-1,2)=subset2(2:length(subset2),2);
            else
                intervals(step+1:length(subset2)+step,1)=subset2(1:length(subset2),1);
                intervals(step+1:length(subset2)+step,2)=subset2(1:length(subset2),2);
            end
            step=length(subset2)+step;
            st=1-(a2-floor(a2)-st);
            if st==2
                st=1;
            end
            if a3~=0
                atc=linspace(0,floor(a3),floor(a3)+1);
                if st~=1
                    atc=atc+st;
                    atc=[0 atc];
                end
                if (floor(a3))~=a3
                    atc=[atc a3];
                end
                subset3=ltrianglint(rangepair(i,3:4),atc,height);
                if st~=1
                    intervals(step,2)=subset3(1,2);
                    intervals(step+1:length(subset3)+step-1,1)=subset3(2:length(subset3),1);
                    intervals(step+1:length(subset3)+step-1,2)=subset3(2:length(subset3),2);
                else
                    st=1;
                    intervals(step+1:length(subset3)+step,1)=subset3(1:length(subset3),1);
                    intervals(step+1:length(subset3)+step,2)=subset3(1:length(subset3),2);
                end
            end
            step=length(subset3)+step;
            st=1-(a3-floor(a3)-st);
            if st==2
                st=1;
            end
        end
elseif size(rangepair,2)==5
        height=num/sum((rangepair(:,4)-rangepair(:,1)+rangepair(:,3)-rangepair(:,2))/2);
        step=1;
        for i=1:size(rangepair,1) 
            st=1;
            height=rangepair(i,5)/((rangepair(i,4)-rangepair(i,1)+rangepair(i,3)-rangepair(i,2))/2);
            a1=height*(rangepair(i,2)-rangepair(i,1))/2;
            a2=height*(rangepair(i,3)-rangepair(i,2));
            a3=height*(rangepair(i,4)-rangepair(i,3))/2;
            if a1~=0
                atc=linspace(0,floor(a1),floor(a1)+1);
                if (floor(a1))~=a1
                    atc=[atc a1];
                end
                subset1=rtriaglint(rangepair(i,1:2),atc,height);
                if i==1
                    intervals(step:length(subset1)+step-1,1)=subset1(:,1);
                    intervals(step:length(subset1)+step-1,2)=subset1(:,2);
                    step=length(subset1);
                else
                    intervals(step+1:length(subset1)+step,1)=subset1(1:length(subset1),1);
                    intervals(step+1:length(subset1)+step,2)=subset1(1:length(subset1),2);
                    step=length(subset1)+step;
                end
            end
            st=1-(a1-floor(a1)-st);
            atc=linspace(0,floor(a2),floor(a2)+1);
            if st==2
                st=1;
            end
            if st~=1
                atc=atc+st;
                atc=[0 atc];
            end
            if (floor(a2))~=a2
                atc=[atc a2];
            end
            subset2=thint(rangepair(i,2:3),atc,height);
            if st~=1
                intervals(step,2)=subset2(1,2);
                intervals(step+1:length(subset2)+step-1,1)=subset2(2:length(subset2),1);
                intervals(step+1:length(subset2)+step-1,2)=subset2(2:length(subset2),2);
            else
                intervals(step+1:length(subset2)+step,1)=subset2(1:length(subset2),1);
                intervals(step+1:length(subset2)+step,2)=subset2(1:length(subset2),2);
            end
            step=length(subset2)+step;
            st=1-(a2-floor(a2)-st);
            if st==2
                st=1;
            end
            if a3~=0
                atc=linspace(0,floor(a3),floor(a3)+1);
                if st~=1
                    atc=atc+st;
                    atc=[0 atc];
                end
                if (floor(a3))~=a3
                    atc=[atc a3];
                end
                subset3=ltrianglint(rangepair(i,3:4),atc,height);
                if st~=1
                    intervals(step,2)=subset3(1,2);
                    intervals(step+1:length(subset3)+step-1,1)=subset3(2:length(subset3),1);
                    intervals(step+1:length(subset3)+step-1,2)=subset3(2:length(subset3),2);
                else
                    st=1;
                    intervals(step+1:length(subset3)+step,1)=subset3(1:length(subset3),1);
                    intervals(step+1:length(subset3)+step,2)=subset3(1:length(subset3),2);
                end
            end
            step=length(subset3)+step;
        end        
end
end
function intervals=rtriaglint(rangepair,tempr1,height)
%range pair [a,b]for triangles with right angle on the right
    base1=(rangepair(2)-rangepair(1));
    tempr=sqrt(tempr1*2*base1/height)+rangepair(1);
    num=length(tempr1)-1;
    intervals=zeros(num,2);
    intervals(:,1)=tempr(1:num);intervals(:,2)=tempr(2:num+1);
end
function intervals=ltrianglint(rangepair,tempr1,height)
%range pair [a b] for triangle with right angle on the left
    base2=(rangepair(2)-rangepair(1));
    totarea=base2*height/2;
    tempr2=totarea-tempr1;
    num=length(tempr1)-1;
    tempr2=-sqrt(tempr2*2*base2/height)+rangepair(2);
    intervals=zeros(num,2);
    intervals(:,1)=tempr2(1:num);intervals(:,2)=tempr2(2:num+1);
end
function xvvel=impulse(rangepair, num)
% rangepair [a,b]*n  is the speed b is the population with the speed
% if range pair is just a number,  then all population should have the same
% speed
% with initial speed directly output velocitys
xvvel=zeros(num,2);
spd=zeros(num,1);
inidir=rand(1,num)*2*pi;
if size(rangepair,1)==2
    step=1;
    for i=1:size(rangepair,2)
        spd(step:rangepair(i,2)+step-1)=ones(1,rangepair(i,2))*rangepair(i,1);
        step=rangepair(i,2)+step;
    end
    xvvel(:,1)=spd.*cos(inidir)';xvvel(:,2)=spd.*sin(inidir)';
elseif size(rangepair,1)==1
    xvvel(:,1)=rangepair*cos(inidir);xvvel(:,2)=rangepair*sin(inidir);
end
end
function xvvel=randomdis(range, num)
%range [a b]
inidir=rand(1,num)*2*pi;
spd=rand(num,1)*(range(2)-range(1))+range(1);
xvvel(:,1)=spd.*cos(inidir)';xvvel(:,2)=spd.*sin(inidir)';    
end
function xvvel=settemp(vel,temp,num)
xvvel=vel;
sumvel=sum(vel.^2,"all");
curtemp=sumvel/num;
factor=sqrt(temp/curtemp);
xvvel=vel*factor;
end
function xvvel=radinrange(rangex,rangey,num)
xloc=randperm(num);
yloc=randperm(num);
vxpro=rand(num,1);
vypro=rand(num,1);
vxr1=rangex(:,2)-rangex(:,1);
vxr2=rangex(:,4)-rangex(:,3);
vxrtot=vxr2+vxr1;
dvx=vxrtot.*vxpro;
ifrx=(dvx<=vxr1);
vx=ifrx.*(dvx+rangex(:,1))+(ifrx==0).*(dvx-vxr1+rangex(:,3));
vyr1=rangey(:,2)-rangey(:,1);
vyr2=rangey(:,4)-rangey(:,3);
vyrtot=vyr2+vyr1;
dvy=vyrtot.*vxpro;
ifry=(dvy<=vyr1);
vy=ifry.*(dvy+rangey(:,1))+(ifry==0).*(dvy-vyr1+rangey(:,3));
xvvel=zeros(num,2);
for i=1:num
    xvvel(i,1)=vx(xloc(i));
    xvvel(i,2)=vy(yloc(i));    
end
end

%% Collision Integral
