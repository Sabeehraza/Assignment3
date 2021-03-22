% part 1 from assignment 1
%using all the infromation from assignment 1 apart from where it is
%indicated for us to change 
set(0,'DefaultFigureWindowStyle','docked')
clear all;
%Given Paramters
T = 300;  %Temp in Kelvin
m = 9.10938356e-31; %mass of the Electron
Cm = 0.26*m; % effective Mass of Electrons
q = -1.60217662e-19; % eharge of an electron
width = 200e-9;
len = 100e-9;
%chaning Vx applied at x to 0.1v given in 1a
Vx = 0.1; 
Vy = 0;
%electron concentration given in part 1d
elecconc = 1e15*100^2;
% thermal Velocity
K = 1.38064852e-23;  %Boltzmann Constatnt
Vt = sqrt(2*K*T/Cm);
elecpop = 30000;
elecnum = 50;
Time_Step = len/Vt/100;
iterations = 200;
% The scattering probabily is given by
scatter_p = 1 - exp(-Time_Step/0.2e-12);
%Using Maxewell-Boltzmann Distribution to Generate random velocities
Distr_MB = makedist('Normal', 0, sqrt(K*T/Cm));
ap = 0; %animation plot 
scatter_p = 1 - exp(-Time_Step/0.2e-12);
Xr = 2e-7;
Yr = 1e-7;
EField = Vx / Xr;
EF = q * EField;
AccEle = EF/Cm;%electron acc
% Calculating the Electric Field 
Electric_Field = Vx/width; 
fprintf(' The Electric field experienced by the electrons is %i\n',Electric_Field);
Force = Electric_Field*q;
Accelaration = Force/Cm; 
ElecFieldY = Vy/len
ForceY = q*ElecFieldY
DX = Force*Time_Step/Cm;  %velocity of particle 
DY = ForceY*Time_Step/Cm;   %velocity of particle 
DX = DX.*ones(elecpop,1);    %velocity of particle 
DY = DY.*ones(elecpop,1);   %velocity of particle 
% Defining the Specular conditions upper and lower boundary;
Specupper = 0;
Speclower = 0;
%randomlyasigning posto a particle 
%giving the particle an x pos
posx = zeros(elecpop, 4);

for ni = 1:elecpop
    
    theta = rand*2*pi;
    
    posx(ni,:) = [width*rand len*rand random(Distr_MB) random(Distr_MB)];
    
end

%traking the path of the particle
path = zeros(iterations, elecnum*2);

%traking the temprature
Temp = zeros(iterations,1);

%current density 
CurrDensity = zeros(iterations,2);

% Path of electron
figure(1);
plot([],[]);
title('Electron path, '); 
xlabel('X pos');
ylabel('Y pos');

%current v time
figure(2);
current_plot =  animatedline('Color','b','LineWidth',5); 
title('Current v Time, ');
xlabel('Time');
ylabel('Current');
grid on;

% Simulating the Monte Carlo Model

for ni = 1:iterations 
    
    %Utilzing the velocities of the electrons to determine positions
    posx(:,3) = posx(:,3) + DX;
    posx(:,4) = posx(:,4) + DY;

     %Boundary Conditions 
    posx(:,1:2) = posx(:,1:2) + Time_Step.*posx(:,3:4);

    nk = posx(:,1) > width;
    posx(nk,1) = posx(nk,1) - width;
    
    nk = posx(:,1) < 0;
    posx(nk,1) = posx(nk,1) + width;
    
    nk = posx(:,2) > len;

    % At the upper of Region Boundary 
    if(Specupper)    
        
        posx(nk,2) = 2*len - posx(nk,2);
        posx(nk,4) = -posx(nk,4);
        
    else 
        
        posx(nk,2) = 100e-9;      
        Z = sqrt(posx(nk,3).^2 + posx(nk,4).^2);
        
        theta = rand([sum(nk),1])*2*pi;
        posx(nk,3) = Z.*cos(theta);
        posx(nk,4) = -abs(Z.*sin(theta));
      
    end
  
    nk = posx(:,2) < 0;
    
    % At the lower region boundary 
      if(Speclower)
                                                          
        posx(nk,2) = -posx(nk,2);
        posx(nk,4) = -posx(nk,4);
        
      else  
                                                           
        posx(nk,2) = 0;                                       
        Z = sqrt(posx(nk,3).^2 + posx(nk,4).^2);
        
        theta = rand([sum(nk),1])*2*pi;
        posx(nk,3) = Z.*cos(theta);
        posx(nk,4) = abs(Z.*sin(theta));
          
    end
 
    nk = rand(elecpop, 1) < scatter_p;           
    posx(nk,3:4) = random(Distr_MB, [sum(nk),2]);  
    
    %Calculating and Recording the Temperature of the electrons
    Temp(ni) = (sum(posx(:,3).^2) + sum(posx(:,4).^2))*Cm/K/2/elecpop;           

    for nk=1:elecnum
        
        %Storing the positions of the electrons on the path 
        path(ni, (2*nk):(2*nk+1)) = posx(nk, 1:2);
        
    end
  
    CurrDensity(ni, 1) = q.*elecconc.*mean(posx(:,3));           
    
    CurrDensity(ni, 2) = q.*elecconc.*mean(posx(:,4));

    addpoints(current_plot, Time_Step.*ni, CurrDensity(ni,1));           

    
    %Plotting the electron trajectories and temeperature variance
    % The animation updates for every 10 iterations
    if(ap && mod(ni,10) == 0)
        figure(1);
        hold off;
        plot(posx(1:elecnum,1), posx(1:elecnum,2), 'x');       
        hold on;
        title(' Electron path, '); 
        xlabel('X pos)');
        ylabel('Y pos');
        pause(0.05);
        
    end
    
end

%Final plotting of the trajectories after all the iterations have been
%completed

figure(1);
title('Electron path, ');
xlabel('X pos');
ylabel('Y pos');
grid on;
hold on;

%Storing the trajectory after completion 
for ni=1:elecnum
    
    plot(path(:,ni*2), path(:,ni*2+1), '+');
    
end

% Plot the Density Map 

%Utilizing the Hist3 command to generate the density map 
elecconc = hist3(posx(:,1:2),[200 100])';

 %Utilzing the bins from the historgram to plot the the Temperature and Density map 
histbins = 10;    

[xelec yelec] = meshgrid(round(-histbins/2):round(histbins/2), round(-histbins/2):round(histbins/2));
L=exp(-xelec.^2/(2*1^2)-yelec.^2/(2*1^2));
L=L./sum(L(:));

figure(3);
elecconc = conv2(elecconc,L,'same');
elecconc = elecconc/(len./size(elecconc,1)*width./size(elecconc,2));
surf(conv2(elecconc,L,'same'));
title('Electron Density, ');
xlabel('X position');
ylabel('Y position');

% Temperature Map

XTemp = zeros(ceil(width/1e-9),ceil(len/1e-9));
YTemp = zeros(ceil(width/1e-9),ceil(len/1e-9));
TotalTemp = zeros(ceil(width/1e-9),ceil(len/1e-9));


for ni=1:elecpop
   
    xelec = floor(posx(ni,1));
    yelec = floor(posx(ni,2));   
    if(xelec==0)
        xelec = 1;
    end
    if(yelec==0) 
        yelec= 1;
    end
    
    YTemp(xelec,yelec) = YTemp(xelec,yelec) + posx(ni,3)^2;
    XTemp(xelec,yelec) = XTemp(xelec,yelec) + posx(ni,4)^2;
    TotalTemp(xelec,yelec) = TotalTemp(xelec,yelec) + 1;
    
end


%temp can be calculated as follows 

temp = (XTemp + YTemp)*Cm/K/2/TotalTemp;
temp= temp;  

%Generate The temperature plot
histbins = 10;
[xelec , yelec] = meshgrid(round(-histbins/2):round(histbins/2), round(-histbins/2):round(histbins/2));
L=exp(-xelec.^2/(2*1^2)-yelec.^2/(2*1^2));
L=L/sum(L(:));
figure(4);
imagesc(conv2(temp,L,'same'));  
view(0,90)
title('Temperature Map, ');
xlabel('X pos');
ylabel('Y pos');

% Part 2 Using Assignment 2 

L = 2;
W = 3;
V0 = 1;    

dx = 0.2;
dy = 0.2;

nx = 100*L;  
ny = 100*W;

WBox = 0.4;     
LBox = 0.4;     

%sigma inside an outside of the box
sigin = 1;        
sigout = 1e-2;     

C = zeros(ny,nx);   %conductivity 

for j = 1:ny
    
    for g = 1:nx
        
        if(g >= nx*WBox && g <= nx-nx*WBox && (j >= ny-ny*LBox || j <= ny*LBox)) 
           
            % If inside the box, then sigma = 1e^-2
           
            C(j,g) = sigout;  
       
        else                                        
            %outside the box
            C(j,g) = sigin;
            
        end
   
    end

end



%  Creating the G matrix and B vector for the GV = F solution 

 G = sparse(nx*ny);
 F = zeros(nx*ny,1);

for g = 1:nx
    
    for j = 1:ny
        
       %node equation  
        n = j + (g - 1)*ny;
     
       %nodes at g and j
       
        nxm = j + (g - 2)*ny;
        nxp = j + g*ny;
        nym = (j - 1) + (g - 1)*ny;
        nyp = (j + 1) + (g - 1)*ny;

        if(g == 1 || g == nx) 
            
            
            G(n,n) = 1;       %Left Side Set Vx 
   
        elseif (j == 1)    %Evalutation at the bottom region 
            
            UY = (C(j,g)+C(j+1,g))/2;
            UX = (C(j,g)+C(j,g+1))/2;
            UXD = (C(j,g)+C(j,g-1))/2;
            
            G(n,n) = -(UY + UX + UXD);
            G(n,nyp) = UY;
            G(n,nxp) = UX;
            G(n,nxm) = UXD;
     
       % Evaluation at upper region   
       
       elseif (j == ny)
       
              YDYE = (C(j,g)+C(j-1,g))/2;
              UDXE = (C(j,g)+C(j,g+1))/2;
              DDXE = (C(j,g)+C(j,g-1))/2;
              
              G(n,n) = -(YDYE + UDXE + DDXE);
              G(n,nym) = YDYE;
              G(n,nxp) = UDXE;
              G(n,nxm) = DDXE;
       
        else 
   
           
            UY = (C(j,g)+C(j+1,g))/2;
            D_Y = (C(j,g)+C(j-1,g))/2;
            UX = (C(j,g)+C(j,g+1))/2;
            UXD = (C(j,g)+C(j,g-1))/2;
            
            
            G(n,n) = -(UY + D_Y + UX + UXD);
            G(n,nyp) = UY;
            G(n,nym) = D_Y;
            G(n,nxp) = UX;
            G(n,nxm) = UXD;
            
            
        end
        
        
    end
    
    
end





for g = 1:nx
    
    for j = 1:ny
     
        
        
        %Node Mapping Equation 
        n = j + (g - 1)*ny;
       
        
        if (g == 1) %Indicating a shift towards left so the value must be set equal to V_0
            
          
            F(n) = V0;
            
            
        end
        
        
        
    end
    
end

% Utilizing  GV = F to solve the equation 

V = G\F;


for g = 1:nx
    
    for j = 1:ny
        
        % Node mapping to put entries into the correct place
        n = j + (g - 1)*ny;
        
        Vmap(j,g) = V(n);
        
        
    end
end



% Plotting the Vx V across the region
figure(5)
surf(Vmap)

xlabel('L (um)')
ylabel('W (um)')
title({'Vx V(x,y) Plot, '})



[Ex,Ey] = gradient(-Vmap);   %Electric Field of the regions can be determined from using the gradient




% Plotting the Electric Field over the region 


figure(6)
quiver(Ex,Ey)
xlabel('L (um)')
ylabel('W (um)')
title({'Electric Field in the Region, '})



%Given Paramters
T = 300;  %Semicondctor Temperature 
Cm = 9.10938356e-31; %Rest mass of the Electron
Cm = 0.26*Cm; %Given Effective Mass of Electrons
q = -1.60217662e-19; % Charge on electron


%Given Nominal Dimensions of Semiconductor 200 nm x 100 nm
width = 200e-9;
len = 100e-9;

% Updated Conditions include the Vx being applied in x direction to be
% 0.1v
Vx = 0.1; 

% No changes to the Vx being applied in Y direction
Vy = 0;


% The electron concentration is given as
elecconc = 1e15*100^2;

% Calculation of Thermal Velocity
K = 1.38064852e-23;  %Boltzmann Constatnt
Vt = sqrt(2*K*T/Cm);


elecpop = 10000;
elecnum = 100;

% Setting the Step Size
Time_Step = len/Vt/100;
iterations = 200;

% The scattering probabily is given by
scatter_p = 1 - exp(-Time_Step/0.2e-12);

%Using Maxewell-Boltzmann Distribution to Generate random velocities
Distr_MB = makedist('Normal', 0, sqrt(K*T/Cm));



ap = 0;


% The scattering probabily is given by
scatter_p = 1 - exp(-Time_Step/0.2e-12);


% Calculating the Electric Field 
Electric_Field = Vx/width; 
fprintf(' The Electric field experienced by the electrons is %i\n',Electric_Field);

% Calculating the Force
Force = Electric_Field*q;
fprintf(' The Force experienced by the electrons is %i\n',Force);

% Calculating the Accelaration 
Accelaration = Force/Cm; 
fprintf(' The accelaration of the electrons is %i\n',Accelaration);

%The Electric Field and the Force in the Y direction calculated as 
ElecFieldY = Vy/len
ForceY = q*ElecFieldY


% Calculating velocity at each time step of the electrons
DX = Force*Time_Step/Cm;
DY = ForceY*Time_Step/Cm;
DX = DX.*ones(elecpop,1);
DY = DY.*ones(elecpop,1);


% Defining the Specular conditions at the Top and the Bottom boundary of
% the region
Specupper = 0;
Speclower = 0;



%Initialization of the posof Random Particles
%Setting up x and positions
posx = zeros(elecpop, 4);



%Keeping track of trajectories
path = zeros(iterations, elecnum*2);


%Recording the temeperatures
Temp = zeros(iterations,1);


st_size = 1e-9;
boxes = st_size.*[80 120 0 40; 80 120 60 100];
spec_boxes = [0 1];

%Initialization of the posof Random Particles
%Setting up x and positions
for ni = 1:elecpop
    
    theta = rand*2*pi;
    
    posx(ni,:) = [width*rand len*rand random(Distr_MB) random(Distr_MB)];
    
end



% Simulating the Monte Carlo Model

for ni = 1:iterations                        % Performs the simulation with Random Velocities with updating positions. 
    
    
    %Utilzing the velocities of the electrons to determine positions
    posx(:,3) = posx(:,3) + DX;
    posx(:,4) = posx(:,4) + DY;
    
    
     %Boundary Conditions 
    posx(:,1:2) = posx(:,1:2) + Time_Step.*posx(:,3:4);

    nk = posx(:,1) > width;
    posx(nk,1) = posx(nk,1) - width;
    
    nk = posx(:,1) < 0;
    posx(nk,1) = posx(nk,1) + width;
    
    nk = posx(:,2) > len;

    
    % At the Top of Region Boundary 
    if(Specupper)                                  % for specular condition 
        
        posx(nk,2) = 2*len - posx(nk,2);
        posx(nk,4) = -posx(nk,4);
        
        
    else 
        
        posx(nk,2) = 100e-9;                                        %Diffusive condition has been met
        Z = sqrt(posx(nk,3).^2 + posx(nk,4).^2);
        
        theta = rand([sum(nk),1])*2*pi;
        posx(nk,3) = Z.*cos(theta);
        posx(nk,4) = -abs(Z.*sin(theta));
        
        
    end
    
    
    nk = posx(:,2) < 0;
    
    
    % At the bottom of the region boundary 
      if(Speclower)
                                                          % for specular condition 
        posx(nk,2) = -posx(nk,2);
        posx(nk,4) = -posx(nk,4);
        
        
      else  
                                                           
        posx(nk,2) = 0;                                       %Diffusive condition has been met
        Z = sqrt(posx(nk,3).^2 + posx(nk,4).^2);
        
        
        theta = rand([sum(nk),1])*2*pi;
        posx(nk,3) = Z.*cos(theta);
        posx(nk,4) = abs(Z.*sin(theta));
        
        
    end
    
    
    %Figuring out if the particles have moved into to the box. Updates the
  %posand restores the location of the particles
    
    for nk= 1:elecnum
        bottle = box(posx(nk,1:2), boxes);
        
        
        % Checking for the collision with a box and determining the
        % location of the box of collision. 
        
        while(bottle ~= 0)
            
            dist_X = 0;                  %Finding and updating the X position
            
            X_updated = 0;
            
            
            if(posx(nk,3) > 0)
                
                dist_X = posx(nk,1) - boxes(bottle,1);
                X_updated = boxes(bottle,1);
                
            else
                
                dist_X = boxes(bottle,2) - posx(nk,1);
                X_updated = boxes(bottle,2);
                
                
            end

            dist_Y = 0;                  %Finding and updating the Y position
            Y_updated = 0;
            
            if(posx(nk,4) > 0)
                
                dist_Y = posx(nk,2) - boxes(bottle, 3);
                Y_updated = boxes(bottle, 3);
                
            else
                
                dist_Y = boxes(bottle, 4) - posx(nk,2);
                Y_updated = boxes(bottle, 4);
                
            end

            if(dist_X < dist_Y)
                
                posx(nk,1) = X_updated;
                
                if(~spec_boxes(bottle))
                    
                    sgn = -sign(posx(nk,3));
                    Z = sqrt(posx(nk,3).^2 + posx(nk,4).^2);
                    
                    theta = rand()*2*pi;
                    posx(nk,3) = sgn.*abs(Z.*cos(theta));
                    posx(nk,4) = Z.*sin(theta);
                    
                    
                else 
                    
                    %For specular condition
                    
                    posx(nk,3) = -posx(nk,3);
                    
                end
                
                
            else
                
                
                posx(nk,2) = Y_updated;
                if(~spec_boxes(bottle))
                    
                    sgn = -sign(posx(nk,4));
                    Z = sqrt(posx(nk,3).^2 + posx(nk,4).^2);
                    theta = rand()*2*pi;
                    
                    posx(nk,3) = Z.*cos(theta);
                    posx(nk,4) = sgn.*abs(Z.*sin(theta));
                    
                else 
                    
                    %For speuclar condition
                    
                    posx(nk,4) = -posx(nk,4);
                    
                end
            end
            

            bottle = box(posx(nk,1:2), boxes);
            
            
        end
        
    end
    
     nk = rand(elecpop, 1) < scatter_p;           % Random dsitribution of particles within the scattering probability limit
     posx(nk,3:4) = random(Distr_MB, [sum(nk),2]);          %Scattering particles using the exponential scattering probability 
    
    
    
    %Calculating and Recording the Temperature of the electrons
    Temp(ni) = (sum(posx(:,3).^2) + sum(posx(:,4).^2))*Cm/K/2/elecpop;           
    
   
    
    for nk=1:elecnum
        
        %Storing the positions of the electrons on the trajector 
        path(ni, (2*nk):(2*nk+1)) = posx(nk, 1:2);       
    end
  
    %Plotting the electron trajectories
    % The animation updates for every 10 iterations
    if(ap && mod(ni,10) == 0)
        figure(7);
        hold off;
        plot(posx(1:elecnum,1), posx(1:elecnum,2), 'o');       
        hold on;
        
        
          for nk=1:size(boxes,1)          %Plotting the rectangular boxes 
            
           plot([boxes(nk, 1) boxes(nk, 1) boxes(nk, 2) boxes(nk, 2) boxes(nk, 1)],...
               [boxes(nk, 3) boxes(nk, 4) boxes(nk, 4) boxes(nk, 3) boxes(nk, 3)], 'k-');
          
          end
        
        title('Electron Trajectories, '); 
        xlabel('X position)');
        ylabel('Y position');
        pause(0.05);
        
    end
    
end


%Final Plotting of the boxes after completion of iterations
for nk=1:size(boxes,1)
    
   plot([boxes(nk, 1) boxes(nk, 1) boxes(nk, 2) boxes(nk, 2) boxes(nk, 1)],...
       [boxes(nk, 3) boxes(nk, 4) boxes(nk, 4) boxes(nk, 3) boxes(nk, 3)], 'k-');
   
end

%Final plotting of the trajectories after all the iterations have been
%completed

figure(7);
title('Electron Trajectories, ');
xlabel('X position');
ylabel('Y position');

grid on;
hold on;

%Storing the trajectory after completion 
for ni=1:elecnum
    
    plot(path(:,ni*2), path(:,ni*2+1), '.');
    
end

% Part 3 Combining the Simulators to extract simple paramters

% a) 
% Plot the Density Map 

%Utilizing the Hist3 command to generate the density map 
elecconc = hist3(posx(:,1:2),[200 100])';

 %Utilzing the bins from the historgram to plot the the Temperature and Density map 
histbins = 10;    

[xelec yelec] = meshgrid(round(-histbins/2):round(histbins/2), round(-histbins/2):round(histbins/2));
L=exp(-xelec.^2/(2*1^2)-yelec.^2/(2*1^2));
L=L./sum(L(:));

figure(8);
elecconc = conv2(elecconc,L,'same');  
elecconc = elecconc/(len./size(elecconc,1)*width./size(elecconc,2));
surf(conv2(elecconc,L,'same'));
title('Electron Density at volatge of 0.8V , ');
xlabel('X pos');
ylabel('Y pos');

% Setting the Length and Width paramters of the rectangular region
%Utilizng the Length and width to configure the elements ranges
nx = 100*L;
ny = 100*W;

for  numIter = 1:5
     C = zeros(ny,nx); % Map of Conductivity 
     WBox = 0.4*(1+(numIter/20));      
     LBox = 0.4*(1+(numIter/20));      
   
    for j = 1:ny
        
       
         for g = 1:nx
           
            
            if(g >= nx*WBox && g <= nx-nx*WBox && (j >= ny-ny*LBox || j <= ny*LBox))               
                C(j,g) = sigout;
     
            else
 
                C(j,g) = sigin;
        
            end
            
        end
        
    end
    
    
    
    % creating a sparse G matrix
    G = sparse(nx*ny,nx*ny);
    F = zeros(nx*ny,1);

    for g = 1:nx
        
        for j = 1:ny

                                           
            n = j + (g - 1)*ny;              
            
            nxm = j + (g - 2)*ny;
            nxp = j + g*ny;
            nym = (j - 1) + (g - 1)*ny;
            nyp = (j + 1) + (g - 1)*ny;
  
            if(g == 1 || g == nx) 
                G(n,n) = 1;
      % at lower region 
            
            elseif (j == 1)   
                
                UY = (C(j,g)+C(j+1,g))/2;
                UX = (C(j,g)+C(j,g+1))/2;
                UXD = (C(j,g)+C(j,g-1))/2;

                G(n,n) = -(UY + UX + UXD);
                G(n,nyp) = UY;
                G(n,nxp) = UX;
                G(n,nxm) = UXD;
    
            %at upper region 
            
            elseif (j == ny) 

                D_Y = (C(j,g)+C(j-1,g))/2;
                UX = (C(j,g)+C(j,g+1))/2;
                UXD = (C(j,g)+C(j,g-1))/2;

                G(n,n) = -(D_Y + UX + UXD);
                G(n,nym) = D_Y;
                G(n,nxp) = UX;
                G(n,nxm) = UXD;     
                
            else  
                
                UY = (C(j,g)+C(j+1,g))/2;
                D_Y = (C(j,g)+C(j-1,g))/2;
                UX = (C(j,g)+C(j,g+1))/2;
                UXD = (C(j,g)+C(j,g-1))/2;

                G(n,n) = -(UY + D_Y + UX + UXD);
                G(n,nyp) = UY;
                G(n,nym) = D_Y;
                G(n,nxp) = UX;
                G(n,nxm) = UXD;
                

            end
            
        end
        
    end
    
    
   
    for g = 1:nx
        
        for j = 1:ny
            
          
            n = j + (g - 1)*ny;   
 
            if (g == 1) 
                
                F(n) = V0;
                
                
            end

        end
        
    end

    V = G\F;
 
    Vmap = 0;
    
    for g = 1:nx
        
        for j = 1:ny
             
            n = j + (g - 1)*ny;
     
            Vmap(j,g) = V(n);

    end

    [Ex,Ey] = gradient(-Vmap);  %using gradient from assignment 2
 
    boxSL(numIter) = sum(C(:,1).*Ex(:,1));
    boxSR(numIter) = sum(C(:,nx).*Ex(:,nx));
    
    
end

% c) 

figure(9)
plot(linspace(1,5,5),boxSR,'r:')
hold on
plot(linspace(1,5,5),boxSR,'r:')

hold off
title({'Current vBox Width, '})
xlabel('Box Size Changes')
ylabel({'I (A)'})

end
%function from my assignment 1 
function nbox = box(pos, boxes)
    nbox = 0;
    for i=1:size(boxes,1)
        if(pos(1) > boxes(i,1) && pos(1) < boxes(i,2) && pos(2) > boxes(i,3) && pos(2) < boxes(i,4))
            nbox = i;
            return;
        end
    end
end
