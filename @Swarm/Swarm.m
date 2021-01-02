classdef Swarm < handle
    % SWARM - This class represents an ensemble of dynamic agents of type
    % "Drone"
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Swarm general properties:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % drones:
    %           vector of Drone objects
    % nb_agents:
    %           size of the above vector
    % equivalent_drone:
    %           for path planner, drone at the barycenter ...
    %           of the swarm for command computations
    % pos_ned:

    properties
        drones % a vector of Drone objects
        nb_agents % size of the above vector
        equivalent_drone % for path planner, drone at the barycenter ...
                         % of the swarm for command computations
        algorithm SwarmAlgorithm
        collisions_history
        GS_required % is ground station required or not
        GS
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods
        %%%%%%%%%% Constructor %%%%%%%%%%%%
        function self = Swarm()
            self.drones = [];
            self.nb_agents = 0;
            self.collisions_history = [];
            self.GS=struct;
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function add_GS(self, mx,my, Az, El,frequency,Noise_power,bandwidth,SNR)
            % this function creates a Ground station array with Mx X My
            % elements
            self.GS.N0 = Noise_power;
            self.GS.pu = SNR;
            self.GS.bw = bandwidth;
            self.GS.Mx = mx;
            self.GS.My = my;
            ang= 0:0.1:pi;
            
            self.GS.c = physconst('LightSpeed');
            self.GS.freq = frequency;
            lambda = self.GS.c/self.GS.freq;
            delx = lambda/2;
            dely = lambda/2;
            X = zeros(1,mx*my);
            Y = zeros(1,mx*my);
            Z = zeros(1,mx*my);
            phi = zeros(1,mx*my);
            theta = zeros(1,mx*my);
            psi = zeros(1,mx*my);
            for p = 1:mx
                for q = 1:my
                    i = (q-1)*mx+p;
                    X(i) = (p-1)*delx;
                    Y(i) = (q-1)*dely;
                    Z(i) = 0;
                    phi(i)=ang(randi(length(ang)));
                    theta(i)=ang(randi(length(ang)));
                    psi(i)=ang(randi(length(ang)));
                end
            end
            self.GS.antenna_pos=[X;Y;Z];
            self.GS.antenna_ang=[phi;theta;psi];
            self.GS.sCA = phased.ConformalArray('ElementPosition',self.GS.antenna_pos,'ElementNormal',[Az;El],'Element',dipoleCrossed);
            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function g_kl=get_channelVector(self,freq)
            c = self.GS.c;
            lambda = c/freq;
            M = length(self.GS.antenna_pos(1,:));
            K = self.nb_agents;
            pathloss=self.get_pathloss(freq);
            d_kl = self.get_antenna_dist();
            g_kl=zeros(K,M);
            %h_kl=ones(K,M);
            h_kl=self.get_PLF(freq);
            for k=1:K
                for l=1:M                    
                    g_kl(k,l) = h_kl(k,l)*sqrt(pathloss(k,l))*exp(-1i*2*pi*d_kl(k,l)/lambda);
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %gives the uplink throughput
        function S=get_uplinkThroughput(self,Bc,V,fc,tau_dl,tau_ctrl)
            B=self.GS.bw;
            c=self.GS.c;
            Sk=self.get_uplink_capacity(fc);
            S=zeros(self.nb_agents,1);
            for k=1:self.nb_agents
                S(k)=B*(1-((2*V*fc*(k+tau_dl+tau_ctrl))/(Bc*c)))*Sk(k);
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Gives the uplink capacity of Kth drone with GS
        function Sk=get_uplink_capacity(self,freq)
            g_kl=self.get_channelVector(freq);
            p_uk=self.get_uplinkTransmitpower(freq);
            Sk=zeros(self.nb_agents,1);
            
            for k=1:self.nb_agents
                mod_gk=norm(g_kl(k,:));
                ici=0;
                for j=1:self.nb_agents
                    g_kh=(g_kl(k,:).')' ;%g_k is a column vector,g_kh is a row vector
                    if k~=j
                        g_j=g_kl(j,:).';%g_j should is a column vector
                        ici=ici + p_uk(j)*(abs(g_kh*g_j))^2;
                    end
                    
                end
                n0=self.GS.N0;
                bw=self.GS.bw;
                Sk(k)=log2(1+((p_uk(k)*(mod_gk^4))/(ici + n0*bw*mod_gk^2)));
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function sigstrength=get_signalstrength(self,power_gs,antennapos)
           pathloss=get_pathloss(antennapos);
           sigstrength=(ones(self.nb_agents) .* power_gs)-pathloss;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % give the GS antenna array gain for a given Azimuth and Elevation
        function g = get_GS_gain(self,Az,El)
            if self.GS_required
                ga = phased.ArrayGain('SensorArray',self.GS.sCA);
                g = ga(self.GS.freq,[Az;El]);
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % gives the drone position in a ENU frame
        function xyz = get_xyz(self)
            Pos_ned = self.get_pos_ned();
            %NED to ENU conversion
            xyz = Pos_ned .* [1;1;-1];
            t= xyz(1,:);
            xyz(1,:)=xyz(2,:);
            xyz(2,:)=t;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function angle= get_GS_antenna_angles(self)
            % Gives GS antennas orientation with respect to drone frame
            % It returns a 3xKxL array
            % angle(:,K,L) returns a column vector
            % [phi^{+}_{kl};theta^{+}_{kl};psi^{+}_{kl}] of angles made by
            % Lth antenna in GS with Kth drone
            angle=zeros(3,self.nb_agents,length(self.GS.antenna_pos(1,:)));
            Pos=self.get_pos_ned();
            for k =1:self.nb_agents
                for l=1:length(self.GS.antenna_pos(1,:))
                    drone = self.drones(k);
                    dronepos=Pos(:,k);
                    antennapos=self.GS.antenna_pos(:,l);
                    droneorient=[drone.attitude(1),drone.attitude(2),drone.attitude(3)];
                    ang=getangle_in_DF(dronepos,antennapos,droneorient);
                    angle(1,k,l)=ang(1);
                    angle(2,k,l)=ang(2);
                    angle(3,k,l)=ang(3);
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function angle= get_D_antenna_angles(self)
            % Gives Drone orientation with respect to GS antenna frame
            % It returns a 3xKxL array
            % angle(:,K,L) returns a column vector
            % [phi^{x}_{kl};theta^{x}_{kl};psi^{x}_{kl}] of angles made by 
            % Kth drone with Lth antenna in GS
            angle=zeros(3,self.nb_agents,length(self.GS.antenna_pos(1,:))); 
            Pos=self.get_pos_ned();
            for k =1:self.nb_agents
                for l=1:length(self.GS.antenna_pos(1,:))
                    %drone = self.drones(k);
                    dronepos=Pos(:,k);
                    antennapos=self.GS.antenna_pos(:,l);
                    antennaorient=[self.GS.antenna_ang(1,l),self.GS.antenna_ang(2,l),self.GS.antenna_ang(3,l)];
                    ang=getangle_in_GS(dronepos,antennapos,antennaorient);
                    angle(1,k,l)=ang(1);
                    angle(2,k,l)=ang(2);
                    angle(3,k,l)=ang(3);
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function h_kl=get_PLF(self,freq)
            %returns polarization loss factor
            c=self.GS.c;
            h_kl=zeros(self.nb_agents,length(self.GS.antenna_pos(1,:)));
            for k=1:self.nb_agents
                for l=1:length(self.GS.antenna_pos(1,:))
                    h_kl(k,l)=get_polarization_loss(self.GS.antenna_ang(:,l)',self.drones(k).attitude',self.drones(k).pos_ned,self.GS.antenna_pos(:,l),freq);
                    
                end
            end
            function h_kl=get_polarization_loss(antennaorient,droneorient,dronepos,antennapos,f0)
                dlen=1;
                rotmp=eul2rotm(droneorient,'ZYX');
                rotmx=eul2rotm(antennaorient,'ZYX');
                %c=physconst('LightSpeed');
                p=dronepos-antennapos;
                x1=mtimes(rotmx,p);
                dkl=norm(p);
                thetax_kl=acos(x1(3)/dkl);
                psix_kl=acos(x1(2)/dkl);
                x1=mtimes(rotmp,-p);
                Rxp=(rotmx.')*rotmp;
                thetap_kl=acos(x1(3)/dkl);
                psip_kl=acos(x1(2)/dkl);
                constmat=[0,0;0,1;1,0];
                x=p(1);
                y=p(2);
                z=p(3);
                magxy=norm([x,y]);
                magxz=norm([x,z]);
                firstmat=1/dkl*[ -x*z/magxy, -y*z/magxy,    magxy; ...
                                 -x*y/magxz,      magxz, -z*y/magxz ];
                T=firstmat*Rxp*constmat;
                
                function fout=F(theta,f)
                    if mod(theta,pi)==0
                        fout=0;
                    else
                        fout=(cos((pi*dlen/(c/f))*cos(theta))-cos(pi*dlen/(c/f)))/sin(theta);
                    end
                end
                
                E_ltheta=1;
                E_lpsi=1;
                E_ktheta=1;
                E_kpsi=1;
                El=[(E_ltheta*F(thetax_kl,f0));(E_lpsi*F(psix_kl,f0))];
                Ek=[(E_ktheta*F(thetap_kl,f0));(E_kpsi*F(psip_kl,f0))];
                
                h_kl=El'*T*Ek;
                
                
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function pathloss=get_pathloss(self,freq)
            c=self.GS.c;
            lambda=c/freq;
            dist=self.get_antenna_dist();
            pathloss=zeros(self.nb_agents,length(self.GS.antenna_pos(1,:)));
            for k=1:self.nb_agents
                for l=1:length(self.GS.antenna_pos(1,:))
                    pathloss(k,l)=(lambda/(4*pi*dist(k,l)))^2;
                end
            end
            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function dist= get_antenna_dist(self)
           Pos=self.get_pos_ned();
           dist=zeros(self.nb_agents,length(self.GS.antenna_pos(1,:)));
           for k =1:self.nb_agents
               for l=1:length(self.GS.antenna_pos(1,:))
                   dist(k,l)=norm(Pos(:,k)-self.GS.antenna_pos(:,l));
               end
           end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function uplink=get_uplinkTransmitpower(self,freq)
            
            n0=self.GS.N0;
            bw=self.GS.bw;
            pu =self.GS.pu;
            Pu=1;
            ki_kl=abs(self.get_PLF(freq)).^2;
            uplink=zeros(self.nb_agents,1);
            for k=1:self.nb_agents
                pathloss=self.get_pathloss(freq);
                pow=pu*n0*bw/mean(pathloss(k,:).*ki_kl(k,:));
                uplink(k)=min([pow,Pu]);
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function gain=get_HW_Dipole_AntennaGain(self,angle)
            gain=zeros(1,self.nb_agents);
            for i =1:self.nb_agents
                gain(i)=cos(pi*cos(angle)/2)/sin(angle);
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function add_drone(self, drone_type, p_drone, p_battery, p_sim, p_physics, map)
            self.nb_agents = self.nb_agents + 1;
            drone = Drone(drone_type, p_drone, p_battery, p_sim, p_physics, map);
            self.drones = [self.drones; drone];     
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function add_n_drones(self, drone, n)
            for i = 1:n
                self.add_drone(drone);
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function init_rand_pos(self, map_size)

            for i = 1:self.nb_agents
                self.drones(i).init_rand_pos(map_size);
            end

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function set_pos(self, pos)

            for i = 1:self.nb_agents
                self.drones(i).set_pos(pos(:, i));
            end

        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function set_vel(self, vel)

            for i = 1:self.nb_drones
                self.drones(i).set_vel(vel(:, i));
            end

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function Pos_ned = get_pos_ned(self)
            % Return positions of the agent of the swarm in a matrix shape
            % of size 3 x nb_agents
            %
            %        agent_1   agent_2   ...   agent_N
            %   pn
            %   pe
            %   pd

            Pos_ned = zeros(3, self.nb_agents);

            for i = 1:self.nb_agents
                drone = self.drones(i);
                Pos_ned(:, i) = drone.pos_ned;
            end

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function Vel_ned = get_vel_ned(self)
            % Return velocities of the agents of the swarm in a matrix shape
            % of size 3 x nb_agents
            %        agent_1   agent_2   ...   agent_N
            %   vn
            %   ve
            %   vd
            
            Vel_ned = zeros(3, self.nb_agents);

            for i = 1:self.nb_agents
                drone = self.drones(i);

                phi = drone.attitude(1);
                theta = drone.attitude(2);
                psi = drone.attitude(3);
                Rbi = Rb2i(phi, theta, psi);

                Vel_ned(:, i) = Rbi * drone.vel_xyz;
            end

        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function set_state(self, state)
            Pos_ned = state(repmat([true true true false false false], ...
                self.nb_drones,1));
            Pos_ned = reshape(Pos_ned,3,[]);
            Vel_xyz = state(repmat([false false false true true true], ...
                self.nb_drones,1));
            Vel_xyz = reshape(Vel_xyz,3,[]);
            
            self.set_pos(Pos_ned);
            self.set_vel(Vel_xyz);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function state = get_state(self)
            Pos_ned = self.get_pos_ned();
            Vel_ned = self.get_vel_ned();
            state = [Pos_ned; Vel_ned];
            state = state(:);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function Path_len = get_path_len(self)

            Path_len = zeros(1, self.nb_agents);

            for i = 1:self.nb_agents
                drone = self.drones(i);

                Path_len(1, i) = drone.path_len;
            end

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function set_velocity_commands(self, command)

            for i = 1:self.nb_agents
                drone = self.drones(i);
                drone.command(1) = 0;
                drone.command(2:4) = command(2:4);
            end

        end
 
        function set_vel_commands(self, commands)

            for i = 1:self.nb_agents
                drone = self.drones(i);
                drone.command(1) = 0;
                drone.command(2:4) = commands(:, i);
            end

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function Qt = get_Qt(self)
            Qt = zeros(1, self.nb_agents);

            for i = 1:self.nb_agents
                Qt(i) = self.drones(i).Qt;
            end

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function Q = get_Q(self)
            Q = zeros(1, self.nb_agents);

            for i = 1:self.nb_agents
                Q(i) = self.drones(i).Q;
            end

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function colors = get_colors(self)
            colors = zeros(3, self.nb_agents);

            for i = 1:self.nb_agents
                colors(:, i) = self.drones(i).color;
            end

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function update_state(self,wind,time,freq)
            %uplink=zeros(self.nb_agents,1);
            uplink=self.get_uplinkTransmitpower(freq);
            for i = 1:self.nb_agents
                self.drones(i).update_state(wind,time,uplink(i));
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [vel_commands, collisions] = update_command(self, p_swarm, r_coll, dt)

            % Select the swarm algorithm and call the associated update
            if self.algorithm == "vasarhelyi"
                [vel_commands, collisions] = self.compute_vel_vasarhelyi(p_swarm, r_coll, dt);
            elseif self.algorithm == "olfati_saber"
                [vel_commands, collisions] = self.compute_vel_olfati_saber(p_swarm, r_coll, dt);
            end
            if isempty(self.collisions_history)
                self.collisions_history = collisions;
            else
                self.collisions_history = [self.collisions_history; collisions];
            end
            self.set_vel_commands(vel_commands);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function path_planner_swarm(self, path_type, time)
            % Creates an equivalent drone which will receive swarm
            % commands
            self.equivalent_drone = get_barycenter(self);
            self.equivalent_drone.plan_path(path_type, time);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function equivalent_drone = get_barycenter(self)
            pos = zeros(3, 1);
            vel = zeros(3, 1);

            for i = 1:self.nb_agents
                pos = pos + self.drones(i).pos_ned;
                vel = vel + self.drones(i).vel_xyz;
            end

            pos = pos / self.nb_agents;
            vel = vel / self.nb_agents;
            equivalent_drone = Drone(self.drones(1).drone_type, ...
                self.drones(1).p_drone, self.drones(1).p_battery, ...
                self.drones(1).p_sim, self.drones(1).p_physics, ...
                self.drones(1).map);
            equivalent_drone.set_pos(pos);
            equivalent_drone.set_vel(vel);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function path_planner_individual(self, path_type, time)
            % Each agent creates its waypoints independently
            for i = 1:self.nb_agents
                self.drones(i).plan_path(path_type, time);
            end

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function path_manager_individual(self, time)
            % Each agent creates its path independently
            for i = 1:self.nb_agents
                self.drones(i).path_manager_wing(time);
            end

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function path_follower_individual(self, time)
            % Each agent follows its path independently
            for i = 1:self.nb_agents
                self.drones(i).follow_path(time);
            end

        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function pos_ned_history = get_pos_ned_history(self)
            for i = 1:self.nb_agents
                pos_ned_history(:, (3 * (i - 1) + 1) : (3 * (i - 1) + 3)) = self.drones(i).pos_ned_history;
            end

        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function vel_xyz_history = get_vel_xyz_history(self)
            vel_xyz_history = [];
            for i = 1:self.nb_agents
                vel_xyz_history(:, (3 * (i - 1) + 1) : (3 * (i - 1) + 3)) = self.drones(i).vel_xyz_history;
            end

        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % // TODO: add this function to the SwarmViewer
        fig_handle = draw_agents_energy(self, time, period, fig_handle, axes_lim);


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        record_state(self, time, T, period, is_active_cyl, ...
            collisions, map, dirname);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        record_video(self, time, T, period, fig_handle, path);
    end

end
 