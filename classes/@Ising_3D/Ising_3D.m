classdef Ising_3D
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        L
        N
        J double
        h double
        beta double
        taus 
        maxtau
        from_frozen 
        all_signs
        all_nrgs
        all_neighbors
        total_steps
        burnin_steps
        sim_steps
        avging_steps
        run_num

        X
        Y
        Z

        E
        E_data
        M double
        M_data

        F_spns
        spec
        inv_
        modes
        Emodes
        corr
        Ecorr
        avg_modes
        avg_Emodes
        avg_corr
        avg_Ecorr
        num_avging
        Sc
        Cc
        ScE
        CcE
        Sc_data
        Cc_data
        ScE_data
        CcE_data
        wolff_n

        spin_signs_histories
        px
        px1
        py
        pxx1
        px1y
        pxx1y
        eyes
        %Te
        %Te_data
        %avg_Te
        %Te_Ys

        base_dir
        save_dir


    end

    methods
        function obj = Ising_3D(L,beta,run_num,from_frozen)
            %Ising Construct an instance of this class
            %   Detailed explanation goes here
            if pwd == "C:\Users\seanm\OneDrive\Documents\MATLAB"
                obj.base_dir = "P:/Research/Information_Flocking/matlab/bin";
            else
                obj.base_dir = "/scratch/skelty/Information_Flocking/matlab/bin";
            end
            %obj.base_dir = "/scratch/skelty/Information_Flocking/matlab/bin"; 
            obj.save_dir = "3D_simulations";

            obj.L = double(L);
            obj.N = double(L^3);
            %obj.J = double(J);
            %obj.h = h;
            obj.J = double(1);
            obj.h = 0;
            obj.beta = beta;
            obj.taus = [1,2];
            obj.maxtau = max(obj.taus);
            obj.burnin_steps = 1000;
            obj.burnin_steps = 0;
            obj.sim_steps = 2000;
            %obj.Te_Ys = (2:obj.N);
            obj.run_num = run_num;
            obj.avging_steps = 5;
			rng(obj.run_num);
            obj.from_frozen = from_frozen;
            %obj.spins(1:obj.N) = Node_(0); 

			    
            obj = obj.check_directory(obj.save_dir);
            %obj.base_dir = "../../../bin/simulation_runs";
		    %filename = obj.base_dir + "/"+obj.save_dir+"/L"+obj.L+"/J"+obj.J+"/h"+obj.h+"/beta.44_starting_grid.mat";

            obj = obj.initialize_model();
        end

        function obj = initialize_sign(obj,i,j,k)
            %initialize_sign Sets the sign of the node
            %   Detailed explanation goes here
            
            sign = randi([0,1],1,1);
            if sign == 0
                sign = -1;
            end
            obj.all_signs(i,j,k) = double(sign);
        end

        function obj = initialize_model(obj)
            %initialize_model Initializes 2-D square Lattice for Ising
            %Model
            %   Detailed explanation goes here
            obj.all_signs = zeros([obj.L,obj.L,obj.L],'double');
			%base_dir  = "../../../bin";

			%base_dir = "/scratch/skelty/Information_Flocking/matlab/bin";
            neighbor_file = obj.base_dir + "/grid_neighbors_3D/"+obj.L+".mat";

		    %filename = obj.base_dir + "/L"+obj.L+"/J"+obj.J+"/h"+obj.h+"/beta"+obj.beta+"/starting_grid.mat";
			%fprintf("neighbor file: %s\n\n\n",neighbor_file);
            if isfile(neighbor_file)
                file_ = load(neighbor_file);
                obj.all_neighbors = file_.sav_neighb;
                if size(obj.from_frozen,2) == 1
                    for i = 1:obj.L
                       for j = 1:obj.L 
                           for k = 1:obj.L
                               if obj.from_frozen == true
                                   obj.all_signs(i,j,k) = 1;
                               else
                                   obj = obj.initialize_sign(i,j,k);
                               end
                           end
                       end
                    end
                else
                    disp('Generating from Previous beta State');
                    obj.all_signs = obj.from_frozen;
                end


            else
                disp('Neighbors do not exist in file. Generating and saving')
                %obj.all_neighbors = zeros([obj.L,obj.L,obj.L,6,3],'double');
                %size(obj.all_neighbors)

    
                for node = 1:obj.N
                    if obj.from_frozen == true
                        obj.all_signs(node) = 1;
                    else
                        obj = obj.initialize_sign(node);
                    end
                        [i,j,k] = ind2sub(size(obj.all_signs),node);
                        dummy = [sub2ind(size(obj.all_signs),mod(i,obj.L)+1,j,k),...
                               sub2ind(size(obj.all_signs),mod(i-2,obj.L)+1,j,k),...
                               sub2ind(size(obj.all_signs),i,mod(j,obj.L)+1,k),...
                               sub2ind(size(obj.all_signs),i,mod(j-2,obj.L)+1,k),...
                               sub2ind(size(obj.all_signs),i,j,mod(k,obj.L)+1),...
                               sub2ind(size(obj.all_signs),i,j,mod(k-2,obj.L)+1)]';
                           %dummy = reshape([[mod(i,obj.L)+1,j,k],[mod(i-2,obj.L)+1,j,k],...
                           %    [i,mod(j,obj.L)+1,k],[i,mod(j-2,obj.L)+1,k],...
                           %    [i,j,mod(k,obj.L)+1],[i,j,mod(k-2,obj.L)+1]],[3,6])';
                           %size(dummy)
                           obj.all_neighbors = [obj.all_neighbors dummy];
                           %obj.all_neighbors(i,j,k,:,:) = dummy;

                end
                if size(obj.from_frozen,2) > 1
                    disp('Generating from Previous beta State');
                    obj.all_signs = obj.from_frozen;
                end
                   %end                 

                %end
                
                %fileID = fopen(neighbor_file,'w');
                %fwrite(fileID,obj.all_neighbors,'double');
                %fclose(fileID);

                sav_neighb = obj.all_neighbors;
                save(neighbor_file,"sav_neighb");
                    
            end

            obj = obj.reset_values();
        end

		function obj = reset_values(obj)
            obj = initalize_plotting_params(obj);
            obj.total_steps = 0;
            obj = obj.compute_E();
            obj.E_data = zeros(1,obj.sim_steps);
            obj.E_data(1) = obj.E;
            obj.M_data = zeros(1,obj.sim_steps);
            obj.Sc_data = zeros(1,obj.sim_steps);
            obj.Cc_data = zeros(1,obj.sim_steps);
            obj.ScE_data = zeros(1,obj.sim_steps);
            obj.CcE_data = zeros(1,obj.sim_steps);
            %obj.Te_data = zeros(1,obj.sim_steps);
            %obj.Te_data = zeros(size(obj.taus,2),size(obj.Te_Ys,2),obj.sim_steps);
            obj.num_avging = 0;
            if obj.burnin_steps <= obj.total_steps
                obj = obj.compute_M();
                obj.M_data(1) = obj.M;
                obj = obj.get_autocorrelation();
                obj = obj.compute_Sc();
                obj.Sc_data(1) = obj.Sc;
                obj = obj.compute_Cc();
                obj.Cc_data(1) = obj.Cc;
                obj = obj.get_autocorrelation_E();
                obj = obj.compute_ScE();
                obj.ScE_data(1) = obj.ScE;
                obj = obj.compute_CcE();
                obj.CcE_data(1) = obj.CcE;  
            else
                obj.M = 0;
                obj.Sc = 0;
                obj.Cc = 0;
                obj.ScE = 0;
                obj.CcE = 0;

                %obj.avg_modes = 0;
                %obj.avg_Emodes = 0;
            end
            %obj = obj.update_probs_2();
            %obj.Te = zeros(size(obj.taus,2),size(obj.Te_Ys,2));
            obj.avg_modes = zeros(obj.L,obj.L,obj.L,int16((obj.sim_steps-obj.burnin_steps)/5));
            obj.avg_Emodes = zeros(obj.L,obj.L,obj.L,int16((obj.sim_steps-obj.burnin_steps)/5));
            obj.avg_corr = zeros(obj.L,obj.L,obj.L,int16((obj.sim_steps-obj.burnin_steps)/5));
            obj.avg_Ecorr = zeros(obj.L,obj.L,obj.L,int16((obj.sim_steps-obj.burnin_steps)/5));
            %obj.avg_Te = zeros(size(obj.taus,2),size(obj.Te_Ys,2),int16((obj.sim_steps-obj.burnin_steps)/5));
            obj = obj.update_avg_modes();

		end
        function obj = initalize_plotting_params(obj)
            [obj.X,obj.Y] = meshgrid(1:obj.L,1:obj.L);
            obj.Z = zeros([obj.L,obj.L]);
            for i = 1:obj.L
                for j = 1:obj.L
                    obj.Z(i,j) = obj.all_signs(1,i,j);
                end
            end

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = compute_E(obj)
            %all_signs_up = circshift(obj.all_signs,obj.L);
            %all_signs_down = circshift(obj.all_signs,-obj.L);
            %all_signs_left = circshift(obj.all_signs,1);
            %all_signs_right = circshift(obj.all_signs,-1);
            
            %term1 = -obj.J/2*(obj.all_signs.*all_signs_down + obj.all_signs.*all_signs_left...
            %    +obj.all_signs.*all_signs_right + obj.all_signs.*all_signs_up);
            %term2 = -obj.h*obj.all_signs; 

            
            %obj.all_signs = reshape(obj.all_signs,[obj.L,obj.L])';

            term1 = -obj.J/2*(obj.all_signs.*circshift(obj.all_signs,1,1) + obj.all_signs.*circshift(obj.all_signs,-1,1)...
                + obj.all_signs.*circshift(obj.all_signs,1,2) + obj.all_signs.*circshift(obj.all_signs,-1,2)...
                 + obj.all_signs.*circshift(obj.all_signs,-1,3) + obj.all_signs.*circshift(obj.all_signs,-1,3));
            term2 = -obj.h*obj.all_signs;
            obj.all_nrgs = term1 + term2;
            obj.E = sum(obj.all_nrgs,'all');
            %obj.E

        end
        function obj = compute_M(obj)
            obj.M = 1/obj.N*sum(obj.all_signs,'all');
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = get_autocorrelation(obj)
            spinsigns = obj.all_signs - obj.M;
            %spns = reshape(spinsigns,[obj.L,obj.L])';
            obj.F_spns = fftn(spinsigns);
            obj.spec = conj(obj.F_spns).*obj.F_spns;
            obj.corr = ifftn(obj.spec);
			obj.modes = real(obj.spec);
            %obj.modes = reshape(real(obj.spec),[1,obj.N]);
        end
        function obj = compute_Sc(obj)
            %modes_sqr = power(obj.modes,2);
            modes_sqr = obj.modes;
            pk = modes_sqr/sum(modes_sqr,'all');
            obj.Sc = -1*sum(pk.*log2(pk + eps),'all');
        end
        function obj = compute_Cc(obj)
            %modes_sqr = power(obj.modes,2);
            modes_sqr = obj.modes;
            pk = modes_sqr/max(modes_sqr,[],'all');
            obj.Cc = -1*sum(pk.*log2(pk + eps),'all');
        end

        function obj = get_autocorrelation_E(obj)
            spinE = obj.all_nrgs - mean(obj.all_nrgs,'all');
            %spns = reshape(spinE,[obj.L,obj.L])';
            F_spns = fftn(spinE);
            spec = conj(F_spns).*F_spns;
            obj.Ecorr = ifftn(spec);
			obj.Emodes = real(spec);
            %obj.Emodes = reshape(real(spec),[1,obj.N]);
        end
        function obj = compute_ScE(obj)
            %modes_sqr = power(obj.Emodes,2);
            modes_sqr = obj.Emodes;
            pk = modes_sqr/sum(modes_sqr,'all');
            obj.ScE = -1*sum(pk.*log2(pk + eps),'all');
        end
        function obj = compute_CcE(obj)
            %modes_sqr = power(obj.Emodes,2);
            modes_sqr = obj.Emodes;
            pk = modes_sqr/max(modes_sqr,[],'all');
            obj.CcE = -1*sum(pk.*log2(pk + eps),'all');
        end

        function obj = update_avg_modes(obj)
            if obj.burnin_steps <= obj.total_steps
                if mod(obj.total_steps,obj.avging_steps) == 0
                    obj.num_avging = obj.num_avging + 1;
                    %obj.avg_modes(:,obj.num_avging) = obj.modes;
                    %obj.avg_Emodes(:,obj.num_avging) = obj.Emodes;
                    obj.avg_modes(:,:,:,obj.num_avging) = obj.modes;
                    obj.avg_Emodes(:,:,:,obj.num_avging) = obj.Emodes;
                    obj.avg_corr(:,:,:,obj.num_avging) = obj.corr;
                    obj.avg_Ecorr(:,:,:,obj.num_avging) = obj.Ecorr;
                    %obj.avg_Te(:,:,obj.num_avging) = obj.Te;
                end
            end
        end

        function obj = compute_avg_modes(obj)
            %obj.avg_modes = mean(obj.avg_modes(:,1:obj.num_avging),2);
            %obj.avg_Emodes = mean(obj.avg_Emodes(:,1:obj.num_avging),2);
            obj.avg_modes = mean(obj.avg_modes(:,:,:,1:obj.num_avging),4,"omitnan");
            obj.avg_Emodes = mean(obj.avg_Emodes(:,:,:,1:obj.num_avging),4,"omitnan");
            obj.avg_corr = mean(obj.avg_corr(:,:,:,1:obj.num_avging),4,"omitnan");
            obj.avg_Ecorr = mean(obj.avg_Ecorr(:,:,:,1:obj.num_avging),4,"omitnan");
            %obj.avg_Te = mean(obj.avg_Te(:,:,1:obj.num_avging),3);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %{
        function [x,x1,y,xx1,x1y,xx1y] = set_probs_2(obj)
            hists = obj.spin_signs_histories;
            x = zeros(2,size(obj.taus,2),size(obj.Te_Ys,2));
            x1 = zeros(2,size(obj.taus,2),size(obj.Te_Ys,2));
            y = zeros(2,size(obj.taus,2),size(obj.Te_Ys,2));
            xx1 = zeros(2,2,size(obj.taus,2),size(obj.Te_Ys,2));
            x1y = zeros(2,2,size(obj.taus,2),size(obj.Te_Ys,2));
            xx1y = zeros(2,2,2,size(obj.taus,2),size(obj.Te_Ys,2));
            lastrow = size(obj.spin_signs_histories,1);

            for i = (1:size(obj.taus,2))
                for j = (1:size(obj.Te_Ys,2))
                    xx1y(-1*((hists(lastrow,1)/2) - 1.5),-1*((hists(lastrow-1,1)/2) - 1.5), ...
                        -1*((hists(lastrow-obj.taus(i),obj.Te_Ys(j))/2) - 1.5),i,j) = 1;
                end
            end
            x = sum(xx1y,[2,3]);
            x1 = sum(xx1y,[1,3]);
            y = sum(xx1y,[1,2]);
            xx1 = sum(xx1y,3);
            x1y = sum(xx1y,1);
        end

        function obj = update_probs_2(obj)
            if obj.total_steps == 0
                obj.spin_signs_histories = obj.all_signs;
            else
                obj.spin_signs_histories = cat(1,obj.spin_signs_histories,obj.all_signs);
            end

            if obj.total_steps > obj.maxtau
                %x,x1,y,xx1,x1y,xx1y = set_probxx1y();
                [x,x1,y,xx1,x1y,xx1y] = obj.set_probs_2();
                obj.px = obj.px + x;
                obj.px1 = obj.px1 + x1;
                obj.py = obj.py + y;
                obj.pxx1 = obj.pxx1 + xx1;
                obj.px1y = obj.px1y + x1y;
                obj.pxx1y = obj.pxx1y + xx1y;
                obj.spin_signs_histories = [obj.spin_signs_histories(2:end,:)];
            end
            if obj.total_steps == obj.maxtau
                [obj.px,obj.px1,obj.py,obj.pxx1,obj.px1y,obj.pxx1y] = obj.set_probs_2();
            end
        end

        function cond_prob = conditional_prob_2(obj,big)
            norm_ = double(obj.total_steps+1 - obj.maxtau);
            if big == true
                given_x = obj.px1y/norm_;
                joint_x = obj.pxx1y/norm_;
                dum_given_x = zeros(2,2,2,size(obj.taus,2),size(obj.Te_Ys,2));
                dum_given_x(1,:,:,:,:) = given_x;
                dum_given_x(2,:,:,:,:) = given_x;
            else
                given_x = obj.px1/norm_;
                joint_x = obj.pxx1/norm_;
                dum_given_x = zeros(2,2,1,size(obj.taus,2),size(obj.Te_Ys,2));
                dum_given_x(1,:,:,:,:) = given_x;
                dum_given_x(2,:,:,:,:) = given_x;

            end
            cond_prob = joint_x./dum_given_x;
            cond_prob(isnan(cond_prob)) = 0;

        end

        function obj = transfer_entropy_2(obj)
            if obj.total_steps <= obj.maxtau
                obj.Te = zeros(size(obj.taus,2),size(obj.Te_Ys,2));
                
            else
                px_given_x1y = obj.conditional_prob_2(true);
                px_given_x1 = obj.conditional_prob_2(false);
                dum_px_given_x1 = zeros(2,2,2,size(obj.taus,2),size(obj.Te_Ys,2));
                dum_px_given_x1(:,:,1,:,:) = px_given_x1;
                dum_px_given_x1(:,:,2,:,:) = px_given_x1;
                px_given_x1 = dum_px_given_x1;
                arg_ = px_given_x1y./px_given_x1;
                %px_given_x1y
                %px_given_x1
                %arg_(isnan(arg_)) = 0;
                logarg = log2(arg_);
                logarg(isnan(logarg)) = 0;
                logarg(isinf(logarg)) = 0;

                norm_ = double(obj.total_steps+1 - obj.maxtau);
                pxx1y = obj.pxx1y/norm_;
                obj.Te = sum(pxx1y.*logarg,[1,2,3]);
                

            end

        end
        %}
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %{
        function obj = update_probs(obj)
            if obj.total_steps == 0
                obj.spin_signs_histories = obj.all_signs;
            else
                obj.spin_signs_histories = cat(1,obj.spin_signs_histories,obj.all_signs);
            end

            if obj.total_steps > obj.tau
                %x,x1,y,xx1,x1y,xx1y = set_probxx1y();
                [x,x1,y,xx1,x1y,xx1y] = obj.set_probs();
                obj.px = obj.px + x;
                obj.px1 = obj.px1 + x1;
                obj.py = obj.py + y;
                obj.pxx1 = obj.pxx1 + xx1;
                obj.px1y = obj.px1y + x1y;
                obj.pxx1y = obj.pxx1y + xx1y;
                obj.spin_signs_histories = [obj.spin_signs_histories(2:end,:)];
            end
            if obj.total_steps == obj.tau
                [obj.px,obj.px1,obj.py,obj.pxx1,obj.px1y,obj.pxx1y] = obj.set_probs();
            end
        end


        function [x,x1,y,xx1,x1y,xx1y] = set_probs(obj)
            hists = obj.spin_signs_histories;
            x = zeros(2,obj.N);
            x1 = zeros(2,obj.N);
            y = zeros(2,obj.N);
            xx1 = zeros(2,2,obj.N);
            x1y = zeros(2,2,obj.N);
            xx1y = zeros(4,2,obj.N);
            lastrow = size(obj.spin_signs_histories,1);
            for i = (1:obj.N)
                if hists(lastrow,i) == 1
                    x(:,i) = [1,0];
                else
                    x(:,i) = [0,1];
                end
                if hists(lastrow-1,i) == 1
                    x1(:,i) = [1,0];
                else
                    x1(:,i) = [0,1];
                end
                

                if hists(1,i) == 1
                    y(:,i) = [1,0];
                else
                    y(:,i) = [0,1];
                end
                xx1(:,:,i) = kron(x(:,i),x1(:,i)');
                x1y(:,:,i) = kron(y(:,i),x1(:,i)');
                xx1y(:,:,i) = kron(y(:,i),xx1(:,:,i));
                %xx1y(:,:,i) = x(:,i)*x1(:,i)'*y(1,i)';
                
            end

            %probs = {[x,x1,y,xx1,x1y,xx1y]};
        end

        function cond_prob = conditional_prob(obj,big)
            norm_ = double(obj.total_steps+1 - obj.tau);
            if big == true
                given_x = obj.px1y/norm_;
                joint_x = obj.pxx1y/norm_;
                dum_given_x = zeros(4,2,obj.N);
                for i = (1:obj.N)
                    dum_given_x(:,:,i) = kron(given_x(:,:,i),[1;1]);
                end
            else
                given_x = obj.px1/norm_;
                joint_x = obj.pxx1/norm_;
                dum_given_x = zeros(2,2,obj.N);
                for i = (1:obj.N)
                    dum_given_x(:,:,i) = kron(given_x(:,i),[1,1])';
                end

            end
            %size(joint_x)

            %dum_given_x = zeros(4,2,obj.N);
            %    for i = (1:size(given_x,2))
            %        dum_given_x(:,:,i) = kron([1;1],given_x(:,:,i));
            %    end
            cond_prob = joint_x./dum_given_x;
            cond_prob(isnan(cond_prob)) = 0;

        end

        function obj = transfer_entropy(obj)
            x1 = 1;
            y = obj.Te_Ys;
            if obj.total_steps <= obj.tau
                obj.Te = 0;
                
            else
                px_given_x1y = obj.conditional_prob(true);
                px_given_x1 = obj.conditional_prob(false);
                dum_px_given_x1 = zeros(4,2,obj.N);
                for i = (1:size(px_given_x1,2))
                    dum_px_given_x1(:,:,i) = kron([1;1],px_given_x1(:,:,i));
                end
                px_given_x1 = dum_px_given_x1;
                arg_ = px_given_x1y./px_given_x1;
                %arg_(isnan(arg_)) = 0;
                logarg = log2(arg_);
                logarg(isnan(logarg)) = 0;
                logarg(isinf(logarg)) = 0;

                norm_ = double(obj.total_steps+1 - obj.tau);
                pxx1y = obj.pxx1y/norm_;
                %size(pxx1y)
                %logarg
                obj.Te = sum(pxx1y.*logarg,'all');
                %obj.Te = -1*sum(pxx1y.*logarg);
                

            end

        end
        %}
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%function del_E = get_del_E(obj,i,j,k)
        function del_E = get_del_E(obj,id_)
            flpd_spin = -1*obj.all_signs(id_);
            neighbor_signs = obj.all_signs(obj.all_neighbors(1:6,id_));
            del_E = double(-2*flpd_spin*(obj.J*sum(neighbor_signs) + obj.h));
        end
        %    flpd_spin = -1*obj.all_signs(i,j,k);
        %    neighbor_signs = (1:6);
        %    for neighb = 1:6
        %        doop = obj.all_neighbors(i,j,k,neighb,:);
        %        neighbor_signs(neighb) = obj.all_signs(doop(1),doop(2),doop(3));
        %    end
            %neighbor_signs
        %    del_E = double(-2*flpd_spin*(obj.J*sum(neighbor_signs) + obj.h));
        %end

        function obj = update_node(obj,id_)
            flpd_spin = -1*obj.all_signs(id_);

            %x = mod(id_,obj.L);
            %    if x == 0
            %        x = obj.L;
            %    end
            del_E = obj.get_del_E(id_);
            
            if del_E < 0
                %obj.spins(id_).sign = flpd_spin;
                obj.all_signs(id_) = flpd_spin;
                %obj.Z(x,ceil(id_/obj.L)) = flpd_spin;
                obj.E = obj.E + del_E;
                obj.all_nrgs(id_) = -1*obj.all_nrgs(id_);
                for neighb = 1:6
                    nn = obj.all_neighbors(neighb,id_);
                    %obj.all_nrgs(nn) = obj.all_nrgs(nn) - obj.J*flpd_spin;
                    dummy = -obj.J/2*obj.all_signs(nn)*sum(obj.all_signs(obj.all_neighbors(1:6,nn)),'all');
                    obj.all_nrgs(nn) = dummy;
                end
            else
                p = exp(-1*obj.beta*del_E);
                u = rand;
                if p > u
                    %obj.spins(id_).sign = flpd_spin;
                    obj.all_signs(id_) = flpd_spin;
                    %obj.Z(x,ceil(id_/obj.L)) = flpd_spin;
                    obj.E = obj.E + del_E;
                    obj.all_nrgs(id_) = -1*obj.all_nrgs(id_);
                    for neighb = 1:6
                        nn = obj.all_neighbors(neighb,id_);
                        %obj.all_nrgs(nn) = obj.all_nrgs(nn) - obj.J*flpd_spin;
                        dummy = -obj.J/2*obj.all_signs(nn)*sum(obj.all_signs(obj.all_neighbors(1:6,nn)),'all');
                        obj.all_nrgs(nn) = dummy;
                    end

                end
            end
        end


        function obj = update_properties(obj)
            obj.total_steps = obj.total_steps+1;
            obj.E_data(1+obj.total_steps) = obj.E;
            obj = obj.initalize_plotting_params();
            %obj = obj.update_probs_2();
            if obj.burnin_steps <= obj.total_steps
                obj = obj.compute_M();
                obj.M_data(1+obj.total_steps) = obj.M;
                obj = obj.get_autocorrelation();
                obj = obj.compute_Sc();
                obj.Sc_data(1+obj.total_steps) = obj.Sc;
                obj = obj.compute_Cc();
                obj.Cc_data(1+obj.total_steps) = obj.Cc;
                obj = obj.get_autocorrelation_E();
                obj = obj.compute_ScE();
                obj.ScE_data(1+obj.total_steps) = obj.ScE;
                obj = obj.compute_CcE();
                obj.CcE_data(1+obj.total_steps) = obj.CcE;
                %obj = obj.transfer_entropy_2();
                %obj.Te_data(:,:,1+obj.total_steps) = obj.Te;
                obj = obj.update_avg_modes();
            end
        end

        function obj = update_all_nodes(obj)
            %close;
            nodes = (1:obj.N);
            nodes = nodes(randperm(length(nodes)));
            for i = 1:obj.N
                obj = obj.update_node(nodes(i));
            end
            obj = obj.update_properties();
        end



        function obj = update_all_and_plot(obj,num_steps)
            close;
            colormap(gray);
            for i = obj.total_steps+1:(num_steps + obj.total_steps)
                obj = obj.update_all_nodes();
                obj = obj.plot_2(i);
                %pause(.1)
            end
        end

        function obj = update_all_and_plot_hybrid(obj,num_steps1,num_steps2)
            close;
            colormap(gray);
            for i = obj.total_steps+1:(num_steps1 + obj.total_steps)
                obj = obj.update_all_nodes();
                obj = obj.plot_2(i);
                %pause(.1)
            end

            for i = obj.total_steps+1:(num_steps2 + obj.total_steps)
                obj = obj.wolffAlgorithm();
                obj = obj.plot_2(i);
                %pause(.1)
            end
        end

        function obj = plot_2(obj,i)
            tiledlayout(4,4)
            nexttile([1,2])
            
            plot(linspace(0,i,i+1),obj.E_data(1:obj.total_steps+1))
            title('Energy')
            
            nexttile([1,2])
            plot(linspace(0,i,i+1),obj.M_data(1:obj.total_steps+1))
            title('Magnetization')
            text( -.5, 1, '3D Simulation', 'FontSize', 14', 'FontWeight', 'Bold', ...
            'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' ) ;
            nexttile([2,2])
            imshow(obj.Z, 'InitialMagnification', 1000);
            nexttile([2,2])
            plot(linspace(0,i,i+1),obj.ScE_data(1:obj.total_steps+1))
            title('Configurational Entropy (E)')
            %plot(linspace(0,i,i+1),obj.Te_data(1:obj.total_steps+1))
            %title('Transfer Entropy')
            nexttile([1,2])
            plot(linspace(0,i,i+1),obj.Sc_data(1:obj.total_steps+1))
            title('Configurational Entropy')
            nexttile([1,2])
            plot(linspace(0,i,i+1),obj.Cc_data(1:obj.total_steps+1))
            title('Configurational Complexity')
            
            refreshdata
            drawnow
        end

        function obj = plot(obj,i)
            tiledlayout(4,4)
            nexttile([1,2])
            plot(linspace(0,i,i+1),obj.E_data(1:obj.total_steps+1))
            title('Energy')
            nexttile([1,2])
            plot(linspace(0,i,i+1),obj.M_data(1:obj.total_steps+1))
            title('Magnetization')
            nexttile([2,2])
            plot(linspace(0,i,i+1),obj.ScE_data(1:obj.total_steps+1))
            title('Configurational Entropy (E)')
            %imshow(obj.Z, 'InitialMagnification', 1000);
            nexttile([2,2])
            plot(linspace(0,i,i+1),obj.CcE_data(1:obj.total_steps+1))
            title('Configurational Complexity (E)')
            %plot(linspace(0,i,i+1),obj.Te_data(1:obj.total_steps+1))
            %title('Transfer Entropy')
            nexttile([1,2])
            plot(linspace(0,i,i+1),obj.Sc_data(1:obj.total_steps+1))
            title('Configurational Entropy')
            nexttile([1,2])
            plot(linspace(0,i,i+1),obj.Cc_data(1:obj.total_steps+1))
            title('Configurational Complexity')
             text( 0.5, 0, '3D Simulation', 'FontSize', 14', 'FontWeight', 'Bold', ...
            'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' ) ;
            refreshdata
            drawnow
        end

        function obj = update_all_no_plot(obj,num_steps)
            for i = 1:num_steps
                obj = obj.update_all_nodes();
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function obj = generate_starting_grid(obj,file_ext)
            T_c = .222;
            x1 = 0:45;
            deltax = .002;

            num_steps = 2000;
            if strcmp(file_ext,'wolff_simulation') | strcmp(file_ext,"wolff_simulation_3D")
                num_steps = 500;
                disp('Wolff Sim')
            end
			%y1 = T_c - x1.^2*((T_c-.3)/(30^2));
			%y2 = T_c + x1(2:31).^2*((.6-T_c)/(30^2));
            %betas = cat(2,flip(y1,2),y2);
            %ind_ = find(betas == obj.beta)-1;
            y1 = T_c - x1.*deltax;
			y2 = T_c + x1(2:46).*deltax;
            betas = cat(2,flip(y1,2),y2);
            ind_ = find(betas == obj.beta)-1;
            if ~isempty(ind_) & (ind_ > 0)
                beta = betas(ind_);
                filename = obj.base_dir + "/"+file_ext+"/L"+obj.L+"/J"+obj.J+"/h"+obj.h+"/beta"+beta+"/run"+obj.run_num+".mat";

			    %obj = obj.check_directory_starting_grids();
                if exist(filename)
                    %disp('Starting Grid Exists')
                    obj.all_signs = load(filename).allsigns;
                    obj = obj.reset_values();
                else
                    filename = obj.base_dir + "/"+file_ext+"/L"+obj.L+"/J"+obj.J+"/h"+obj.h+"/beta"+obj.beta+"/run"+obj.run_num+".mat";
                    disp('Generating Starting Grid');
                    obj = obj.check_directory(file_ext);
                    obj = obj.initialize_model();
                    %%{
                    for i = 1:num_steps
                        obj = obj.update_all_nodes();
                    end
                    %%}
                    %obj = obj.update_all_and_plot(1000);
			        allsigns = obj.all_signs;
                    save(filename,"allsigns");
                    obj = obj.reset_values();
                end
            else
                filename = obj.base_dir + "/"+file_ext+"/L"+obj.L+"/J"+obj.J+"/h"+obj.h+"/beta"+obj.beta+"/run"+obj.run_num+".mat";
                disp('Generating Starting Grid');
                obj = obj.check_directory(file_ext);
                obj = obj.initialize_model();
                %%{
                for i = 1:num_steps
                    obj = obj.update_all_nodes();
                end
                %%}
                %obj = obj.update_all_and_plot(1000);
			    allsigns = obj.all_signs;
                save(filename,"allsigns");
                obj = obj.reset_values();
            end

        end


        function obj = runsim_MH(obj)
            tstart = tic;
            tic
            obj.avging_steps = 5;
            obj.burnin_steps = 500;
            %obj.burnin_steps = obj.burnin_steps;;
            num_steps = obj.sim_steps;
            obj.check_directory(obj.save_dir);
            obj = obj.reset_values();
			filename = obj.base_dir + "/"+obj.save_dir+"/L"+obj.L+"/J"+obj.J+"/h"+obj.h+"/beta"+obj.beta+"/run"+obj.run_num+".mat";

            if ~exist(filename,'file')
                obj = obj.update_all_and_plot(obj.burnin_steps);
                obj.sim_steps = num_steps;
                obj = obj.reset_values();
                obj.burnin_steps = 0;
                obj = obj.update_all_no_plot(obj.sim_steps);

                tEnd = toc(tstart);
                toc
                obj = obj.save_props(filename,tEnd);
            else
                dummy = load(filename);
                if ~isfield(dummy,"L")
                    obj = obj.update_all_no_plot(obj.burnin_steps);
                    obj = obj.reset_values(); 
                    obj.burnin_steps = 0;
                    obj = obj.update_all_no_plot(obj.sim_steps);
                    tEnd = toc(tstart);
                    toc
                    obj = obj.save_props(filename,tEnd);
                end    
            end
            
        end

        function obj = runsim_Wolff(obj)
            tstart = tic;
            tic
            file_ext = obj.save_dir;
            obj.avging_steps = 2;
            num_steps = 500;
            %obj.avging_steps = 8;
            %num_steps = 2000;
            obj.sim_steps = num_steps;
            obj.burnin_steps = 500;
            obj.check_directory(file_ext);
			filename = obj.base_dir + "/"+file_ext+"/L"+obj.L+"/J"+obj.J+"/h"+obj.h+"/beta"+obj.beta+"/run"+obj.run_num+".mat";

            obj = obj.reset_values();
            if ~exist(filename,'file')

                %obj = obj.generate_starting_grid(file_ext);
                %obj = obj.wolffAlgorithm();
                
                obj = obj.update_all_no_plot(obj.burnin_steps);
                obj.sim_steps = num_steps;
                obj = obj.reset_values();
                obj.burnin_steps = 0;
                
                for i = (1:num_steps)
                    obj = obj.wolffAlgorithm();
                    %obj = obj.autocorrelation_t(true);
                end
                tEnd = toc(tstart);
                toc
                obj = obj.save_props(filename,tEnd);
            else
                dummy = load(filename);
                if ~isfield(dummy,"L")
                    obj = obj.reset_values();
                    obj.burnin_steps = 0;
                    for i = (1:num_steps)
                        obj = obj.wolffAlgorithm();
                        %obj = obj.autocorrelation_t(true);
                    end
                    tEnd = toc(tstart);
                    toc
                    obj = obj.save_props(filename,tEnd);
                end
            end
        end
        function obj = save_props(obj,filename,end_time)
                runnum = obj.run_num;
                beta = obj.beta;
                L = obj.L;
                %Te_data = reshape(obj.Te_data,size(obj.taus,2)*size(obj.Te_Ys,2),num_steps+1);
				%Te_data = Te_data(:,obj.burnin_steps+1:num_steps+1);
                run_data = [obj.E_data(obj.burnin_steps+1:obj.sim_steps);
					obj.M_data(obj.burnin_steps+1:obj.sim_steps);
                    obj.Sc_data(obj.burnin_steps+1:obj.sim_steps);obj.Cc_data(obj.burnin_steps+1:obj.sim_steps);
                    obj.ScE_data(obj.burnin_steps+1:obj.sim_steps);obj.CcE_data(obj.burnin_steps+1:obj.sim_steps);
                    obj.burnin_steps+1:obj.sim_steps];
				
				taus = obj.taus;
				%Te_Ys = obj.Te_Ys;

                obj = obj.compute_avg_modes();
                modes = obj.avg_modes;
                emodes = obj.avg_Emodes;
                corr = obj.avg_corr;
                ecorr = obj.avg_Ecorr;
                %avg_Te = obj.avg_Te;
				allsigns = obj.all_signs;

				%save(filename,"L_","beta_","runnum","run_data","taus","modes","emodes"...
                %    ,"corr","ecorr","avg_Te","allsigns");
                save(filename,"L","beta","runnum","end_time","run_data","taus","modes","emodes"...
                    ,"corr","ecorr","allsigns");
        end


        function obj = wolffAlgorithm(obj)
            %WOLFFALGORITHM Generate a 3d Ising state using the Wolff Algorithm
            %   INPUT
            %       seed   {LxLxL grid of +1 and -1s}   Starting grid OR
            %              {L, positive integer}      Size of grid
            %       beta   {positive real}      Inverse temperature
            %
            %   Generates a single state of the 3d Ising model on the square lattice
            %   with periodic boundary conditions and ferromagnetic coupling. If the
            %   input SEED is a positive integer then it is the side length of the grid
            %   to be created. Otherwise it is the initial state the algorithm will
            %   work on. In the former case this function will use a lot more timesteps
            %   to ensure that the first state that pops out is close to equilibrium.
            %   Otherwise fewer flips are used to generate a new state.
            %
            acceptP = 1 - exp(-2*obj.beta); 


            
            %tic
            totalFlipped = 0;
            while totalFlipped < obj.N
                stack = randi(obj.N);
                cluster = stack;
                %checked = stack;
                while not(isempty(stack))
                    idx = stack(1);
                    nbd = obj.all_neighbors(:,idx);
                    %nbd = mod([idx+1 idx-1 idx+L idx-L]-1,N)+1;
                    for nn = 1:6
                        if obj.all_signs(idx) == obj.all_signs(nbd(nn))
                        %if state(idx)==state(nbd(nn))
                            if all(nbd(nn)~=cluster)
                            if rand() < acceptP
                                stack = [stack nbd(nn)];
                                cluster = [cluster nbd(nn)];
                            end
                            end
                        end
                    end
                    stack(1) = [];
                end
                totalFlipped = totalFlipped + length(cluster);
                obj.wolff_n = [obj.wolff_n length(cluster)];
                obj.all_signs(cluster) = -obj.all_signs(cluster);
                %obj.Z(cluster) = obj.all_signs(cluster);
                %obj = obj.plot(tt);
            end
            obj = obj.compute_E();
            obj = obj.update_properties();
            %toc
        end



        function obj = check_directory(obj,file_ext)
            %base_dir = "P:\Research\Information_Flocking\matlab\bin\simulation_runs\";
			%base_dir = "/scratch/skelty/Information_Flocking/matlab/bin/simulation_runs/";
            %base_dir = "../../../bin/simulation_runs";
            basedir = obj.base_dir+"/"+file_ext;
			dir1 = basedir+"/L"+string(obj.L);
			dir2 = basedir+"/L"+string(obj.L)+"/J"+string(obj.J);
			dir3 = basedir+"/L"+string(obj.L)+"/J"+string(obj.J)+"/h"+string(obj.h);
			dir4 = basedir+"/L"+string(obj.L)+"/J"+string(obj.J)+"/h"+string(obj.h)+"/beta"+string(obj.beta);
            if ~exist(basedir, 'dir')
               mkdir(basedir)
            end
            if ~exist(dir1, 'dir')
               mkdir(dir1)
            end
            if ~exist(dir2, 'dir')
               mkdir(dir2)
            end
            if ~exist(dir3, 'dir')
               mkdir(dir3)
            end
            if ~exist(dir4, 'dir')
               mkdir(dir4)
            end
        end
		function obj = check_directory_starting_grids(obj)
            %base_dir = "P:\Research\Information_Flocking\matlab\bin\simulation_runs\";
			%base_dir = "/scratch/skelty/Information_Flocking/matlab/bin/simulation_runs/";
            %base_dir = "../../../bin/simulation_runs";
            basedir = obj.base_dir+"/starting_grids";
			dir1 = basedir+"/L"+string(obj.L);
			dir2 = basedir+"/L"+string(obj.L)+"/J"+string(obj.J);
			dir3 = basedir+"/L"+string(obj.L)+"/J"+string(obj.J)+"/h"+string(obj.h);
			dir4 = basedir+"/L"+string(obj.L)+"/J"+string(obj.J)+"/h"+string(obj.h)+"/beta"+string(obj.beta);
            if ~exist(basedir, 'dir')
               mkdir(basedir)
            end
            if ~exist(dir1, 'dir')
               mkdir(dir1)
            end
            if ~exist(dir2, 'dir')
               mkdir(dir2)
            end
            if ~exist(dir3, 'dir')
               mkdir(dir3)
            end
            if ~exist(dir4, 'dir')
               mkdir(dir4)
            end
        end



    end
end
