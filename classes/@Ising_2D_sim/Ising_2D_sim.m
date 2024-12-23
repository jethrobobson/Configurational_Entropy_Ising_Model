classdef Ising_2D_sim
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        Ls
        N
        J double
        h double
        betas
        taus
        Te_Ys
        from_frozen 
        all_signs
        all_neighbors
        total_steps

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
        Sc
        Cc
        Sc_data
        Cc_data

        spin_signs_histories
        px
        px1
        py
        pxx1
        px1y
        pxx1y
        eyes
        Te
        Te_data

        sim_steps
        burnin_steps
        wolff_cutoff
        base_dir
        save_dir
        run_num

    end

    methods
        function obj = Ising_2D_sim(Ls,run_num)
            obj.Ls = double(Ls);
            obj.J = 1;
            obj.h = 0;
            obj.taus = linspace(1,8,8);
            obj.run_num = run_num;
            %obj.base_dir = "P:\Research\Information_Flocking\matlab\bin";
            if pwd == "C:\Users\seanm\OneDrive\Documents\MATLAB"
                obj.base_dir = "P:/Research/Information_Flocking/matlab/bin";
            else
                obj.base_dir = "/scratch/skelty/Information_Flocking/matlab/bin";
            end 
            obj.save_dir = "2D_simulations";
			

            %{
            T_c = .44072;
            x1 = 0:45;
            deltax = .002;
            y1 = T_c - x1.*deltax;
			y2 = T_c + x1(2:46).*deltax;
            obj.betas = flip(cat(2,flip(y1,2),y2),2);
            %obj.betas = cat(2,y0(1:29),obj.betas);
            
            %}


            n = 30;
            T_c = .44072;

            x1 = 0:n;
            deltax = .001;
            y1 = T_c - x1.*deltax;
			y2 = T_c + x1(2:n+1).*deltax;
			betas_0 = flip(cat(2,flip(y1,2),y2),2);
            y1 = T_c - .0005 - x1.*deltax;
			y2 = T_c - .0005 + x1(2:n+1).*deltax;
            betas = flip(cat(2,flip(y1,2),y2),2);

            n2 = 20;
            x1 = 0:n2;
            deltax = .003;
            %y1 = betas(1) + x1(2:n2).*deltax;
			y1 = betas_0(1) + x1(2:n2).*deltax;

            n2 = 80;
            x1 = 0:n2;
            deltax = .001;
            y2 = betas_0(size(betas_0,2)) - x1(2:n2+1).*deltax;
			%y2 = betas(size(betas,2)) - x1(2:n2+1).*deltax;
            betas = cat(2,flip(y1,2),betas);
            obj.betas = cat(2,betas,y2);
            %obj.wolff_cutoff = [obj.betas(5),.43];
			obj.wolff_cutoff = [obj.betas(20),.43];

            %cold_reg = 1./linspace(.25,1.75,14);
            %obj.betas = cat(2,cold_reg,obj.betas);

            %hot_reg = 1./linspace(2.9,4,12);
            %obj.betas = cat(2,obj.betas,hot_reg);

            obj.save_dir = "2D_simulations";

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function obj = runsim_Wolff(obj)
            tic;
            k = 0;
			fprintf("L = %d\n",obj.L)
			fprintf("Wolff Algorithm")
			fprintf("Running Simulation: Run %d\n",obj.run_num);
            for beta = obj.betas

                mod_ = Ising_2D(obj.L,beta,obj.taus,obj.run_num,true);
                mod_ = mod_.runsim_Wolff();
                fprintf("%d out of %d: %.5f ",k+1,size(obj.betas,2),beta);
                if mod(k,10) == 0
                    fprintf('X ');
                else
                    fprintf('* ');
                end
                k = k+1;
            end
            fprintf('\n');
            toc;
            fprintf('End Simulation');
        end

        function obj = runsim_MH(obj,run_num)
            tic;
            k = 0;
			fprintf("L = %d\n",obj.L)
			fprintf("Glauber Dynamics")
			fprintf("Running Simulation: Run %d\n",run_num);
            for beta = obj.betas

                mod_ = Ising_2D(obj.L,beta,obj.taus,run_num,true);
                mod_ = mod_.runsim_MH();
                fprintf("%d out of %d: beta = %.5f ",k+1,size(obj.betas,2),beta);
                if mod(k,10) == 0
                    fprintf('X ');
                else
                    fprintf('* ');
                end
                k = k+1;
            end
            fprintf('\n');
            toc;
            fprintf('End Simulation');
        end

        function obj = runsim_hybrid(obj)
            
            k = 0;

            for L = obj.Ls
                k = 0;
                fprintf("L = %d\n",L)
			    fprintf("Running Simulation: Run %d\n\n",obj.run_num);
                startL = tic;
                for beta = obj.betas
                    fprintf('\n');
                    starting_grid = generate_starting_grid(obj,L,beta);
                    %disp(size(starting_grid));
                    mod_ = Ising_2D(L,beta,obj.run_num,starting_grid);
                    %if (beta < obj.wolff_cutoff(1)) && (beta > obj.wolff_cutoff(2))
                    if (beta < obj.wolff_cutoff(1))
                        fprintf(" :Wolff: ")
                        mod_ = mod_.runsim_Wolff();
                    else
                        fprintf(" :Metropolis-Hastings: ")
                        mod_ = mod_.runsim_MH();
                    end
    
                    fprintf("%d out of %d: beta = %.5f ",k+1,size(obj.betas,2),beta);
                    if mod(k,10) == 0
                        fprintf('X ');
                    else
                        fprintf('* ');
                    end
                    k = k+1;
                end
            fprintf('\n\n Total Elapsed Time for L=%d:',L);

            toc(startL);
            fprintf('End Simulation\n\n');
            end
        end



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function allsigns = generate_starting_grid(obj,L,beta)
            file_ext = obj.save_dir;
            num_steps = 2000;

            ind_ = find(obj.betas == beta)-1;
            if ~isempty(ind_) & (ind_ > 0)
                beta = obj.betas(ind_);
                filename = obj.base_dir + "/"+file_ext+"/L"+L+"/J"+obj.J+"/h"+obj.h+"/beta"+beta+"/run"+obj.run_num+".mat";

			    %obj = obj.check_directory_starting_grids();
                
                if exist(filename)
                    %disp('Starting Grid Exists')
                    allsigns = load(filename).allsigns;
                    %obj = obj.reset_values();
                else
                    filename = obj.base_dir + "/"+file_ext+"/L"+L+"/J"+obj.J+"/h"+obj.h+"/beta"+beta+"/run"+obj.run_num+".mat";
                    disp('Generating Starting Grid');
                    mod_ = Ising_2D(L,beta,obj.run_num,true);
                    mod_ = mod_.check_directory(file_ext);
                    mod_ = mod_.initialize_model();
                    %%{
                    for i = 1:num_steps
                        mod_ = mod_.update_all_nodes();
                    end
                    %%}
                    %obj = obj.update_all_and_plot(1000);
			        allsigns = mod_.all_signs;
                    save(filename,"allsigns");
                    %mod_ = mod_.reset_values();
                end
            else
                filename = obj.base_dir + "/"+file_ext+"/L"+L+"/J"+obj.J+"/h"+obj.h+"/beta"+beta+"/run"+obj.run_num+".mat";
                disp('Generating Starting Grid (From Simulation)');
                tic;
                mod_ = Ising_2D(L,beta,obj.run_num,true);
                mod_ = mod_.check_directory(file_ext);
                mod_ = mod_.initialize_model();
                %%{
                for i = 1:num_steps
                    mod_ = mod_.update_all_nodes();
                end
                %%}
                %obj = obj.update_all_and_plot(1000);
		        allsigns = mod_.all_signs;
                save(filename,"allsigns");
                toc
                %mod_ = mod_.reset_values();
            end

        end

    end
end
