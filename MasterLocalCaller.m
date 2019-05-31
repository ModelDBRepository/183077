%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  %%% Master Caller %%%					  	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MasterLocalCaller()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CURRENT MAIN METHOD %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clear; clc; %clear variable space and screen
rehash; %refresh all data files
over = tic;

% % Run Test Function % %
test([]);
% % % % % % % % % % % % %

t = toc(over);
fprintf('Total Time Used: %.1fs\n',t);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Testing Functions %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%% BEGIN SINGLETON TEST %%%%%%%%%%%%%%%%%%%
    function test(~)
        % Create neuron and set its properties
        disp('Creating Neuron')
        init = tic;
        neuron = TorusNeuronMod();
        ss = 0.025; % Step size in ms (McGillivray 2012)
        
        %%% Neuron target parameters %%%
        neuron.sigmaB = 0.42 ;
        neuron.Ws = 0.7915 ;
        neuron.I_bias = -9.39 ;
        neuron.N_xi = 1.0 ; 
        neuron.f_h = 0.2998 ; 
        neuron.f_ca = 0.737 ;  

        % Other commonly altered Parameters
        tag = '070715-test'; % Name tag to attach to file output
        saveStruct = true;    % boolean to save anything
        saveAllTraces = false; % boolean to save currents/conds 
            % (depends on retI/G) and Vm or not
            
        % Set and initialize other properties
        threshold = -25;         % For spike time detection
        q = 100;                 % For VPD (Vonderschen, 2011)
        numTrials = 20;          % Number of simulations per stimulus type
        nChirps = 5;             % Number of chirps   
        alpha = 0.01;            % Alpha of the SI score
        siglev = 0;              % The needed PSTH level for CSI significance 
        numProcs = 1;            % Number of processors to use
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% File Input Paths %%%%%
        % Read in and format input ELL data  
        IN_DATA = 'Data/InputData_unnorm_07-24-14.mat';
        IN_SYN = 'Data/InputSynapses_unnorm_07-24-14.mat';
        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        

        % Make Pool
        poolobject = gcp('nocreate');
        if isempty(poolobject) && ~(numProcs==1)
            parpool(numProcs);
        end

        % Set the random seed for the RNG 
        rng(1);   
        
        % Check for output file existence
        if saveStruct
            if exist( sprintf('sRes_%s',tag) ,'file')
                disp('Output file already exists!');
                return;
            end
        end
        
        % Read in  input data structures
            inputdata = load(IN_DATA);
            a = inputdata.a;
            c = inputdata.c;
        
        % Build synapses from ELL data (1-4 = small chirps, 5 = big chirps)
        if ~exist(IN_SYN,'file')
            for i = 1:length(a)
                fprintf('\tOn Synapse %d\n',i);
                disp(length(a(i).avgEpsth));
                s = sprintf('Chirp%d-E_avg',i); 
                synapseStruct(i).Esyn = AlphaSynapse(s,a(i).avgEpsth,a(i).chirp,'E'); %#ok<AGROW>
                fprintf('\t\tFinished E-Synapse %d\n',i);            
                s = sprintf('Chirp%d-I_avg',i); 
                synapseStruct(i).Isyn = AlphaSynapse(s,a(i).avgIpsth,a(i).chirp,'I'); %#ok<AGROW>
                fprintf('\t\tFinished I-Synapse %d\n',i);            
            end          
            save(IN_SYN,'synapseStruct');
        else
            insyn = load(IN_SYN);
            synapseStruct = insyn.synapseStruct;
        end
        
        % Preparations and preallocations for parallel processing
        for i = length(a):-1:1 % c = chirps, a = chirps+noise. Need a for plotting.
            % Make the neuron clone
            neurons(i) = neuron.deepCopy();
            % Add synapses (i.e. specific for a given chirp)
            neurons(i).addAlphaSynapse(synapseStruct(i).Esyn);
            neurons(i).addAlphaSynapse(synapseStruct(i).Isyn);
            psths{i} = zeros(round(length(neurons(i).getAlphaSynapse(1).conved)/10)/ss,1);
            for j = 1:numTrials
                stsCa{i,j} = [];  
                gcells{i,j} = [];
                Vms{i,j} = [];
            end
        end
        
        elapsed = toc(init);
        fprintf('Initialization Time Used: %.1fs\n',elapsed);
        
        % Run simulation (multiThreadable)
        for i = 1:length(c) % c = chirps, a = chirps+noise
            chirpt = tic;
            neuronCopy = neurons(i);
            % Prepare necessary variables
            tmax = round(length(neuronCopy.getAlphaSynapse(1).conved)/10);
            % Run the multiple trials for the given input
            if numProcs > 1
                parfor j = 1:numTrials
                    % Run simulation
                    [~,v,gs] = NSUtils.FastNeuronEulerMaruyamaRetGs(neuronCopy,ss,tmax); 
                    % Procure psth and binaries
                    [psth,~,sts] = NSUtils.getPsthStsAndBinaries(v,threshold,ss);
                    stsCa{i,j} = sts;  %#ok<NASGU,PFOUS>
                    psthsTemp{i,j} = psth;
                    % Save currents, if needed
                    if saveAllTraces
                        gcells{i,j} = gs;
                        Vms{i,j} = v;
                    end
                end
            else
                for j = 1:numTrials
                    % Run simulation
                    [~,v,gs] = NSUtils.FastNeuronEulerMaruyamaRetGs(neuronCopy,ss,tmax); 
                    % Procure psth and binaries
                    [psth,~,sts] = NSUtils.getPsthStsAndBinaries(v,threshold,ss);
                    stsCa{i,j} = sts;  %#ok<NASGU,PFOUS>
                    psthsTemp{i,j} = psth;
                    % Save currents, if needed
                    if saveAllTraces
                        gcells{i,j} = gs;
                        Vms{i,j} = v;
                    end
                end                
            end
            % Average psths
            for k = 1:numTrials; psths{i} = psthsTemp{i,k} + psths{i}; end
            psths{i} = psths{i}/numTrials;
            tt = toc(chirpt);
            fprintf('\tTime Used for Chirp %d: %.1fs\n',i,tt);
        end
        
        postsim = tic;
        
        % Compute Simulation stats (CSI and VPD)
        stats = NSUtils.computePostSimulationStats...
            (psths,nChirps,ss,stsCa,q,numTrials,alpha,siglev);
        
        % Save results
        disp('Saving Results');
        if saveAllTraces && saveStruct
            save(sprintf('sRes_%s',tag),'neuron','stats','stsCa','psths','gcells','Vms');
        elseif saveStruct
            save(sprintf('sRes_%s',tag),'neuron','stats','stsCa','psths'); %#ok<UNRCH>
        end
        
        % Plot the results
        disp('Plotting Results');
        graphNoise = false;
        NSUtils.plotMultiTrialOutput(neuron,ss,synapseStruct,stsCa,psths,graphNoise);
        disp(stats);
        outt = toc(postsim);
        fprintf('Post Simulation Time Used: %.1fs\n',outt);
        
    end    
%%%%%%%%%%%%%%%%%%% END SINGLETON TEST %%%%%%%%%%%%%%%%%%%

end

