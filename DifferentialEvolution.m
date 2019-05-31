%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                           % MASTER FUNCTION %                             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DifferentialEvolution()

                    %%%%%
    %%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%
                    %%%%%
                    
    % Input and Output File Parameters %
    IN_SYN = 'Data/InputSynapses_unnorm_07-24-14.mat';
    TS_neuron_SC_file = 'Data/Cell43/cell43_062514_TS_SC.mat';
    TS_neuron_BC_file = 'Data/Cell43/cell43_062514_TS_BC.mat';
    tag = '012815-fh0ft0'; % String to add to end of file
    
    % Algorithm Parameters %
    F = 0.5;        % The differential weight (or factor of differentiation)
    CR = 0.9;       % The crossover rate probability of mutation
    seed = 9;       % The random seed governing the random number generator
    N = 300;          % Number of individuals in the population
    numTrials = 20;  % Number of simulation trials per stimulus type
    numChirps = 5;  % Number of chirp stimulations in total
    maxIts = 150;     % Maximum number of iterations to run the DE rounds
    vpdq = 100;     % The q-value (temporal precision) used to get VPD 
    alpha = 0.01;   % Selective invariance (SI) score parameter
    neuron0 = TorusNeuronMod(); % A "blank slate" simulated neuron
    ss = 0.025;      % Step size in ms (McGillivray 2012)
    threshold = -25; % For spike time detection
    siglev = 50; % Requisite PSTH level for CSI significance
    numProcs = 1;         % Number of processors to run on
    
    % Biological Parameter Constraints %
    limits.sigmaB = [0.15 , 0.85];    % EI balance factor limits
    limits.f_h = [0 , 0];           % Ih current strength factor limits
    limits.f_ca = [0 , 0];          % Ica current strength factor limits
    limits.I_bias = [2 , -7];    % Ibias current level limits
    limits.N_xi = [1 , 1];      % Gaussian white noise intensity limits
    limits.Ws = [0 , 1.5];            % Synaptic current weight limits
    
                    %%%%%
    %%%%%%%%%%%%% ALGORTIHM %%%%%%%%%%%%%
                    %%%%%
    
    % Make Pool
    poolobject = gcp('nocreate');
    if isempty(poolobject) && ~(numProcs==1)
        parpool(numProcs);
    end
                    
    %%% Initialize population, input structures, and other parameters %%%
    
    % Create Real TS neuron
    % Note, the below are only needed if one wishes to override the fitness
    %   function with a different one based on comparison to a real neuron.
%    scRTS = load(TS_neuron_SC_file);
%    bcRTS = load(TS_neuron_BC_file);
     RTSN = []; % loadTSneuron(scRTS,bcRTS,numChirps);    

    % Build input data
    insyn = load(IN_SYN);
    synapses = insyn.synapseStruct;    
    % Initialize miscellaneous parameters 
    pop = InitializePopulation(limits, N);
    pop = parEvaluatePop(pop,[],synapses,vpdq,numChirps,numTrials,alpha,siglev);
    iterCounter = 0;
    rng(seed);
        
    % Until a goal is reached (maxiters reached, pop unchanged, invar > k):
    while iterCounter < maxIts
        % Clone the current population
            newPop = pop;
        
        % Construct fitness-biased probability distribution by creating an
        % array of probability values for choosing targets 
            probValArray = zeros(1,length(newPop));

            % Get Max fitness value (normalizer)
            statss = [newPop.stats];
            siss = [statss.SI];
            maxFitness = max(exp(-siss)); % Note worst fitness

            % Compute unnormalized probabilities
            for u = 1 : length(newPop)
                fitness = exp( -newPop(u).stats.SI );
                probValArray(u) = exp( (-fitness)/maxFitness ); 
            end

            % Normalize
            k = sum(probValArray);
            probValArray = probValArray / k;

        % Run in parallel across all agents:
        parfor aa = 1 : N    
            % Choose the current target
            target = newPop(aa);
            newAgent = target;
            % Pick three other agents, with a fitness-biased probability             
            	n = length(probValArray);
                borders = cumsum(probValArray);
                targs = [];
                while(length(targs) < 3)
                    r = rand;
                    a = -1;
                    for k = 1 : n
                        if r < borders(k); a = k; break; end
                    end
                    if a == -1; a = floor((n)*rand+1); end

                    if ~any(targs==round(a))
                        targs = [a targs]; %#ok<AGROW>
                    end
                end
                c = num2cell(pop(targs));
                [x1,x2,x3] = c{:};
            
            % Differentiation: Combine the three into a new agent
            	names = fieldnames(limits); % Cell array of field names
                Nv = length(names);
                for s = 1 : Nv
                    name = names{s};
                    newAgent.(name) = x1.(name) + F * ( x2.(name) - x3.(name) );
                    mx = max(limits.(name));
                    mn = min(limits.(name));
                    
                    % Handle boundary conditions by resampling
                    while ( (newAgent.(name) > mx) || (newAgent.(name) < mn) )
                        % Choose 3 agents (uniform distribution,NR)
                        ddd = 1 : N;
                        y = datasample(ddd,3,'Replace',false);
                        c = num2cell( pop(y) );
                        [x1,x2,x3] = c{:};
                        newAgent.(name) = x1.(name) + F * ( x2.(name) - x3.(name) );
                    end    
                    
                end
                newAgent.isEvaluated = false;
                newAgent.stats = struct;
            
            % Recombination: randomly mutate parts of the new agent
                for p = 1 : Nv
                    name = names{p};
                    r = rand;
                    if r > CR
                        newAgent.(name) = target.(name);
                    end
                end
                newAgent.isEvaluated = false;
                newAgent.stats = struct;
            
            % Evaluate the new agent
            %evaluate(agent,RTSN,synapses,q,nChirps,numTrials,alpha)
                % Preallocations
                psths = cell(1,numChirps);
                stsCa = cell(numChirps,numTrials);

                % Run Simulations
                for j = 1 : numChirps
                    neuron = neuron0.deepCopy();
                    for namae = 1 : Nv % Set the vals from the agent
                        targ = names{namae}; 
                        neuron.(targ) = newAgent.(targ);
                    end    
                    neuron.addAlphaSynapse(synapses(j).Esyn);
                    neuron.addAlphaSynapse(synapses(j).Isyn);
                    tmax = round(length(neuron.getAlphaSynapse(1).conved)/10);
                    psthTemp = zeros(tmax/ss,1);
                    % Run the multiple trials for the given input
                    for k = 1:numTrials % For all trials
                        % Run simulation
                        [~,v] = NSUtils.FastNeuronEulerMaruyama(neuron,ss,tmax); 
                        % Procure psth and binaries
                        [psth,~,sts] = NSUtils.getPsthStsAndBinaries(v,threshold,ss);
                        stsCa{j,k} = sts;
                        psthTemp = psthTemp + psth;
                    end
                    psths{j} = psthTemp/numTrials;
                    neuron.pruneSynapses();
                end

                % Compute post-sim statistics
                % May prefer to use matched trial vpd for speed
                stats = NSUtils.computePostSimulationStats...
                    (psths,numChirps,ss,stsCa,vpdq,numTrials,alpha,siglev);

                % Compare to real TS neuron (with full VPD, not matched trial)
                % Can use as a different fitness measure
%                 summ = 0; counter = 0;
%                 n = numChirps; m = numTrials;
%                 for i = 1 : n % across chirps
%                     for j = 1 : m % across trials
%                         sts1 = stsCa{i,j};
%                         for k = 1 : n % across all RTSN chirps
%                             for p = 1 : m % across all trials
%                                 sts2 = RTSN(k).rcell{p};
%                                 d = NSUtils.spkd(sts1,sts2,vpdq);
%                                 summ = summ + d;
%                                 counter = counter + 1;
%                             end
%                         end
%                     end
%                 end
%                 stats.vpd_RTSNavg = summ / counter;

                % Set newAgent stats and return
                newAgent.stats = stats;
                newAgent.isEvaluated = true;
            
            % Store Results
            newPop(aa) = newAgent;
            
        end
        
        % Perform Selection
        for i = 1 : length(pop)
            oldSI = pop(i).stats.SI;
            newSI = newPop(i).stats.SI;
            oldFitness = exp(-oldSI);
            newFitness = exp(-newSI);
            if ( (newFitness) < (oldFitness) )
                pop(i) = newPop(i);
            end
        end    
        
        % Save and Write current pop to file  
        savePop(pop,RTSN,iterCounter,tag); 
        iterCounter = iterCounter + 1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                             % SUBFUNCTIONS %                           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Spawn the initial population of parameter vectors %%%%%
function pop = InitializePopulation(limits, N) 
    names = fieldnames(limits); % Cell array of field names
    Nv = length(names);         % Size of DE individual parameter vector
    pop = [];   % The empty population variable
    
    % Safety checks
    if N < 3; error('Need at least 3 population members!'); end
    
    % For each population member:
    for i = N : -1 : 1 
        % For every parameter, generate a random value for it
        for j = 1 : Nv 
            name = names{j};
            mx = max(limits.(name));    % Max limit
            mn = min(limits.(name));    % Min limit
            r = mn + (mx-mn)*rand;      % Random number between max and min
            pop(i).(name) = r;          % Set value with random number
        end    
        % Set a flag to show the target has not been evaluated yet
        pop(i).isEvaluated = false;
        pop(i).stats = struct;
    end
end

%%%%% Evaluates an agent (i.e. a parameter vector) %%%%%
function newAgent = evaluate(agent,RTSN,synapses,q,nChirps,numTrials,alpha,siglev)
    % The new agent will be the old one, just with extra stats data
    newAgent = agent;
    
    % Safety Check 
    if agent.isEvaluated; return; end
    
    % Create neuron, sim parameters, and data storage
    neuron0 = TorusNeuronMod();
    ss = 0.025;      % Step size in ms (McGillivray 2012)
    threshold = -20; % For spike time detection
    
    % Preallocations
    psths = cell(1,nChirps);
    stsCa = cell(nChirps,numTrials);
    
    % Run Simulations
    for j = 1 : nChirps
        neuron = neuron0.deepCopy();
        neuron.addAlphaSynapse(synapses(j).Esyn);
        neuron.addAlphaSynapse(synapses(j).Isyn);
        tmax = round(length(neuron.getAlphaSynapse(1).conved)/10);
        psthTemp = zeros(round(length(neuron.getAlphaSynapse(1).conved)/10)/ss,1);
        % Run the multiple trials for the given input
        for k = 1:numTrials % For all trials
            % Run simulation
            [~,v] = NSUtils.FastNeuronEulerMaruyama(neuron,ss,tmax); 
            % Procure psth and binaries
            [psth,~,sts] = NSUtils.getPsthStsAndBinaries(v,threshold,ss);
            stsCa{j,k} = sts;
            psthTemp = psthTemp + psth;
        end
        psths{j} = psthTemp/numTrials;
        neuron.pruneSynapses();
    end
    
    % Compute post-sim statistics
    stats = NSUtils.computePostSimulationStats...
        (psths,nChirps,ss,stsCa,q,numTrials,alpha,siglev);
    
    % Compare to real TS neuron (with full VPD, not matched trial one)
%     sum = 0; counter = 0;
%     n = nChirps; m = numTrials;
%     for i = 1 : n % across chirps
%         for j = 1 : m % across trials
%             sts1 = stsCa{i,j};
%             for k = (i+1) : n % across all other chirps
%                 for p = 1 : m % across all trials
%                     sts2 = RTSN(k).rcell{p};
%                     d = NSUtils.spkd(sts1,sts2,q);
%                     sum = sum + d;
%                     counter = counter + 1;
%                 end
%             end
%         end
%     end
%     stats.vpd_RTSNavg = sum / counter;
     
    % Set newAgent stats and return
    newAgent.stats = stats;
    newAgent.isEvaluated = true;
end

%%%%% Performs the differentiation operator %%%%%
function newAgent = differentiate(x1, x2, x3, F, limits)
	names = fieldnames(limits); % Cell array of field names
    Nv = length(names);
	
    for s = 1 : Nv
		name = names{s};
		newAgent.(name) = x1.(name) + F * ( x2.(name) - x3.(name) );
        mx = max(limits.(name));
        mn = min(limits.(name));
        if newAgent.(name) > mx; 
            newAgent.(name) = mx;
        elseif newAgent.(name) < mn; 
            newAgent.(name) = mn;
        end
    end
    newAgent.isEvaluated = false;
    newAgent.stats = struct;
end

%%%%% Performs the recombination operator %%%%%
function newAgent = recombine(newAgent, target, CR, limits)
	names = fieldnames(limits); % Cell array of field names
    Nv = length(names);
    for p = 1 : Nv
		name = names{p};
		r = rand;
        if r > CR
			newAgent.(name) = target.(name);
        end
    end
    newAgent.isEvaluated = false;
    newAgent.stats = struct;
end


%%%%% Creates the array of probability values for choosing targets %%%%%
function probValArray = makeProbValArray(newPop)
    probValArray = zeros(1,length(newPop));
    
    % Get Max fitness value (normalizer)
    statss = [newPop.stats];
    vs = [statss.vpd_RTSNavg];
    maxFitness = max(vs); %Note worst fitness
    
    % Compute unnormalized probabilities
    for u = 1 : length(newPop)
        fitness = newPop(u).stats.vpd_RTSNavg;
        probValArray(u) = exp( (-fitness)/maxFitness ); 
    end
    
    % Normalize
    k = sum(probValArray);
    probValArray = probValArray / k;
    
end

%%%%% Picks agents from the population based on biased probdist %%%%%
function [x1,x2,x3] = pickAgents(pop, probValArray)
	n = length(probValArray);
	borders = cumsum(probValArray);
	targs = [];
	
	while(length(targs) < 3)
		r = rand;
		a = -1;
		for k = 1 : n
			if r < borders(k); a = k; break; end
		end
		if a == -1; a = floor((n)*rand+1); end
		
		if ~any(targs==round(a))
			targs = [a targs]; %#ok<AGROW>
		end
	end
	
	c = num2cell(pop(targs));
	[x1,x2,x3] = c{:};
end

%%%%% Saves the current generation %%%%%
function savePop(pop,RTSN,iterCounter,tag) %#ok<INUSL>
    generation = iterCounter; %#ok<NASGU>
    population = pop;  %#ok<NASGU>
    % Avg and max vpd_RTSNavg and SI
    %popstats.minVpdRtsn = inf;
    popstats.maxSI = 0;
    %sumVpdRtsn = 0; 
    sumSI = 0;
    for p = 1 : length(pop)
        targ = pop(p);
%        sumVpdRtsn = sumVpdRtsn + targ.stats.vpd_RTSNavg;
        sumSI = sumSI + targ.stats.SI;
%         if targ.stats.vpd_RTSNavg < popstats.minVpdRtsn
%             popstats.minVpdRtsn = targ.stats.vpd_RTSNavg;
%         end    
        if targ.stats.SI > popstats.maxSI
            popstats.maxSI = targ.stats.SI;
        end    
    end    
%    popstats.avgVPD_RTSNavg = sumVpdRtsn / length(pop);
    popstats.avgSI = sumSI / length(pop);
    s = sprintf('GenDE-%d-%s',iterCounter,tag);
    save(s,'generation','population','RTSN','popstats');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                           % HELPER FUNCTIONS %                           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Helper for building the combined synapse structure %%%%%
function synapseStruct = buildSynapses(chirps) %#ok<DEFNU>
    for i = 1:length(chirps)
        s = sprintf('Chirp%d-E_avg',i); 
        synapseStruct(i).Esyn = ...
            AlphaSynapse(s,chirps(i).avgEpsth,chirps(i).chirp,'E'); %#ok<AGROW>
        s = sprintf('Chirp%d-I_avg',i); 
        synapseStruct(i).Isyn = ... 
            AlphaSynapse(s,chirps(i).avgIpsth,chirps(i).chirp,'I'); %#ok<AGROW>
    end  
end

%%%%% Loads the data for the real Torus neuron into a struct %%%%%
function ts = loadTSneuron(scRTS,bcRTS,numChirps)
    ts = [];
    for i = numChirps : -1 : 1
        if round(i)==5
            ts(i).rcell = bcRTS.rcell;
            ts(i).psth = bcRTS.psth;
        else
            ts(i).rcell = scRTS.rcells(i).rcell;
            ts(i).psth = scRTS.psth(i).psthA;
        end
    end
end

%%%%% Helper function to evaluate the full population %%%%%
function newpop = parEvaluatePop(pop,RTSN,synapses,vpdq,nChrps,nTrls,alpha,siglev)
    newpop = pop;
    parfor ii = 1:length(pop)
        agent = pop(ii);
        newpop(ii) = evaluate(agent,RTSN,synapses,vpdq,nChrps,nTrls,alpha,siglev);
    end
end


