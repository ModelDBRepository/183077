classdef AlphaSynapse < handle

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         ALPHA SYNAPSE PROPERTIES                  	 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	properties(Access = public)
		% Default Synapse: AMPAR mediated excitatory synapse %
		tau_syn = 20; %ms, as in Khosravi-Hashemi, 2011. 
		g_syn_m = 0.13; % Maximal AMPAR conductance value (mS/cm^2) [jensen,2001]
		E_syn = 0; % AMPAR reversal potential
		%W_syn = 1; % Synaptic weight (strength) factor % Accounted for in TSNM class
        conved = []; % Convolved firing pattern
        type = ''; % Neuron type (E or I, usually)
        % Pre-calculated, interpolated vector of current values %
        preCalcInterpVals = [];
    end
	properties(GetAccess = public,SetAccess = private)
        % Used for direct convolution of PSTH with Alpha Function % 
            psth = [];
        % Identifier for the synapse %
            name = []; 
        % Scaling factor for psth %
            conductanceScalingFactor = 0.0005;
        % Chirp stim %
            chirp = [];
        % Step Sizes %    
            ss1 = 0.1; % Step size of input data
            ss2 = 0.025; % Step size of simulation
    end
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     ALPHA SYNAPSE NON-STATIC METHODS                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
	methods(Static = false, Access = public)
		% % % % Constructor % % % % 
		function t = AlphaSynapse(s,psthIn,chirp,type)
            t.name = s; % name of synapse (e.g. E cell input)
			if isempty(psthIn);error('No psth given to constructor');end
            if isempty(type);error('No type given to constructor');end
			t.psth = psthIn;
            t.conved = t.produceConductanceByConvolution();
            t.preCalcInterpVals = t.precomputeInterpolatedConds(t.ss1,t.ss2);
            t.chirp = chirp;
            t.type = type;
        end
		
        % Precalculates the interpolated values of the synaptic
        % conductance
        function interpedVals = precomputeInterpolatedConds(as,ss1,ss2)
            times1 = ss1:ss1:length(as.conved)*ss1;
            origVals = as.conved;
            times2 = ss2:ss2:length(as.conved)*ss1; %1600004
            interpedVals = interp1(times1,origVals,times2,'linear');
            interpedVals(isnan(interpedVals)) = 0; % FIXES EDGE CASES
        end
            
		% Compute synaptic current at time t
		function isyn = I_syn(this,v,t)
			isyn = this.preCalcInterpVals(round(t/this.ss2)) * ...
                (v - this.E_syn);
		end
			
        % Feeding in a psth with 0.1ms time steps
        % Alpha function by Rall 1967 
        % Convolution as in Khosravi-Hashemi, 2011
		function g_syn = produceConductanceByConvolution(this)
            % Produce Alpha function
            v = (1:length(this.psth));
            alphaVec = (v/this.tau_syn).*exp(1 - v/this.tau_syn);
            % Convolve with psth
            gSynConv = conv(this.psth,alphaVec,'full');
            g_syn = this.g_syn_m * gSynConv(1:length(this.psth))* ...
                this.conductanceScalingFactor;
        end
		
		% Simple linear interpolation by averaging between the two values  
		% The times will come in milliseconds
		% Our convolved array is in tenths of milliseconds
		% Hence we multiply by ten to get the proper array value
% 		function gs = calcGsyn(this,t)
%             t1 = floor(t*10);
%             t2 = ceil(t*10);
%             if t1==t2 || t1==0; 
%                 gs = this.conved(t2);
%                 return
%             end
%             leftVal = this.conved(t1);
% 			rightVal = this.conved(t2);
%             deltaT = (t*10) - t1;
%             slope = (rightVal-leftVal)/(t2-t1);
%             gs = leftVal + (slope*deltaT);
%         end
        
        % % % % HELPER FUNCTIONS % % % %
        % Returns an array of conductances for an array of times
        function gSs = gSynArray(this,t)
            gSs = zeros(length(t),1);
            for i = 1:length(t)
                gSs(i) = this.calcGsyn(t(i));
            end
        end
        % Returns the name of this structure
        function n = getName(this)
            n = this.name;
        end
        % Get length of the stored PSTH
        function l = getPsthLength(this)
            l = length(this.psth);
        end  
        % Get the chirp stim this synapse represents
        function c = getChirp(this)
            c = this.chirp;
        end
        % Deep copies the current synapse
        function newAS = deepCopy(this)
            warning off MATLAB:structOnObject % Suppress warnings
            if ~isa(this,'AlphaSynapse')
                 error('Unexpected Variable Type');
            end               
            newAS = feval(class(this));
            p = fieldnames(struct(this));   % Procure all field names (incl. private)
            for i = 1:length(p)
                newAS.(p{i}) = this.(p{i}); % Set new neuron properties equivalently
            end
            warning on MATLAB:structOnObject % Turn warnings back on
        end
    end
end