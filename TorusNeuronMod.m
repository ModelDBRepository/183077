classdef TorusNeuronMod < handle % Modify class state by reference
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%                            NEURON PROPERTIES                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	properties(GetAccess = public, SetAccess = public)
	% % Neuron Biophysical Parameters % %
        % Parameters from Khosravi-Hashemi, 2012 and 2011 
            % Includes leak, delayed K-rectifier, and Na
            % Experimental verification from Chacron and Fortune, 2010 %%
            Cm = 1.0;   % microFarads/cm^2
            gNa = 30;   % microS
            gK = 10;    % microS
            gL = 0.18;  % microS
            ENa = 60;   % mV 
            EK = -85;   % mV
            EL = -65;   % mV (leak)
        % Describes Ih current from HCN channels %
            % Based on Migliore, 2012; Poolos, 2002; Neymotin, 2013 %
            % Note: we ignore the t & v independent Ilk current as of now
            f_h = 1;    % scaling factor to vary Ih strength (Neymotin goes between 0 and 2)
            Eh = -30;   % mV
            g_h = 7;    % maximal conductance density (This value from Migliore 2012 Table 1)
        % Calcium Channel (T-type) Parameters %%
            f_ca = 1; 
            ECa = 120;
            gT_Ca = 0.32; % microS
            tau_hCa = 30;
        % Bias Current Level %
            I_bias = -1.3; % -1.3, Khosravi-Hashemi, 2012; -4.3, Khosravi-Hashemi, 2011
        % Stochastic Gaussian White Noise Current Parameters %
            % See Khosravi-Hashemi, 2011; Guo, 2011
            gaussMean = 0;      % Mean value for stochastic process
            gaussVariance = 0.8;  % Variance of process (from McGillivray, 2012)
            N_xi = 1;   % Scaling factor for process (Gaussian Intensity)
        % Synaptic properties
            Ws = 1;     % The overall synaptic weight (scales Isyn after summing)
            sigmaB = 0.5; % Balance factor between E and I cell input strength (1=all E)
    end
    
    properties(GetAccess = public, SetAccess = private)
        alphaSynapses = []; % List of Input alpha synapses
        vi;                 % mV (resting). Always reset.
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%                    PUBLIC NON-STATIC METHODS                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	methods(Static = false, Access = public)
        
    % MAIN VOLTAGE DERIVATIVE FUNCTION %
        function dvdt = voltageDerivative(this,v,n,h,hca,t)
            
            % Depolarizing Sodium Current %
            am = TorusNeuronMod.alphaM(v);
            minf =  am/(am + TorusNeuronMod.betaM(v));
            I_Na = this.gNa*(minf^3)*(0.85-n)*(v-this.ENa);
            
            % Delayed rectifier Potassium Current %
            I_K = this.gK*(n^4)*(v-this.EK);
            
            % HCN Channel-mediated Hyperpolarization Activated Current %
            I_H = (this.g_h*this.f_h*h) * (v-this.Eh);
            
            % Leak Current %
            I_L = this.gL*(v-this.EL);
            
            % T-type Calcium Current %
            sinf = 1/(1+exp(-1*(v+69)/7.8));
            I_Ca = this.gT_Ca*(sinf^3)*hca*(v-this.ECa)*this.f_ca;
            
            % Synaptic Input Current %
            I_syns = zeros(length(this.alphaSynapses),1); 
            as = this.alphaSynapses;
            for i = 1:length(as)
                if strcmpi(as(i).type,'E')
                    I_syns(i) = (2*this.sigmaB) * as(i).I_syn(v,t);
                else
                    I_syns(i) = (2*(1-this.sigmaB)) * as(i).I_syn(v,t);
                end
            end 
            I_syn = this.Ws * sum(I_syns);
            
            % Sum Currents (adds negative sign here)
            I_total = -I_Na - I_K - I_H - I_L - I_Ca - I_syn + this.I_bias;
            
            % Get final dv/dt
            dvdt = I_total/this.Cm;
        end
        
    % % INITIATION HELPER % %
        % Returns the initial (Steady-state) vector describing the system %%
        function iv = getInitialVector(this) 
            this.resetRestingVi();
            v = this.vi;
            n = this.nInf(v);
            h = this.hInf(v);
            hCa = this.hInf_Ca(v);
            iv = [v,n,h,hCa]';
        end
        
    % % HELPER FUNCTIONS FOR SYNAPTIC INPUT % %
        function addAlphaSynapse(this,alphasyn)
            if ~isa(alphasyn,'AlphaSynapse'); 
                error('Adding non-alpha-synapse to neuron!');
            end
            this.alphaSynapses = [this.alphaSynapses alphasyn];
        end
        function as = getAlphaSynapse(this,num)
            as = this.alphaSynapses(num);
        end    
        function pruneSynapses(this)
            this.alphaSynapses = [];
        end
        
	% % CONSTRUCTOR % %
		function t = TorusNeuronMod(varargin) %Permits 0-arg constructor calls
           if ~any([nargin==0,nargin==3])
               error('Improper argument number to neuron constructor')
           end
           if nargin == 3 
                cm = varargin{1};Gs=varargin{2};Es=varargin{3};
                if length(Gs)~=3 || length(Es)~=3
                    error('G and E vector length too short')
                end
                t.Cm = cm; 
                t.gNa = Gs(1); t.gK = Gs(2); t.gL = Gs(3);
                t.ENa = Es(1); t.EK = Es(2); t.EL = Es(3);
                t.resetRestingVi();
           end
        end
        
    % % DEEP COPY METHOD % %
        function newNeuron = deepCopy(this)
            warning off MATLAB:structOnObject % Suppress warnings
            if ~isa(this,'TorusNeuronMod')
                error('Unexpected Variable Type');
            end
            newNeuron = feval(class(this)); % Instantiate class
            p = fieldnames(struct(this));   % Procure all field names (incl. private)
            for i = 1:length(p)
                newNeuron.(p{i}) = this.(p{i}); % Set new neuron properties equivalently
            end
            for j = 1:length(this.alphaSynapses)
                newNeuron.alphaSynapses(j) = this.alphaSynapses(j).deepCopy();
            end
            warning on MATLAB:structOnObject % Turn warnings back on
        end
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%                    PRIVATE NON-STATIC METHODS                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods(Static = false, Access = private)
  
        % Reset estimated proper initial membrane voltage by running short sim %
        % No noise and no synaptic input. Basic Euler method to solve the
        % DE.
        function resetRestingVi(this) 
            % Ignore synaptic conductances. 
            ss = 0.025; % in ms
            endTime = 150; % in ms
            ts = ss:ss:endTime;
            N = length(ts);
            v0 = -55;          % Initial guess for v
            y0 = [v0;this.nInf(v0);this.hInf(v0);this.hInf_Ca(v0)];
            y = [ y0 zeros(4,N-1) ];
            for i = 2 : N
               % t = ts(i);
                v = y(1,i-1); %vec(1);
                n = y(2,i-1); %vec(2);
                h = y(3,i-1); %vec(3);
                hca = y(4,i-1); %vec(4);
                %vp = [0;0;0;0];
                %%%% Get dv/dt %%%%
                % Depolarizing Sodium Current %
                amh = 0.1*(v+40.7);
                am = (amh)/(1-exp(-amh));
                bm = 4*exp(-0.05*(v+49.7));
                minf =  am/(am + bm);
                I_Na = this.gNa*(minf^3)*(0.85-n)*(v-this.ENa);
                % Delayed rectifier Potassium Current %
                I_K = this.gK*(n^4)*(v-this.EK);
                % HCN Channel-mediated Hyperpolarization Activated Current %
                I_H = (this.g_h*this.f_h*h) * (v-this.Eh);
                % Leak Current %
                I_L = this.gL*(v-this.EL);
                % T-type Calcium Current %
                sinf = 1/(1+exp(-1*(v+69)/7.8));
                I_Ca = this.gT_Ca*(sinf^3)*hca*(v-this.ECa)*this.f_ca;
                % Synaptic Input Current %
                %ind = round(t/ss);
                %sb = this.sigmaB;
                I_syn = 0; %this.Ws * (2*(v-esAmpa))*((asynValsE(ind)*(sb))+(asynValsI(ind)*(1-sb)));
                % Sum Currents (adds negative sign here)
                I_total = -I_Na - I_K - I_H - I_L - I_Ca - I_syn + this.I_bias;
                % Get final dv/dt
                dvdt = I_total/this.Cm;
                %vp(1) = dvdt;
                %%%% End get dv/dt %%%%

                %%%% Get N' %%%%
                an = (0.01*(v+40.7))/(1-exp(-0.1*(v+40.7)));
                bn = 0.125*exp(-0.0125*(v+50.7));
                ninf = an/(an+bn);
                taun = 0.05/(an+bn);
                ndot = (ninf-n)/taun;
                %vp(2) = ndot;
                %%%% End Get N' %%%%

                %%%% Get h' %%%%
                tauh = exp(0.033*(v+75))/( 0.011*( 1 + exp(0.083*(v+75)) ) );
                v_hm = -73; % half maximal voltage in time constrant of h
                hinf = 1/(1+exp(0.151*(v-v_hm)));
                hdot = (hinf-h)/tauh;
                %vp(3) = hdot;
                %%%% End Get h' %%%%

                %%%% Get T' (Ca) %%%%
                q = sqrt(0.25+exp( (v+82)/6.3 ) );
                hinfca = 1/(0.5+q);
                hCadot = 2*(hinfca - hca)/this.tau_hCa;
                %vp(4) = hCadot;
                %%%% End Get T' (Ca) %%%%

                %%%% Combine the vectors %%%%
                vp = [dvdt; ndot; hdot; hCadot];
                
                %%%%%%%%%%%%%%%%%%%%%% Take the step %%%%%%%%%%%%%%%%%%%%%%
                y(:,i) = y(:,i-1) + (ss * vp);
            end
            this.vi = y(1,end);
        end
        
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%                    PUBLIC STATIC METHODS                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
	methods (Static = true, Access = public)
        
    % PRIMARY VECTOR OF DERIVATIVES FOR ODE SOLVER TO USE %     
    function vp = vprime(t,vec,this)
        % Preparation %
        %[v,n,h,hca]=deal(vec(1),vec(2),vec(3),vec(4));
        v = vec(1);
        n = vec(2);
        h = vec(3);
        hca = vec(4);
        vp = zeros(4,1);
        
        %%%% Get dv/dt %%%%
        % Depolarizing Sodium Current %
        amh = 0.1*(v+40.7);
        am = (amh)/(1-exp(-amh));
        bm = 4*exp(-0.05*(v+49.7));
        minf =  am/(am + bm);
        I_Na = this.gNa*(minf^3)*(0.85-n)*(v-this.ENa);
        % Delayed rectifier Potassium Current %
        I_K = this.gK*(n^4)*(v-this.EK);
        % HCN Channel-mediated Hyperpolarization Activated Current %
        I_H = (this.g_h*this.f_h*h) * (v-this.Eh);
        % Leak Current %
        I_L = this.gL*(v-this.EL);
        % T-type Calcium Current %
        sinf = 1/(1+exp(-1*(v+69)/7.8));
        I_Ca = this.gT_Ca*(sinf^3)*hca*(v-this.ECa)*this.f_ca;
        % Synaptic Input Current %
       % as = this.alphaSynapses;
        % Get Isyns. ASSUMES ONLY TWO SYNAPSES! AND ASSUMES BOTH ARE SAME AMPAR! %
        Is = 0; %#ok<NASGU>
        if strcmpi(this.alphaSynapses(1).type,'E')
            ind = round(t/this.alphaSynapses(1).ss2);
            sb = this.sigmaB;
            esAmpa = this.alphaSynapses(1).E_syn;
            Is = (2*(v-esAmpa))*(this.alphaSynapses(1).preCalcInterpVals(ind)*(sb)...
               + this.alphaSynapses(2).preCalcInterpVals(ind)*(1-sb));
            %Is = (this.alphaSynapses(1).preCalcInterpVals(ind)*(v-esAmpa)*(2*sb)...
            %    + this.alphaSynapses(2).preCalcInterpVals(ind)*(v-esAmpa)*(2*(1-sb)));
            %I_syns(1) = (2*this.sigmaB) * as(1).I_syn(v,t);
            %I_syns(2) = (2*(1-this.sigmaB)) * as(2).I_syn(v,t);
        else
            ind = round(t/this.alphaSynapses(1).ss2);
            sb = this.sigmaB;
            esAmpa = this.alphaSynapses(1).E_syn;
            Is = this.alphaSynapses(2).preCalcInterpVals(ind)*(v-esAmpa)*(2*sb)...
                + this.alphaSynapses(1).preCalcInterpVals(ind)*(v-esAmpa)*(2*(1-sb));
        end
        I_syn = this.Ws * Is;
        % Sum Currents (adds negative sign here)
        I_total = -I_Na - I_K - I_H - I_L - I_Ca - I_syn + this.I_bias;
        % Get final dv/dt
        dvdt = I_total/this.Cm;
        vp(1) = dvdt;
        %%%% End get dv/dt %%%%
        
        %%%% Get N' %%%%
        an = (0.01*(v+40.7))/(1-exp(-0.1*(v+40.7)));
        bn = 0.125*exp(-0.0125*(v+50.7));
        ninf = an/(an+bn);
        taun = 0.05/(an+bn);
        ndot = (ninf-n)/taun;
        vp(2) = ndot;
        %%%% End Get N' %%%%
    
        %%%% Get h' %%%%
        tauh = exp(0.033*(v+75))/( 0.011*( 1 + exp(0.083*(v+75)) ) );
        v_hm = -73; % half maximal voltage in time constrant of h
        hinf = 1/(1+exp(0.151*(v-v_hm)));
        hdot = (hinf-h)/tauh;
        vp(3) = hdot;
        %%%% End Get h' %%%%
    
        %%%% Get T' (Ca) %%%%
        q = sqrt(0.25+exp( (v+82)/6.3 ) );
        hinfca = 1/(0.5+q);
        hCadot = 2*(hinfca - hca)/this.tau_hCa;
        vp(4) = hCadot;
        %%%% End Get T' (Ca) %%%%
    end
    % END LOW OVERHEAD VECTOR DERIVATIVE EVALUATOR %
    
    % (Older Version) %
    function vp = vprimeOld(t,vec,cell)
        [v,n,h,hca]=deal(vec(1),vec(2),vec(3),vec(4));
        vp = zeros(4,1);
        % Compute vector of derivatives
        vp(1) = cell.voltageDerivative(v,n,h,hca,t);
        vp(2) = cell.nPrime(v,n);
        vp(3) = cell.hPrime(v,h);
        vp(4) = cell.hCaPrime(v,hca,cell);
    end

	% SECONDARY DIFFERENTIAL FUNCTIONS (GATING AND ACTIVATION VARIABLES) %
    % From Khosravi-Hashemi, 2012 and 2011 %%
        function ndot = nPrime(v,n)
            ndot = (TorusNeuronMod.nInf(v)-n)/TorusNeuronMod.tau_n(v);
        end
        function ninf = nInf(v)
            an = TorusNeuronMod.alphaN(v);
            ninf = an/(an+TorusNeuronMod.betaN(v));
        end
        function taun = tau_n(v)
            taun = 0.05/(TorusNeuronMod.alphaN(v)+TorusNeuronMod.betaN(v));
        end
        function alphan = alphaN(v)
            alphan = (0.01*(v+40.7))/(1-exp(-0.1*(v+40.7)));
        end
        function betan = betaN(v)
            betan = 0.125*exp(-0.0125*(v+50.7));
        end
        function minf = mInf(v)
            am = TorusNeuronMod.alphaM(v);
            minf = am/(am+TorusNeuronMod.betaM(v));
        end
        function alpham = alphaM(v)
            alpham = (0.1*(v+40.7))/(1-exp(-0.1*(v+40.7)));
        end
        function betam = betaM(v)
            betam = 4*exp(-0.05*(v+49.7));
        end
    % Largely from Neymotin 2013 (See also Poolos 2002, Migliore 2012) %%
        function hdot = hPrime(v,h) 
            hdot = (TorusNeuronMod.hInf(v)-h)/TorusNeuronMod.tau_h(v);
        end
        function hinf = hInf(v)
            v_hm = -73; % half maximal voltage in time constrant of h
            hinf = 1/(1+exp(0.151*(v-v_hm)));
        end
        function tauh = tau_h(v)
           tauh = exp(0.033*(v+75))/( 0.011*( 1 + exp(0.083*(v+75)) ) ); 
        end
    % Calcium Channel Equations %%
        function hCadot = hCaPrime(v,hCa,neuronInstance)
            hCadot = 2*(TorusNeuronMod.hInf_Ca(v) - hCa)/neuronInstance.tau_hCa;
        end
        function sinf = sInf(v)
            sinf = 1/(1+exp(-1*(v+69)/7.8));
        end
        function hinfca = hInf_Ca(v)
            q = sqrt(0.25+exp( (v+82)/6.3 ) );
            hinfca = 1/(0.5+q);
        end
    end
    %#ok<*DSPS>
end   


    
