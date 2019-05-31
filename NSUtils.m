classdef NSUtils 
	methods (Static = true, Access = public)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        % % % % % % % % % % SIMULATION UTILITIES % % % % % % % % %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Master calling function for calculating the CSI and VPDs of the
        % just-run simulation. Note that it should compute VPD_full.
        % The RMSE will also be computed by xcorr and circshift
        function stats = computePostSimulationStats...
                (psths,nChirps,ss,stsCa,q,numTrials,alpha,siglev)
            [allcsias,allcsifs,avgFRc] = NSUtils.computeAvgCSIS(psths,nChirps,ss,siglev); 
            avgvpd = NSUtils.computeTotalAvgVPD(stsCa,q,nChirps,numTrials);
            stats.csias = allcsias; % CSIavgs
            stats.csifs = allcsifs; % CSImaxfrs
            stats.vpd_mt = avgvpd;  % Average cross-trial VPD
            stats.avgFRc = avgFRc;  % Average FR during the window after the chirp
            % Compute SI score
            SI = mean(allcsifs) - (alpha * avgvpd);
            if SI < 0; SI = 0; end
            stats.SI = SI;
            % Store params
            stats.params.q = q;
            stats.params.alphaSI = alpha;
            % Compute pairwise RMSE average
            rmse = 0; counter = 0;
            for c1 = 1 : nChirps
                for c2 = (c1+1) : nChirps
                    p1 = psths{c1};
                    p2 = psths{c2};
                    rmse = NSUtils.calcOptRMSE(p1,p2);
                    counter = counter + 1;
                end
            end
            stats.rmse_avgpw = rmse / counter;
        end
        
        % Input: Two PSTHs
        % Output: the RMSE between them, after optimal circ-shifting
        function rmse = calcOptRMSE(psth1,psth2)
            N = length(psth1);
            % Compute cross-correlation
            [cc,lags] = xcorr(psth1,psth2);
            shift = lags(cc==max(cc));
            % Circular shift the data (psth2) to match
            newpsth2 = circshift(psth2,shift);
            % Compute the RMSE now
            rmse = sqrt( sum( (psth1(:)-newpsth2(:)).^2 ) / N );
        end   

        % Compute the cross-trial chirp selectivities.
        % Computes the CSIAs, which use the average FRs across the beat and
        % chirp, and the CSIFs, which use the max FRs across the beat and
        % chirp. 
        function [allcsias,allcsifs,cFRavg] = computeAvgCSIS(psths,numChirps,ss,siglev)
            n = numChirps;
            allcsias = zeros(n,1);
            allcsifs = zeros(n,1);
            cFRavg = 0;
            offset = 5; % Because inexactness of the spike time (i.e.
               % to before the chirp time), ignore offset ms before it,
               % else the timing detection may be off as well,
               % especially for the BC. Spike bleeding post-smoothing also
               % requires this.
            %%% Prepare input variables %%%
            % % % % Fixed Input variables % % % % 
            TimeWindow = 100;  % ms (for a 5hz stimului), after the chirp onset
            endTime = 1600;   % The time length of the stimulus in ms
            % % % % Compute necessary variables % % % %
            ChirpTimeOnset = endTime/2 - offset; % How many ms into the stim the chirp occurs
            % % % % Compute actual indices for chirp window % % % %
            t = ss:ss:endTime;
            ind1 = find(abs(t-ChirpTimeOnset)<0.0001);
            ind2 = ind1+(TimeWindow/ss); % Time of ind1 + window (in time points)
            for i = 1 : n
                % Get PSTH %
                psth = psths{i};
                
                %%% Update the FRavg during the post-chirp window %%%
                cFRavg = ( cFRavg*(i-1) + mean(psth(ind1:ind2)) )/i;
                
                %%% Compute individual CSIA %%%
                % % % % Get R_beat % % % %
                % The average FR of the neuron across the whole PSTH, except
                % for the window after the chirp onset (i.e. for TimeWindow ms
                % after ChirpTimeOnset) in two segments (pre/post chirp).
                % Using averages weighted by the number of points in the two
                % intervals, so as not to overwight mean2, albeit slightly.
                wmean1 = mean(psth(1:ind1-1))*length(psth(1:ind1-1));
                wmean2 = mean(psth(ind2+1:end))*length(psth(ind2+1:end));
                rb = (wmean1+wmean2)/(length(psth)-length(ind1:ind2));
                % % % % Get R_Chirp % % % %
                rc = mean(psth(ind1:ind2));
                % % % % Get CSI and Store It % % % %
                if ( (rc+rb)==0 )
                    csia = 0; % No firing at occurred, so there is no selectivity
                else
                    csia = ( (rc-rb) / (rc+rb) );
                end
                allcsias(i,1) = csia;
                
                %%% Compute individual CSIF %%%
                % % % % Get R_beat % % % %
                % The max FR of the neuron anywhere outside the TimeWindow,
                % looking both before and after the chirp:chirp+TW period 
                mfr1 = max(psth(1:ind1-1));
                mfr2 = max(psth(ind2+1:end));
                rb = max(mfr1,mfr2);
                % % % % Get R_Chirp % % % %
                rc = max(psth(ind1:ind2));
               % fprintf('rb=%f,rc=%f\n',rb,rc);
                % % % % Get CSI and Store It % % % %
                if ( (rc+rb)==0 )
                    csif = 0; % No firing at occurred, so there is no selectivity
                else
                    csif = ( (rc-rb) / (rc+rb) );
                end
                
                % % Check if the csi passes our significance requirements % %
                if rc < siglev; csif = 0; end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                allcsifs(i,1) = csif;
            end
        end
        
        % Computes the all-trial average Victor-Purpura distance.
        %   In detail, it calculates:
        %     vpdavg = (1/n) sum_(forall chirps) sum_(forall trials) d( m(i,j),m(k,p) )
        %   where m(u,v) is the spike-time vector for chirp u & trial v, d
        %   is the VPD function, and n is the number of combinations checked.
        % In essence, we compute the VPD between all possible pairs of
        %   chirps and trials, sum them, and take the average.
        % In case the binaries hold responses to non-chirp
        %   input as well, the function asks for the number of chirps and
        %   assumes they are the first numChirps inputs (i.e. binaries).
        % Input:
        %   sts = spike times, as a |nchirps|x|ntrials| cell array (stsCa)
        %   q = the VPD q parameter
        %   numchirps = the total number of chirps to look across (in case
        %       the sts contains more than just the stims we want to analyze)
        %   numTrials = the number of trials for which each stim was run
        function avg = computeTotalAvgVPD(sts,q,numChirps,numTrials)
            n = numChirps; m = numTrials;
            sum = 0; counter = 0;
            for i = 1 : n % across chirps
                for j = 1 : m % across trials
                    sts1 = sts{i,j};
                    for k = i : n % across all other chirps
                        Lij = 1;
                        if i==k; Lij = j+1; end
                        for p = Lij : m
                            sts2 = sts{k,p};
                            d = NSUtils.spkd(sts1,sts2,q);
                            sum = sum + d;
                            counter = counter + 1;
                        end
                    end
                end
            end
            avg = sum / counter;
        end
        
        % Computes the total VPD (then averages it) but only considers
        % matched trials (i.e. computes d across all chirps c_i, but only
        % compares to the same-number trial, i.e:
        %   d_avg = (1/n) sum_{c in chirps} sum_{trial j} sum_{other
        %       chirps} d( M(ci,j), M(ci,j) )
        % Necessary when computational time is problematic.
        function avg = computeMatchedTrialAvgVPD(sts,cost,numChirps,numTrials)
            n = numChirps; m = numTrials;
            sum = 0; counter = 0;
            for i = 1 : n % across chirps
                for j = 1 : m % across trials
                    tli = sts{i,j};
                    for k = (i+1) : n % across all other chirps
                        tlj = sts{k,j};
                        %%% RUN VPD CALCULATOR %%%
                        %d = NSUtils.spkd(sts1,sts2,q);
                        nspi=length(tli);
                        nspj=length(tlj);
                        if cost==0
                           d=abs(nspi-nspj);
                           sum = sum + d;
                           counter = counter + 1;
                           continue;
                           %return
                        elseif cost==Inf
                           d=nspi+nspj;
                           sum = sum + d;
                           counter = counter + 1;
                           continue;
                           %return
                        end
                        scr=zeros(nspi+1,nspj+1);
                        %INITIALIZE MARGINS WITH COST OF ADDING A SPIKE
                        scr(:,1)=(0:nspi)';
                        scr(1,:)=(0:nspj);
                        if nspi && nspj 
                           for ii=2:nspi+1
                              for jj=2:nspj+1
                                 scr(ii,jj)=min([scr(ii-1,jj)+1 scr(ii,jj-1)+1 ...
                                     scr(ii-1,jj-1)+cost*abs(tli(ii-1)-tlj(jj-1))]);
                              end
                           end
                        end
                        d=scr(nspi+1,nspj+1); 
                        %%% END VPD CALCULATOR %%%
                        sum = sum + d;
                        counter = counter + 1;
                    end
                end
            end
            avg = sum / counter;            
        end
        
        %%%%%% Begin External Code %%%%%%
        % The Victor-Purpura Distance between a spike train pair
        % Procured on 5/26/14 at:
        % http://www-users.med.cornell.edu/~jdvicto/spkdm.html
        % Note: this is not the most computationally efficient method
        function d=spkd(tli,tlj,cost)
            % d=spkd(tli,tlj,cost) calculates the "spike time" distance
            % (Victor & Purpura 1996) for a single cost
            % tli: vector of spike times for first spike train
            % tlj: vector of spike times for second spike train
            % cost: cost per unit time to move a spike
            %  Copyright (c) 1999 by Daniel Reich and Jonathan Victor.
            %  Translated to Matlab by Daniel Reich from FORTRAN code by Jonathan Victor.
            nspi=length(tli);
            nspj=length(tlj);
            if cost==0
               d=abs(nspi-nspj);
               return
            elseif cost==Inf
               d=nspi+nspj;
               return
            end
            scr=zeros(nspi+1,nspj+1);
            %INITIALIZE MARGINS WITH COST OF ADDING A SPIKE
            scr(:,1)=(0:nspi)';
            scr(1,:)=(0:nspj);
            if nspi & nspj %#ok<AND2>
               for i=2:nspi+1
                  for j=2:nspj+1
                     scr(i,j)=min([scr(i-1,j)+1 scr(i,j-1)+1 ...
                         scr(i-1,j-1)+cost*abs(tli(i-1)-tlj(j-1))]);
                  end
               end
            end
            d=scr(nspi+1,nspj+1);        
        end
        %%%%%% End External Code %%%%%%        
        

    % Simple method to streamline access to the psth, spike times,
        % and binaries after running a simulation.
        function [psth,bins,sts]=getPsthStsAndBinaries(vals,thresh,ss)
            vals = vals';
            vs = vals(:,1);
            mt = length(vs)*ss;
            times = ss:ss:mt;
            bins = NSUtils.transformHHrespToBinaries(vs,thresh);
            psth = NSUtils.constructPsthFromBinaries(bins);
            sts = NSUtils.convBinariesToSpikeTimes(bins,times);
        end
        
    % Converts a series of binary spike times (implicit time) to actual 
    % spike times in a vector. 
    function sts = convBinariesToSpikeTimes(binaries,timeSeries)
        sts = timeSeries(binaries==1); % Logical vector indexing
    end

    % Applies a simple thresholding to find spikes
        % Returns a series of 1s and 0s, where each 1 represents a spike
        % at the time of its index. 
        function b = transformHHrespToBinaries(voltageTimeSeries,threshold)
            b = zeros(length(voltageTimeSeries),1);
            belowThresh = 1;
            for i = 1:length(voltageTimeSeries)
                v = voltageTimeSeries(i);
                if (v > threshold) & (belowThresh==1) %#ok<AND2>
                    belowThresh = 0;
                    b(i) = 1;
                elseif (v < threshold) & (belowThresh==0) %#ok<AND2>
                    belowThresh = 1;
                end
            end    
        end    
        
        % Computes the PSTH given by a binary array
        % Currently using fixed bin size. 
        % Based on Bretschneider 2006 description
        % Gives an estimated instantaneous firing rate per millisecond
        %   for every point in the binaries
        function psth = constructPsthFromBinaries(binaries)
            delta = 400; % 10 ms -> binaries are in 0.025ms steps -> 400 points
            ss = 0.025; % Step size in ms (note: points per ms = 10)
            N = ceil(length(binaries)/delta); % Number of bins
            spikesPerBin = zeros(N,1);
            i = 1; currBin = 1;
            while i < length(binaries)
                endNum = i + delta - 1;
                if ( (i + delta) > length(binaries) ) % Last bin may be partial delta
                    endNum = length(binaries);
                end
                numSpikes = sum(binaries(i:endNum));
                spikesPerBin(currBin) = numSpikes;
                i = i + delta; 
                currBin = currBin + 1;
            end
            x = (spikesPerBin / (delta*ss))*1000; % Normalize to spikes per s
            r = repmat(x,1,delta)';
            temp = r(:);
            psth = temp(1:length(binaries));
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        % % % % % % % % GENERAL PLOTTING UTILITIES % % % % % % % % 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Plots across multi-trial output of the 7 current stims:
        %   5sc, 1bc, 2n. 
        function plotMultiTrialOutput(cell,ss,synapses,stsCa,psths,graphNoise)
            figure; % Plot the chirps and the response
            numStims = 5;
            numPlts = 5;
            k = size(stsCa);
            numTrials = k(2); % Get the number of trials from the sts cell array matrix
            tmax = round(length(synapses(1).Esyn.conved)/10); % Assumes 1st one is chirp!          
            tsR = 0.1:0.1:tmax;
            tsS = ss:ss:tmax;
            mn = min(tsR);
            mx = max(tsR);
            ssr = tsR(2)-tsR(1);
            for i = 1:numStims
                % Plot chirp
                chirp = synapses(i).Esyn.getChirp();
                subplot(numPlts,numStims,i);
                plot(chirp); xlim([ssr,mx/ssr]);
                % Plot Convolved synaptic activity (unscaled)
                s1e = synapses(i).Esyn; 
                s1i = synapses(i).Isyn;
                eNo = s1e.conved; %.preCalcInterpVals; %s1e.gSynArray(tsR);
                iNo = s1i.conved; %.preCalcInterpVals; %s1i.gSynArray(tsR);
                subplot(numPlts,numStims,i+numStims);
                plot(eNo); xlim([mn,mx/ssr]);
                subplot(numPlts,numStims,i+2*numStims);
                plot(iNo); xlim([mn,mx/ssr]);
                % Plot TS neuron raster
                subplot(numPlts,numStims,i+3*numStims);
                for j = 1:numTrials
                    sts = stsCa{i,j};
                    hold on;
                    plot(sts,ones(length(sts),1)*j,'*');
                    xlim([mn/10,mx/ssr/10]); % Its in real time
                    ylim([0,numTrials+1]); % Little extra room
                end
                % Plot TS neuron psth
                subplot(numPlts,numStims,i+4*numStims);
                plot(psths{i}); xlim([min(tsS),max(tsS)/(tsS(2)-tsS(1))]);
            end
            if ~graphNoise; return; end % If GN is not requested, end now
            figure; % Plot the noise and the response
            numNoises = 2;
            numPlts = 4;
            counter = 1;
            for i = (numStims+1):length(synapses)
                % Get times
                tmax = round(length(synapses(1).Esyn.conved)/10);
                tsR = 0.1:0.1:tmax;
                tsS = ss:ss:tmax;  
                ss1 = tsR(2)-tsR(1);
                % Plot noise
                subplot(numPlts,numNoises,counter);
                chirp = synapses(i).Esyn.getChirp();
                plot(chirp); xlim([ss1,tmax]);
                % Plot synapses
                s1e = synapses(i).Esyn; 
                s1i = synapses(i).Isyn;
                eNo = s1e.gSynArray(tsR);
                iNo = s1i.gSynArray(tsR);
                subplot(numPlts,numNoises,counter+numNoises);
                plot(eNo); xlim([ss1,tmax]);
                subplot(numPlts,numNoises,counter+2*numNoises);
                plot(iNo); xlim([ss1,tmax]);
                % Plot TS neuron psth
                subplot(numPlts,numNoises,counter+3*numNoises);
                plot(psths{i}); xlim([min(tsS),max(tsS)/(tsS(2)-tsS(1))]);
                counter = counter + 1;
            end    
            %title
%             ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 ...
%                 1],'Box','off','Visible','off','Units','normalized',...
%                 'clipping' , 'off');
            s = sprintf('f_h:%.2f,f_ca:%.2f,I_bias:%.2f,Ws:%.2f,sigmaB:%.2f,N_xi:%.2f',...
                cell.f_h,cell.f_ca,cell.I_bias,cell.Ws,cell.sigmaB,cell.N_xi);
            text(0.5, 1,sprintf('\bf %s',s),...
                'HorizontalAlignment','center','VerticalAlignment', 'top');
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        % % % % % % % % % MATHEMATICAL UTILITIES % % % % % % % % %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Simple Numerical integration scheme for stochastic differential
        %   equations based on the Euler-Maruyama scheme with a Gaussian
        %   white noise term in the voltage.
        % General Refs: Maruyama 1955; Kloeden & Platen 1992
        % Neuroscience Methods Refs: Manwani & Koch 1999
        % Neuroscience Use Refs: McGillivray, 2012; Guo, 2011; 
        function [ts,y] = FastNeuronEulerMaruyama(neuron,ss,tmax) 
            
            % % % % % % % % % CHECK INPUT % % % % % % % % %
                if ~isa(neuron,'TorusNeuronMod')
                    error('Attempted EulerMaruyama without proper cell type.');
                end
                if ~( (ss > 0) && (ss < 1) && (length(ss)==1) )
                    error('Unacceptable Step Size');
                end
                if ~( (length(tmax)==1) && (tmax > 1) )
                    error('Unacceptable max time');
                end    
            % % % % % % % % % PARAMETERS % % % % % % % % %
                % Retrieve parameters from neuron
                gMu = neuron.gaussMean;
                gSigma = neuron.gaussVariance;
                D = neuron.N_xi;   
                y0 = neuron.getInitialVector();
                Cm = neuron.Cm;
                %dvecByDt = @(t,v) neuron.vprime(t,v,neuron);
                % Fixed parameters
                ts = ss:ss:tmax;     % times of evaluation
                
            % % % % % % % % % ALGORITHM % % % % % % % % %
                % Generate and Precalculate variables
                N = length(ts);
                xi = D * normrnd(gMu,gSigma,N,1);
                y = [ y0 zeros(4,N-1) ];
                sqrtdt = sqrt(ss);
                if strcmpi(neuron.alphaSynapses(1).type,'E')
                    asynValsE = neuron.alphaSynapses(1).preCalcInterpVals;
                    asynValsI = neuron.alphaSynapses(2).preCalcInterpVals;
                else
                    asynValsE = neuron.alphaSynapses(2).preCalcInterpVals;
                    asynValsI = neuron.alphaSynapses(1).preCalcInterpVals;                     
                end    
                esAmpa = neuron.alphaSynapses(1).E_syn;
                
                % Run algorithm
                for i = 2:N
                    normNoise = [sqrtdt * (xi(i)/Cm); 0; 0; 0]; 
                    %step = ss * neuron.vprime( ts(i), y(:,i-1), neuron );
                    
                    %%%%%%%%%%%%%% START: Compute Step %%%%%%%%%%%%%%
                    t = ts(i);
                    %vec = y(:,i-1);
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
                    I_Na = neuron.gNa*(minf^3)*(0.85-n)*(v-neuron.ENa);
                    % Delayed rectifier Potassium Current %
                    I_K = neuron.gK*(n^4)*(v-neuron.EK);
                    % HCN Channel-mediated Hyperpolarization Activated Current %
                    I_H = (neuron.g_h*neuron.f_h*h) * (v-neuron.Eh);
                    % Leak Current %
                    I_L = neuron.gL*(v-neuron.EL);
                    % T-type Calcium Current %
                    sinf = 1/(1+exp(-1*(v+69)/7.8));
                    I_Ca = neuron.gT_Ca*(sinf^3)*hca*(v-neuron.ECa)*neuron.f_ca;
                    % Synaptic Input Current %
                    ind = round(t/ss);
                    sb = neuron.sigmaB;
                    I_syn = neuron.Ws * (2*(v-esAmpa)) * ...
                        ( (asynValsE(ind)*(sb))+(asynValsI(ind)*(1-sb)) );
                    % Sum Currents (adds negative sign here)
                    I_total = -I_Na - I_K - I_H - I_L - I_Ca - I_syn + neuron.I_bias;
                    % Get final dv/dt
                    dvdt = I_total/neuron.Cm;
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
                    hCadot = 2*(hinfca - hca)/neuron.tau_hCa;
                    %vp(4) = hCadot;
                    %%%% End Get T' (Ca) %%%%
                    
                    %%%% Combine the vectors %%%%
                    vp = [dvdt; ndot; hdot; hCadot];
                    
                    %%%%%%%%%%%%%% END: Compute Step %%%%%%%%%%%%%%
                    step = ss * vp;
                    y(:,i) = y(:,i-1) + step + normNoise;
                end
        end
        
         % Simple Numerical integration scheme for stochastic differential
        %   equations based on the Euler-Maruyama scheme with a Gaussian
        %   white noise term in the voltage.
        % General Refs: Maruyama 1955; Kloeden & Platen 1992
        % Neuroscience Methods Refs: Manwani & Koch 1999
        % Neuroscience Use Refs: McGillivray, 2012; Guo, 2011; 
        function [ts,y,Is] = FastNeuronEulerMaruyamaRetIs(neuron,ss,tmax) 
            
            % % % % % % % % % CHECK INPUT % % % % % % % % %
                if ~isa(neuron,'TorusNeuronMod')
                    error('Attempted EulerMaruyama without proper cell type.');
                end
                if ~( (ss > 0) && (ss < 1) && (length(ss)==1) )
                    error('Unacceptable Step Size');
                end
                if ~( (length(tmax)==1) && (tmax > 1) )
                    error('Unacceptable max time');
                end    
            % % % % % % % % % PARAMETERS % % % % % % % % %
                % Retrieve parameters from neuron
                gMu = neuron.gaussMean;
                gSigma = neuron.gaussVariance;
                D = neuron.N_xi;   
                y0 = neuron.getInitialVector();
                Cm = neuron.Cm;

                % Fixed parameters
                ts = ss:ss:tmax;     % times of evaluation
                
            % % % % % % % % % ALGORITHM % % % % % % % % %
                % Generate and Precalculate variables
                N = length(ts);
                xi = D * normrnd(gMu,gSigma,N,1);
                y = [ y0 zeros(4,N-1) ];
                sqrtdt = sqrt(ss);
                if strcmpi(neuron.alphaSynapses(1).type,'E')
                    asynValsE = neuron.alphaSynapses(1).preCalcInterpVals;
                    asynValsI = neuron.alphaSynapses(2).preCalcInterpVals;
                else
                    asynValsE = neuron.alphaSynapses(2).preCalcInterpVals;
                    asynValsI = neuron.alphaSynapses(1).preCalcInterpVals;                     
                end    
                esAmpa = neuron.alphaSynapses(1).E_syn; % Nernst potential
                
                % Prepare to store copies of the currents as well
                Is.I_syn = zeros(1,N);
                Is.I_Na = zeros(1,N);
                Is.I_K = zeros(1,N);
                Is.I_H = zeros(1,N);
                Is.I_L = zeros(1,N);
                Is.I_Ca = zeros(1,N);
                
                % Run algorithm
                for i = 2:N
                    normNoise = [sqrtdt * (xi(i)/Cm); 0; 0; 0]; 
                    
                    %%%%%%%%%%%%%% START: Compute Step %%%%%%%%%%%%%%
                    t = ts(i);
                    v = y(1,i-1); %vec(1);
                    n = y(2,i-1); %vec(2);
                    h = y(3,i-1); %vec(3);
                    hca = y(4,i-1); %vec(4);

                    %%%% Get dv/dt %%%%
                    % Depolarizing Sodium Current %
                    amh = 0.1*(v+40.7);
                    am = (amh)/(1-exp(-amh));
                    bm = 4*exp(-0.05*(v+49.7));
                    minf =  am/(am + bm);
                    I_Na = neuron.gNa*(minf^3)*(0.85-n)*(v-neuron.ENa);
                    % Delayed rectifier Potassium Current %
                    I_K = neuron.gK*(n^4)*(v-neuron.EK);
                    % HCN Channel-mediated Hyperpolarization Activated Current %
                    I_H = (neuron.g_h*neuron.f_h*h) * (v-neuron.Eh);
                    % Leak Current %
                    I_L = neuron.gL*(v-neuron.EL);
                    % T-type Calcium Current %
                    sinf = 1/(1+exp(-1*(v+69)/7.8));
                    I_Ca = neuron.gT_Ca*(sinf^3)*hca*(v-neuron.ECa)*neuron.f_ca;
                    % Synaptic Input Current %
                    ind = round(t/ss);
                    sb = neuron.sigmaB;
                    I_syn = neuron.Ws * (2*(v-esAmpa))*((asynValsE(ind)*(sb))+(asynValsI(ind)*(1-sb)));
                    % Sum Currents (adds negative sign here)
                    I_total = -I_Na - I_K - I_H - I_L - I_Ca - I_syn + neuron.I_bias;
                    % Get final dv/dt
                    dvdt = I_total/neuron.Cm;
                    %%%% End get dv/dt %%%%

                    %%%% Get N' %%%%
                    an = (0.01*(v+40.7))/(1-exp(-0.1*(v+40.7)));
                    bn = 0.125*exp(-0.0125*(v+50.7));
                    ninf = an/(an+bn);
                    taun = 0.05/(an+bn);
                    ndot = (ninf-n)/taun;
                    %%%% End Get N' %%%%

                    %%%% Get h' %%%%
                    tauh = exp(0.033*(v+75))/( 0.011*( 1 + exp(0.083*(v+75)) ) );
                    v_hm = -73; % half maximal voltage in time constant of h
                    hinf = 1/(1+exp(0.151*(v-v_hm)));
                    hdot = (hinf-h)/tauh;
                    %%%% End Get h' %%%%

                    %%%% Get T' (Ca) %%%%
                    q = sqrt(0.25+exp( (v+82)/6.3 ) );
                    hinfca = 1/(0.5+q);
                    hCadot = 2*(hinfca - hca)/neuron.tau_hCa;
                    %%%% End Get T' (Ca) %%%%
                    
                    %%%% Store the currents %%%%
                    Is.I_syn(1,i) = I_syn;
                    Is.I_Na(1,i) = I_Na;
                    Is.I_K(1,i) = I_K;
                    Is.I_H(1,i) = I_H;
                    Is.I_L(1,i) = I_L;
                    Is.I_Ca(1,i) = I_Ca;
                    
                    %%%% Combine the vectors %%%%
                    vp = [dvdt; ndot; hdot; hCadot];
                    
                    %%%%%%%%%%%%%% END: Compute Step %%%%%%%%%%%%%%
                    
                    step = ss * vp;
                    y(:,i) = y(:,i-1) + step + normNoise;
                end
        end
        
        % Simple Numerical integration scheme for stochastic differential
        %   equations based on the Euler-Maruyama scheme with a Gaussian
        %   white noise term in the voltage.
        % General Refs: Maruyama 1955; Kloeden & Platen 1992
        % Neuroscience Methods Refs: Manwani & Koch 1999
        % Neuroscience Use Refs: McGillivray, 2012; Guo, 2011; 
        function [ts,y,gs] = FastNeuronEulerMaruyamaRetGs(neuron,ss,tmax) 
            
            % % % % % % % % % CHECK INPUT % % % % % % % % %
                if ~isa(neuron,'TorusNeuronMod')
                    error('Attempted EulerMaruyama without proper cell type.');
                end
                if ~( (ss > 0) && (ss < 1) && (length(ss)==1) )
                    error('Unacceptable Step Size');
                end
                if ~( (length(tmax)==1) && (tmax > 1) )
                    error('Unacceptable max time');
                end    
            % % % % % % % % % PARAMETERS % % % % % % % % %
                % Retrieve parameters from neuron
                gMu = neuron.gaussMean;
                gSigma = neuron.gaussVariance;
                D = neuron.N_xi;   
                y0 = neuron.getInitialVector();
                Cm = neuron.Cm;

                % Fixed parameters
                ts = ss:ss:tmax;     % times of evaluation
                
            % % % % % % % % % ALGORITHM % % % % % % % % %
                % Generate and Precalculate variables
                N = length(ts);
                xi = D * normrnd(gMu,gSigma,N,1);
                y = [ y0 zeros(4,N-1) ];
                sqrtdt = sqrt(ss);
                if strcmpi(neuron.alphaSynapses(1).type,'E')
                    asynValsE = neuron.alphaSynapses(1).preCalcInterpVals;
                    asynValsI = neuron.alphaSynapses(2).preCalcInterpVals;
                else
                    asynValsE = neuron.alphaSynapses(2).preCalcInterpVals;
                    asynValsI = neuron.alphaSynapses(1).preCalcInterpVals;                     
                end    
                esAmpa = neuron.alphaSynapses(1).E_syn; % Nernst potential
                
                % Prepare to store copies of the currents as well
                gs.g_syn = zeros(1,N);
                gs.g_Na = zeros(1,N);
                gs.g_K = zeros(1,N);
                gs.g_H = zeros(1,N);
                gs.g_L = zeros(1,N);
                gs.g_Ca = zeros(1,N);
                
                % Run algorithm
                for i = 2:N
                    normNoise = [sqrtdt * (xi(i)/Cm); 0; 0; 0]; 
                    
                    %%%%%%%%%%%%%% START: Compute Step %%%%%%%%%%%%%%
                    t = ts(i);
                    v = y(1,i-1); %vec(1);
                    n = y(2,i-1); %vec(2);
                    h = y(3,i-1); %vec(3);
                    hca = y(4,i-1); %vec(4);

                    %%%% Get dv/dt %%%%
                    % Depolarizing Sodium Current %
                    amh = 0.1*(v+40.7);
                    am = (amh)/(1-exp(-amh));
                    bm = 4*exp(-0.05*(v+49.7));
                    minf =  am/(am + bm);
                    GNA = neuron.gNa*(minf^3)*(0.85-n);
                    I_Na = GNA*(v-neuron.ENa);
                    % Delayed rectifier Potassium Current %
                    GK = neuron.gK*(n^4);
                    I_K = GK*(v-neuron.EK);
                    % HCN Channel-mediated Hyperpolarization Activated Current %
                    GH = (neuron.g_h*neuron.f_h*h);
                    I_H = GH * (v-neuron.Eh);
                    % Leak Current %
                    I_L = neuron.gL*(v-neuron.EL);
                    % T-type Calcium Current %
                    sinf = 1/(1+exp(-1*(v+69)/7.8));
                    GTCA = neuron.gT_Ca*(sinf^3)*hca*neuron.f_ca;
                    I_Ca = GTCA*(v-neuron.ECa);
                    % Synaptic Input Current %
                    ind = round(t/ss);
                    sb = neuron.sigmaB;
                    GSYN = 2 * neuron.Ws * ( (asynValsE(ind)*(sb)) + (asynValsI(ind)*(1-sb)) );
                    %I_syn = neuron.Ws * (2*(v-esAmpa))*((asynValsE(ind)*(sb))+(asynValsI(ind)*(1-sb)));
                    I_syn = (v-esAmpa) * GSYN;
                    % Sum Currents (adds negative sign here)
                    I_total = -I_Na - I_K - I_H - I_L - I_Ca - I_syn + neuron.I_bias;
                    % Get final dv/dt
                    dvdt = I_total/neuron.Cm;
                    %%%% End get dv/dt %%%%

                    %%%% Get N' %%%%
                    an = (0.01*(v+40.7))/(1-exp(-0.1*(v+40.7)));
                    bn = 0.125*exp(-0.0125*(v+50.7));
                    ninf = an/(an+bn);
                    taun = 0.05/(an+bn);
                    ndot = (ninf-n)/taun;
                    %%%% End Get N' %%%%

                    %%%% Get h' %%%%
                    tauh = exp(0.033*(v+75))/( 0.011*( 1 + exp(0.083*(v+75)) ) );
                    v_hm = -73; % half maximal voltage in time constant of h
                    hinf = 1/(1+exp(0.151*(v-v_hm)));
                    hdot = (hinf-h)/tauh;
                    %%%% End Get h' %%%%

                    %%%% Get T' (Ca) %%%%
                    q = sqrt(0.25+exp( (v+82)/6.3 ) );
                    hinfca = 1/(0.5+q);
                    hCadot = 2*(hinfca - hca)/neuron.tau_hCa;
                    %%%% End Get T' (Ca) %%%%
                    
                    %%%% Store the conductances %%%%
                    gs.g_syn(1,i) = GSYN;
                    gs.g_Na(1,i) = GNA;
                    gs.g_K(1,i) = GK;
                    gs.g_H(1,i) = GH;
                    gs.g_Ca(1,i) = GTCA;
                    
                    %%%% Combine the vectors %%%%
                    vp = [dvdt; ndot; hdot; hCadot];
                    
                    %%%%%%%%%%%%%% END: Compute Step %%%%%%%%%%%%%%
                    
                    step = ss * vp;
                    y(:,i) = y(:,i-1) + step + normNoise;
                end
        end       
        
        % Takes a folder of data structures in any format
        % Returns a column vector of strings of target filenames
        function list = getFilesInDir(folder)
            list = {};
            pat = '.mat';
            fs = dir(folder);
            for i = 1 : length(fs)
                targ = fs(i);
                if targ.isdir==0
                    s = targ.name;
                    %  Determines whether string s ends with pattern pat
                    sl = length(s);
                    pl = length(pat);
                    b = (sl >= pl && strcmp(s(sl-pl+1:sl), pat)) || isempty(pat);
                    % End external code
                    if b==1
                        list{end+1} = [folder,'/',s]; %#ok<AGROW>
                    end
                end
            end
        end
        
        
    end
end