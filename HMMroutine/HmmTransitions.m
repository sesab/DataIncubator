%---------------------
% HMM_postfit
%---------------------

% POST FIT PROCESSING

function Sequences=HmmTransitions(Posteriors,Options)
 
% ADMISSIBLE states = states for which p>Detection for at least MinDur consecutive bins
%
% 1) keep only states whose probability pstates>Detection for at least MinDur ms ("admissible states"):
% 2) each admissible state is a point in the gnunits-dimensional
% space on which we run Kmeans
% MinDur=percentage of total number of steps
% 
% hmm_postfit.(area).(sess).(trig{Ev(j)})(ntrials(j)).sequence  
% 1st row: start, end of admissible states in
%           sequence as index of pstates entries
% 2nd row: state number for each event (same as in t.p.m.)
% 3rd row: incremental # of admissible states for Kmeans   
% initialize variables


% % % % % % % TONE_SEQ=0;
% % % % % % % if ~isempty(varargin)
% % % % % % %     TONE_SEQ=1;
% % % % % % %     tone=varargin{1};
% % % % % % % end
% % % % % % % AdjustT=0.2;
% % % % % % % if any(strcmp(fieldnames(HmmParam),'AdjustT'))
% % % % % % %     AdjustT=HmmParam.AdjustT; 
% % % % % % % end
% % % % % % % trig=fieldnames(Spikes); 
% % % % % % % Ev=1:numel(trig);
% % % % % % % BinSize=HmmParam.BinSize;


Detection=Options.Detection;
BinSize=Options.BinSize;
ntrials=numel(Posteriors);
if any(strcmp(fieldnames(Options),'MinDur'))
    MinDur=Options.MinDur;
    MinDur=max([1 MinDur*round(size(Posteriors(1).gamma,2))]);
else
    MinDur=max([1 0.01*round(size(Posteriors(1).gamma,2))]);
end
Sequences=repmat(struct('sequence',[]),1,ntrials);
% keeps track of total time in each state
for trial=1:ntrials
%         win=[HmmParam.windel(1)+AdjustT HmmParam.windel(2)];
% %         win=[HmmParam.windel(1) HmmParam.windel(2)];
%         if TONE_SEQ && ~isempty(regexp(trig{Ev(j)},'self','once'))
%             win(1)=win(1)+tone.(trig{Ev(j)})(trial);
%         end
    % remove first AdjustT s in pstates
    pstates_rate=Posteriors(trial).gamma;
    win=[1, size(pstates_rate,2)];
%         pstates_rate=pstates(:,round(AdjustT/BinSize)+1:end);
    NumStates=size(pstates_rate,1);
    temp_sequence=[]; 
    for st_cnt=1:NumStates
        % 0's where this holds, 1's otherwise
        x=~(pstates_rate(st_cnt,:)>Detection);
        % find all start and end bins and duration of each interval
        dsig = diff([1 x 1]);
        startIndex = find(dsig < 0);
        endIndex = find(dsig > 0)-1;
        duration = endIndex-startIndex+1;
        % keep only events with duration>MinDur
        stringIndex = find(duration >= MinDur);
        % revision: keep 1st states spilling over from
        % pre-delivery, but then ignore their initial
        % transition in HMM3_rates later
%                         % REMOVE 1ST STATES if it's on in the first bin
%                         % (it means it's spilling over from pre-delivery)
%                         if ~isempty(stringIndex)
%                             if startIndex(stringIndex(1))==1
%                                 stringIndex(1)=[];
%                             end
%                         end
        % create sequence data
        % only states lasting more than MinDur ms get included in the sequence 
        tot_time=[]; % collect time interval when admissible states are pstates>Detection;
        if ~isempty(stringIndex)
            feat_temp=zeros(numel(stringIndex),5);
            for ind_st=stringIndex
%                             assigned=0;
                time=[]; temp_startIndex=[]; temp_endIndex=[];
%                             name_cnt=name_cnt+1; % incremental index steps by 1 to include this admissible state
                temp_startIndex=startIndex(ind_st)-1;
                temp_endIndex=endIndex(ind_st)-1;
                time=(temp_endIndex(1)-temp_startIndex(1)+1)*BinSize;
                % collect all start and end times for all states in temp_sequence
                % 1st row=start (s); 
                % 2nd row=end (s); 
                % 3rd row=duration
                % 4th row=state as in pstate
                temp_sequence=[temp_sequence [temp_startIndex*BinSize+win(1);...
                    temp_endIndex*BinSize+win(1); time; st_cnt]];
%                             tot_time=[tot_time; time];

            end
        end
    end
    if ~isempty(temp_sequence)
        [~, I]=sort(temp_sequence(1,:));
        temp_sequence=temp_sequence(:,I);
        % 1st row: start, end of admissible states in sequence
        % 2nd row: state number for each event (same as in t.p.m.)
        % 3rd row: incremental # of admissible states for Kmeans
        Sequences(trial).sequence=temp_sequence;
        %hmm_postfit.(area).(sess).(trig{Ev(j)})(trial).scoreLL=score_LL; % store who's got the max LL 
        % states_data: 
        % field .rates: rows = incremental # of admissible states firing rates
        %               cols: firing rates for each units
%                     states_data.(area).(sess).rates=[states_data.(area).(sess).rates; temp_rates]; 
        % .features field contains: 
        % rows = (incremental # of) admissible states
        % 1st col: state onset
        % 2nd col: state length (s)
        % 3rd col: trig j=1:NumEv to identify which trials.
        % 4th col: 0 if onset of state is before first event (cue for ExpT, taste for UT)
        %          1 if onset of state is in [cue,delivery]; 2 if onset of state is post delivery;
        % 5th col: time of cue onset if state is from ExpT trial; 0 if state is from UT trials; 
%                     states_data.(area).(sess).features=[states_data.(area).(sess).features; feat]; 
    elseif isempty(temp_sequence)
        Sequences(trial).sequence=[];                    
    end


end
