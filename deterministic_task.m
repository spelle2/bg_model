% %% PROBABILISTIC REVERSAL Cools 2006, model as described in baston 2015
clear
close all
clc

%% INITIAL PARAMETER DEFINITION
marker=1;
%define basic parameters 
N_subject=50; % n° subject in the experiment
Dop_tonic=(0.5:0.05:2.1); %value of dopamine for the different group of subjects
Dop_Phasic=[0.7 0.8 0.9 1.0];
gain_drop_dop=0.8;

% we have 2 stimuli and 2 choiche to be taken(red,green)
Ns=2; %n° of neurons in S
Nc=2; %n° of neurons in C

%value of the 2 stimuli [in the case of cools, stimuli represents the highlight of the 2 cards]
S_high=1;
S_small=0;

S1=[S_high; S_small]; % stimulus in which the first card is highlighted
S2=[S_small; S_high]; % stimulus in which the second card is highlighted

% one of the 2 stimuli is reward with a probability of 100,
% the other never rewarded(always punished)

P_high= 1; % prob of reward of the first stimulus
P_small=0; % prob of reward of the other one
P=[P_high P_small];

%define the winning in case of s1 or s2
correct_winner_s1=1;
small_winner_s1=[];

correct_winner_s2=2;
small_winner_s2=[];

% task subdivided in blocks
N_epochs=120; %generally max n° of possible trials in a block

%% different training according to block type
% parameters of the block initialization

% DIFFERENZA SE UN_R O UN_P
switch_types=["unexpected_reward","unexpected_punishment"];
switch_type_in=randi([1,2],1);
switch_type_succession= [switch_types(1),switch_types(2)];

stimoli_iniziali=[S1,S2];
stimoli_iniziali_succ=[stimoli_iniziali(:,1), stimoli_iniziali(:,2)]; 




    

for cont_dop_phasic=1:length(Dop_Phasic)    
    
for k=1:length(Dop_tonic)   %cicle to define the group of subject
    Dop_tonic_value=Dop_tonic(k); 
    Dop_Phasic_value=Dop_Phasic(cont_dop_phasic);
    disp([Dop_tonic_value Dop_Phasic_value])
     % Check if tonic_value is an integer, if so, add '.0' %MS
        if floor(Dop_tonic_value) == Dop_tonic_value
            tonic_str = [num2str(Dop_tonic_value) '.0'];
        else
            tonic_str = num2str(Dop_tonic_value);
        end

        % Check if phasic_value is an integer, if so, add '.0'%MS
        if floor(Dop_Phasic_value) == Dop_Phasic_value
            phasic_str = [num2str(Dop_Phasic_value) '.0'];
        else
            phasic_str = num2str(Dop_Phasic_value);
        end
        for cont_subj=1:N_subject % x each subject
            rng(10+cont_subj)
            switch_type=switch_type_succession(1);   %switch_type_succession(i);
            P=[P_high P_small];
            stimolo_iniziale=stimoli_iniziali_succ(:,1); %stimoli_iniziali_succ(i)
            
            max_trial=120;
            max_trial_corr_rev=14;
            
            cont_stage_rev=0; %n di risposte corrette di seguito
            cont_trial=1; %n di trial
            cont_trial_corr_rev=0; % n di  reversal fatti
            
            name = strcat('W_tot_post_adapt_',num2str(cont_trial));
            noise1=zeros(Nc,N_epochs);
            noise1(1:Nc,1:N_epochs) =  0.08*randn(Nc,N_epochs); %rumore in cortex x subject
            %initial value of synapses before learning
            Wgc = 0.5*diag(ones(Nc,1));   %  weights from cortex to GO
            Wgs = 0.5*(ones(Nc,Ns));       %  weights from stimuli to GO
            Wnc = 0.5*diag(ones(Nc,1));    %  weights from cortex to NOGO
            Wns = 0.5*(ones(Nc,Ns));       %  weights from stimuli to NOGO
            
            %synapses initialization
            Wgc_epocs = zeros(Nc,Nc,N_epochs);
            Wgs_epocs= zeros(Nc,Ns,N_epochs);
            Wnc_epocs = zeros(Nc,Nc,N_epochs);
            Wns_epocs = zeros(Nc,Ns,N_epochs);
            Wgc_epocs(:,:,1) = Wgc;
            Wgs_epocs(:,:,1) = Wgs;
            Wnc_epocs(:,:,1) = Wnc;
            Wns_epocs(:,:,1) = Wns;
            
            %variables for
            vett_correct_answer = zeros(N_epochs,1);
            vett_wrong_asnwer = zeros(N_epochs,1);
            vettore_sw_ok=zeros(N_epochs,2);
            vettore_nsw_ok_rew=zeros(N_epochs,2);
            vettore_nsw_ok_pun=zeros(N_epochs,2);
            cont_ns_rew=zeros(N_epochs,2);
            cont_ns_pun=zeros(N_epochs,2);
            vett_no_response = zeros(N_epochs,1);
            vett_S = zeros(Ns,N_epochs);
            exitC_tot = zeros(Nc,N_epochs);
            idx_reversal=[];
            
            
            max_corr_rev=randi([5,9],1);
            
            %%loop che fra i trial x ogni soggetto
            while ((cont_trial<=120) &&(cont_trial_corr_rev<=max_trial_corr_rev))
                %when the subject perform a certain n of correct answer the
                % probability is switched.
                                
                if cont_stage_rev==max_corr_rev
                    max_corr_rev=randi([5,9],1);
                    P=fliplr(P);
                    cont_trial_corr_rev=cont_trial_corr_rev+1;
                    cont_stage_rev=0;
                    idx_reversal(end+1)=cont_trial;
                    
                    %%NEL LOOP NSUCCESSIVO LO STIMOLO PRESENTATO DEVE ESSERE
                    %%UGUALE --> CONT_TRIAL=IDX_REVERSAL(END)+1
                end
                
                
                Wgc = squeeze(Wgc_epocs(:,:,cont_trial));
                Wgs = squeeze(Wgs_epocs(:,:,cont_trial));
                Wnc = squeeze(Wnc_epocs(:,:,cont_trial));
                Wns = squeeze(Wns_epocs(:,:,cont_trial));
                noiseC = noise1(:,cont_trial);
                num_rand = rand(1,1);
                Correct_winner = [];
                Small_winner = [];
                
                
                %choose stimulus to show and set respectively the right action
                %to perform
                stimulus = 2*randi([1,2],1)-3;
                desired_action=[];
                
                if ~isempty(idx_reversal)
                    if (cont_trial==idx_reversal(end))
                        stimolo_iniziale=flipud(stimolo_iniziale);
                        if (stimolo_iniziale(1,1)==S1(1,1))
                            stimulus=1;
                        else
                            stimulus=-1;
                        end
                    elseif (cont_trial==idx_reversal(end)+1)
                        desired_action=stimulus_before;
                        stimulus=desired_action;
                    end
                end
                
                switch stimulus
                    case 1
                        S = S1;
                        unrew(cont_trial)=1;
                        if num_rand <= P(1,1)
                            Correct_winner = correct_winner_s1;
                        end
                        if num_rand < P(1,2)
                            Correct_winner = correct_winner_s2;
                            Small_winner = small_winner_s1;
                        end
                    case -1
                        S = S2;
                        unrew(cont_trial)=0;
                        if num_rand <= P(1,1)
                            Correct_winner = correct_winner_s2;
                        end
                        if num_rand < P(1,2)
                            Correct_winner = correct_winner_s1;
                            Small_winner = small_winner_s2;
                        end
                end
                %salvo lo stimolo presentato(utile per trial dopo switch)
                stimulus_before=stimulus;
                %i valori di S devono essere mantenuti fra 0 e 1
                S(find(S>1))=1;
                S(find(S<0))=0;
                S_MOSTRATO(:,cont_trial)=S;
                if isempty(Correct_winner)
                    scelta_desiderata(:,cont_trial)=0;
                else
                    scelta_desiderata(:,cont_trial)=Correct_winner;
                end
                %update delle sinapsi
                [Uc,C,Ugo,Go,IGo_DA_Ach,Unogo,NoGo,INoGo_DA_Ach,Ugpe,Gpe,Ugpi,Gpi,Ut,T,Ustn,STN,E,t,Wgc_post,Wgs_post,Wnc_post,Wns_post,r,k_reward,ChI,sw] = BG_model_function_Ach_det(S,Wgc,Wgs,Wnc,Wns,Correct_winner,Small_winner,Dop_tonic_value,noiseC,gain_drop_dop,scelta_desiderata(:,cont_trial),Dop_Phasic_value);
                
                
                scelta_subj(:,cont_trial)=C(:,end);
                exitC=C(end);  %exit cortex ch
                exitGo=Go(end);
                exitNoGo=NoGo(end);
                
                %% DEVO DISTINGUERE NELLE CORRECT ANSWER SE REWARD O PUNISHMENT
                
                
                if r==1 && (sw==1 || sw==2)
                    vett_correct_answer(cont_trial) = 1;
                    cont_stage_rev=cont_stage_rev+1;
                    if ~isempty(idx_reversal)&&(cont_trial==idx_reversal(end)+1)
                        
                        switch switch_type
                            case "unexpected_reward"
                                vettore_sw_ok(cont_trial,1) = 1;
                            case "unexpected_punishment"
                                vettore_sw_ok(cont_trial,2) = 1;
                        end
                    elseif isempty(idx_reversal) || ( ~(cont_trial==idx_reversal(end)) && ~(cont_trial==idx_reversal(end)+1))
                        
                        switch switch_type
                            case "unexpected_reward"
                                if scelta_desiderata(cont_trial)==1
                                    vettore_nsw_ok_rew(cont_trial,1) = 1;
                                    cont_ns_rew(cont_trial,1) = 1;
                                elseif scelta_desiderata(cont_trial)==2
                                    vettore_nsw_ok_pun(cont_trial,1) = 1;
                                    cont_ns_pun(cont_trial,1) = 1;
                                end
                            case "unexpected_punishment"
                                if scelta_desiderata(cont_trial)==1
                                    vettore_nsw_ok_rew(cont_trial,2) = 1;
                                    cont_ns_rew(cont_trial,2) = 1;
                                elseif scelta_desiderata(cont_trial)==2
                                    vettore_nsw_ok_pun(cont_trial,2) = 1;
                                    cont_ns_pun(cont_trial,2) = 1;
                                end
                        end
                        
                        
                        
                    end
                elseif r==-1 && sw==0
                    vett_wrong_asnwer(cont_trial) = 1;
                    cont_stage_rev=0;
                    if ~isempty(idx_reversal)&&(cont_trial~=idx_reversal(end))&&(cont_trial~=1)&&(cont_trial~=idx_reversal(end)+1)
                        switch switch_type
                            case "unexpected_reward"
                                if scelta_desiderata(cont_trial)==1
                                    cont_ns_rew(cont_trial,1) = 1;
                                elseif scelta_desiderata(cont_trial)==2
                                    cont_ns_pun(cont_trial,1) = 1;
                                end
                            case "unexpected_punishment"
                                if scelta_desiderata(cont_trial)==1
                                    cont_ns_rew(cont_trial,2) = 1;
                                elseif scelta_desiderata(cont_trial)==2
                                    cont_ns_pun(cont_trial,2) = 1;
                                end
                        end
                    elseif isempty(idx_reversal)&&(cont_trial~=1)
                        switch switch_type
                            case "unexpected_reward"
                                if scelta_desiderata(cont_trial)==1
                                    cont_ns_rew(cont_trial,1) = 1;
                                elseif scelta_desiderata(cont_trial)==2
                                    cont_ns_pun(cont_trial,1) = 1;
                                end
                            case "unexpected_punishment"
                                if scelta_desiderata(cont_trial)==1
                                    cont_ns_rew(cont_trial,2) = 1;
                                elseif scelta_desiderata(cont_trial)==2
                                    cont_ns_pun(cont_trial,2) = 1;
                                end
                        end
                        
                    end
                    
                else
                    vett_no_response(cont_trial) = 1;
                    cont_stage_rev=0;
                    if ~isempty(idx_reversal)&&(cont_trial~=idx_reversal(end))&&(cont_trial~=1)&&(cont_trial~=idx_reversal(end)+1)
                        
                        switch switch_type
                            case "unexpected_reward"
                                if scelta_desiderata(cont_trial)==1
                                    cont_ns_rew(cont_trial,1) = 1;
                                elseif scelta_desiderata(cont_trial)==2
                                    cont_ns_pun(cont_trial,1) = 1;
                                end
                            case "unexpected_punishment"
                                if scelta_desiderata(cont_trial)==1
                                    cont_ns_rew(cont_trial,2) = 1;
                                elseif scelta_desiderata(cont_trial)==2
                                    cont_ns_pun(cont_trial,2) = 1;
                                end
                        end
                        
                    elseif isempty(idx_reversal)&&(cont_trial~=1)
                        switch switch_type
                            case "unexpected_reward"
                                if scelta_desiderata(cont_trial)==1
                                    cont_ns_rew(cont_trial,1) = 1;
                                elseif scelta_desiderata(cont_trial)==2
                                    cont_ns_pun(cont_trial,1) = 1;
                                end
                            case "unexpected_punishment"
                                if scelta_desiderata(cont_trial)==1
                                    cont_ns_rew(cont_trial,2) = 1;
                                elseif scelta_desiderata(cont_trial)==2
                                    cont_ns_pun(cont_trial,2) = 1;
                                end
                        end
                    end
                    
                end
                
                vett_S(:,cont_trial) = S';
                
                Wgc_epocs(:,:,cont_trial+1) = Wgc_post;
                Wgs_epocs(:,:,cont_trial+1) = Wgs_post;
                Wnc_epocs(:,:,cont_trial+1) = Wnc_post;
                Wns_epocs(:,:,cont_trial+1) = Wns_post;
                
                r_epocs(:,cont_trial)=r;
                
                exitC_tot (:,cont_trial) = exitC;
                exitGo_tot (:,cont_trial) = exitGo;
                exitNoGo_tot (:,cont_trial) = exitNoGo;
                clear Wgc Wgs Wnc Wns
                cont_trial=cont_trial+1;
            end
            tot_sw_ok(cont_subj)=1-(sum(vettore_sw_ok(:,1))/length(idx_reversal));
            tot_nsw_ok_r(cont_subj)=1-(sum(vettore_nsw_ok_rew(:,1))/sum(cont_ns_rew(:,1)));
            tot_nsw_ok_p(cont_subj)=1-(sum(vettore_nsw_ok_pun(:,1))/sum(cont_ns_pun(:,1)));
            switch switch_type
                case "unexpected_reward"
                    tot_reversal_r(cont_subj) = length(idx_reversal);
                case "unexpected_punishment"
                    tot_reversal_p(cont_subj) = length(idx_reversal);
            end
            
            
            
            
            
            
        end

        err_switch=mean(tot_sw_ok);
        sd_switch  = std(tot_sw_ok);
        err_ns_r=mean(tot_nsw_ok_r);
        err_ns_p=mean(tot_nsw_ok_p);
        n_stage_completed=mean(tot_reversal_r);

        vett_finale=[err_switch;err_ns_p;err_ns_r];
       %save (strcat('results_tonic1',num2str(tonic_str),'_phasic',num2str(phasic_str),'.mat'),"err_switch","err_ns_r","err_ns_p","n_stage_completed","vett_finale")
vettore(:,k,cont_dop_phasic)=vett_finale;
matrix_error(k,cont_dop_phasic)=err_switch;
matrix_standard(k,cont_dop_phasic)=sd_switch;
end
end

save Prova vettore matrix_error matrix_standard Dop_tonic Dop_Phasic



