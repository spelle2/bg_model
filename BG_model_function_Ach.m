%% function used to simulate the network during training                                                                                                                                                                                                                                                               
function [Uc,C,Ugo,Go,IGo_DA_Ach,Unogo,NoGo,INoGo_DA_Ach,Ugpe,Gpe,Ugpi,Gpi,Ut,T,Ustn,STN,E,t,Wgc_post,Wgs_post,Wnc_post,Wns_post,r,k_reward,ChI,sw] = BG_model_function_Ach(S,Wgc,Wgs,Wnc,Wns,Correct_winner,Small_winner,Dop_tonic,noiseC,Dop_Phasic,m_value)
% BG_model_function_Ach -----> returns dynamical behaviour of Basal Ganglia structures and cortex
% S                                             stimulus
% Wgc,Wgs,Wnc,Wns                               returns synaptic weights
% C_CORRECT_WINNER                              sets the desired correct response, if possible, depending on the stimulus
% Uc,Ugo,Unogo,Ugpe,Ugpi,Ut,Ustn                returns the input to the sigmoidal function within time of the corresponding brain structures
% C,Go,NoGo,Gpe,Gpi,T,STN,E                     returns activity within time of the corresponding brain structures and energy in the cortex
% IGo_DA_Ach,INoGo_DA_Ach                       returns the input due to Dopa and Ach to Go and NoGo units
% t                                             returns time
% Wgc_post,Wgs_post,Wnc_post,Wns_post           returns synaptic weights AFTER Hebbian learning
% r                                             returns +1 for reward, -1 for punishment, NaN for no feedback
% r_story                                       return the historical record of rewards and punishments
% k_reward                                      returns position of feedback, NaN for no feedback
% ChI                                           returns activity within time of the cholinegic interneuron
% It is different from the other functions since it use the positive part for all the post-synaptic activities for Wgs Wns Wgc and Wnc
%% tempi
tau = 15;          %basal time constant
tauL = 5*tau;      %time constant of lateral inhibition   %tauS = tau/5;
dt = 0.1;          %step
t = (0:dt:800)';   %time [ms]
D = length(t);     %number of samples  (HO MESSO IO IL -1)

%% initialization of the structure
%C: cortex
Nc = size(Wgc,2);   %neurons in the cortex   Nc = 4
C = zeros(Nc,D);    % activity of neurons
Uc = zeros(Nc,D);   % input to the sigmoidal function
Ul = zeros(Nc,D);   % contribution from the lateral inhibition
%Go: striatum, Go
Ng = size(Wgc,1);
Go = zeros(Ng,D);
Ugo = zeros(Ng,D);
%NoGo: striatum, No-Go
Nn = size(Wns,1);
NoGo = zeros(Nn,D);
Unogo = zeros(Nn,D);
%Gpe: globus pallidus pars externa
Gpe = zeros(Nc,D);
Ugpe = zeros(Nc,D);
%Gpi: globus pallidus pars interna
Gpi = zeros(Nc,D);
Ugpi = zeros(Nc,D);
%T: thalamus
T = zeros(Nc,D);
Ut = zeros(Nc,D);
%STN: sub-thalamic nucleus
STN = zeros(1,D);
Ustn = zeros(1,D);
%E: energy (it is an index of conflict in the cortex)
E = zeros(1,D);
%ChI: cholinergic interneuron
ChI = zeros(1,D);
Uchi = zeros(1,D);
%input DA+ACh to Go and NoGo
IGo_DA_Ach = zeros(Ng,D);
INoGo_DA_Ach = zeros(Nn,D);
%% initialization of synapses
%weights from stimulus to cortex
Ns = size(Wns,2);
Wcs = 1*ones(Nc,Ns); 
%lateral inhibition
L0 = 1.2;  
L = -L0*ones(Nc,Nc)+diag(L0*ones(Nc,1));   %tolto autoanello
%weights from thalamus to  cortex
Wct = 4*diag(ones(Nc,1));  %diagonale
%% the following weights are passed as input parameters since subject to learning 
% weights from No-Go to Gpe (inhibitory)
Wen = -2.2*diag(ones(Nn,1));   %diagonal
% weights from Gpe to Gpi (inhibitory)
Wie = -3*diag(ones(Nc,1));   %diagonal
% weights from Go to Gpi (inhibitory)
Wig = -36*diag(ones(Ng,1));   %diagonal
% weights from cortex to thalamus (excitatory)
Wtc = 3*diag(ones(Nc,1));    %diagonale
% weights from Gpi to thalamus (inibitoriy)
Wti = -3*diag(ones(Nc,1));   %diagonal
%weight from energy to STN (excitatory)
Ke = 7;   
%weight from Gpe to STN (inibitory)
Kgpe = -1;
% weight from STN to Gpe (sinapsi excitatory)
Westn = 1;
% weight from STN to Gpi (excitatory)
Wistn = 30;   
%weight from ChI to Go(inhibitory)
wgchi = -1;
%weight from ChI to NoGo (excitatory)
wnchi = 1;
% gain from DA to Go (excitation)
alpha = 0.75; 
% gain from DA to No-Go (inhibition)
beta = -1;
% gain from DA a Chi (inibition)
gamma = -0.5;
%% computation of the dynamics: low-pass filter  + sigmoid
%parameters of sigmoid
a = 4;
U0 = 1.0;
Igpe = 1.0;   %tonic activity Gpe
Igpi = 3;     %tonic activity Gpi 
Ichi = 1.00;  %tonic activity ChI 
r = NaN;   %r = NaN no rew./initialization, r = 1 reward, r = -1 punishment
%r_dop = NaN;
sw = [];
k_reward = NaN; %time of reward o punishment: k_reward = NaN no reward o initialization
%time pattern of phasic dopamine
latency = 100;  %ms, physiological values 50-110 ms (Schultz 1998) %100
klatency = ceil(latency/dt);
duration = 50;   %ms, physiological values <200 ms (Schultz 1998)
kduration = ceil(duration/dt);

%% initial conditions
Ugpe(:,1) = Igpe;
Ugpi(:,1) = Igpi;
Uchi(1) = Ichi+gamma*Dop_tonic;
gain_ADHD = 1;
m = m_value;

%Dop_tonic1 = 0.8;
for k = 1:D
    %% PRODUCTION OF PHASIC DOPAMINE
    if isnan(r)   %nothing r=reward or punishment or nothing.
        dop_phasic = 0;
        Delta=0;
    else          %k_reward initialization
        if k >= k_reward+klatency && k <= k_reward+klatency+kduration   %time interval of dopamine change
            kk = k_reward+klatency;
            if r >= 0             %reward
                r_dop = (1 - Go(C_winner_list,kk-1)) * 2;   % effect on dopamine positiva (reward) 
                delta_Dop = gain_ADHD * Dop_Phasic * r_dop ^ m;            % formula 2 articolo
            elseif r == -1        %punishment
                
                             gain_drop_dop = 0.5;
                             
                r_dop = Go(C_winner_list,kk-1) * 2;  % effect on dopamine negativa
                delta_Dop = gain_drop_dop * Dop_Phasic * (-r_dop ^ m); % formula 3 articolo
            end
            dop_phasic = delta_Dop;
        else      %time interval without dopamine change
            dop_phasic = 0; %non sono nell'intervallo definito dal secondo if
        end
    end
    
    DA = Dop_tonic + dop_phasic + Delta;
    %% Sigmoidi per ogni sezione
    C(:,k) = 1./(1+exp(-a*(Uc(:,k)-U0)));            %corteccia
    Go(:,k) = 1./(1+exp(-a*(Ugo(:,k)-U0)));          %via del GO
    NoGo(:,k) = 1./(1+exp(-a*(Unogo(:,k)-U0)));      %via del NoGO
    Gpe(:,k) = 1./(1+exp(-a*(Ugpe(:,k)-U0)));        %globo pallido pars externa
    Gpi(:,k) = 1./(1+exp(-a*(Ugpi(:,k)-U0)));        %globo pallido pars interna
    ChI(k) = 1./(1+exp(-a*(Uchi(k)-U0)));            %interneuroni colinergici
    T(:,k) = 1./(1+exp(-a*(Ut(:,k)-U0)));            %talamo 
    STN(k) = 1./(1+exp(-a*(Ustn(k)-U0)));            %nucleo subt.

    %% REWARD, PUNISHMENT OR NOTHING : I look for a winner; in case, I must produce reward or punishment
    if isempty(Correct_winner)
        Correct_winner = 0;
    end
    if ~isnan(Correct_winner)
        if isnan(r)   %I still have no result
            C_winner_list = find(C(:,k) >= 0.9);
            lost = isempty (C_winner_list) ;
            if lost==0   %there is al least one winner
                n_C_winner = length(C_winner_list);
                if n_C_winner == 1   %there is just one winner, I give reward or punishment
                    k_reward = k;   %the time of victory
                    if C_winner_list == Correct_winner
                        r = 1;    %reward
                        sw = 1;
                    elseif C_winner_list == Small_winner
                        r = 1;
                        sw = 2;
                    else
                        r = -1; %punishment
                        sw = 0;
                    end
                end
            end
        end
    end

    %% energy in the cortex (index of a conflict)
    for i = 1:Nc
        for j = i:Nc
            E(k) = E(k)+C(i,k)*C(j,k);
        end
    end
    E(k) = E(k)-(sum(C(:,k).^2));

    %% differential equations with Euler method
    Ul(:,k+1) = Ul(:,k)+dt/tauL*(-Ul(:,k)+L*C(:,k));
    Uc(:,k+1) = Uc(:,k)+dt/tau*(-Uc(:,k)+Wcs*S+Ul(:,k)+Wct*T(:,k)+noiseC);
    IGo_DA_Ach(:,k) = alpha*DA*(Go(:,k)-0.35)+wgchi*ChI(k);
    INoGo_DA_Ach(:,k) = (beta*max(DA,0.4)+wnchi*ChI(k))*ones(Nn,1);

    Ugo(:,k+1) = Ugo(:,k)+dt/tau*(-Ugo(:,k)+Wgs*S+Wgc*C(:,k)+IGo_DA_Ach(:,k));     
    Unogo(:,k+1) = Unogo(:,k)+dt/tau*(-Unogo(:,k)+Wns*S+Wnc*C(:,k)+INoGo_DA_Ach(:,k));
    Ugpe(:,k+1) = Ugpe(:,k)+dt/tau*(-Ugpe(:,k)+Wen*NoGo(:,k)+Westn*STN(k)+Igpe);
    Ugpi(:,k+1) = Ugpi(:,k)+dt/tau*(-Ugpi(:,k)+Wig*Go(:,k)+Wie*Gpe(:,k)+Wistn*STN(k)+Igpi);
    Uchi(k+1) = Uchi(k)+dt/tau*(-Uchi(k)+Ichi+gamma*DA);
    Ut(:,k+1) = Ut(:,k)+dt/tau*(-Ut(:,k)+Wti*Gpi(:,k)+Wtc*C(:,k));
    Ustn(k+1) = Ustn(k)+dt/tau*(-Ustn(k)+Ke*E(k)+Kgpe*sum(Gpe(:,k)));   
    if ~isnan(r)
        if k > k_reward+klatency+kduration      % Se è finito il transitorio della dopamina concludo il loop
            C=C(:,1:k_reward+kduration+klatency);
            Go=Go(:,1:k_reward+kduration+klatency);
            NoGo=NoGo(:,1:k_reward+kduration+klatency);
            break;
        end
    end
end  %(chiudo il for k=1:D)

%%  HEBB RULE
delta_Wgc = zeros(Ng,Nc);
delta_Wgs = zeros(Ng,Ns);
delta_Wnc = zeros(Nn,Nc);
delta_Wns = zeros(Nn,Ns);
if ~isnan(r)   %we had reward o punishment
    if ((k_reward+klatency)>D)==0 &&((k_reward+klatency+kduration)>D)==0   %I do not train the net if still varying
        gain=0.5;
        val_Go=Go(:,(k_reward+klatency+kduration));
        val_NoGo=NoGo(:,(k_reward+klatency+kduration));
        Go_th=0.5;
        NoGo_th=0.5;
        C_th=0.5;
        S_th=0.5;

        delta_Wgc = gain*diag((C(:,k_reward+klatency+kduration)-C_th).*(val_Go-Go_th));
        delta_Wnc = gain*diag((C(:,k_reward+klatency+kduration)-C_th).*(val_NoGo-NoGo_th));
        delta_Wgs = gain*(val_Go-Go_th)*(S-S_th)';
        delta_Wns = gain*(val_NoGo-NoGo_th)*(S-S_th)';
    end
    
end %chiudo l'if
Wgc_post = Wgc+delta_Wgc;
Wgs_post = Wgs+delta_Wgs;
Wnc_post = Wnc+delta_Wnc;
Wns_post = Wns+delta_Wns;
%% CONTROL ON SYNAPSES
%synapse saturation
W_min_go=0.4;
W_min_nogo=0.4;
W_max_go = 0.8;
W_max_nogo =0.8;
Wgc_opt_max = W_max_go*diag(ones(Nc,1));   %diagonal
Wgc_opt_min = W_min_go*ones(Nc,Nc);
Wnc_opt_max = W_max_nogo*diag(ones(Nc,1));   %diagonal
Wnc_opt_min = W_min_nogo*ones(Nc,Nc);
Wgs_opt_max = W_max_go*ones(Nc,Ns);   %full
Wgs_opt_min = W_min_go*ones(Nc,Ns);
Wns_opt_max = W_max_nogo*ones(Nc,Ns);   %full
Wns_opt_min = W_min_nogo*ones(Nc,Ns);


Wgc_post = min(max(Wgc_post,Wgc_opt_min),Wgc_opt_max);
Wgs_post = min(max(Wgs_post,Wgs_opt_min),Wgs_opt_max);
Wnc_post = min(max(Wnc_post,Wnc_opt_min),Wnc_opt_max);
Wns_post = min(max(Wns_post,Wns_opt_min),Wns_opt_max);
end



