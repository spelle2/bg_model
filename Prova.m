%% Probabilistic Reversal (R.Cools et al. 2001)

clc, clear, close all

%N_training_vett = 50*ones(3,1);  % Considero un numero maggiore di soggetti
noise1_value = 0.15;
Dop_tonic_vett=[1.2];     % Vettore che contiene i valori di dopamina tonica di ogni gruppo di soggetti
Dop_Phasic=[0.8];
m_vettore=[1:0.5:1.9,2:0.1:2.5,3];
N_training_vett = 50*ones(length(Dop_tonic_vett),1);  % Considero un numero maggiore di soggetti

for kk=1:length(Dop_tonic_vett)        % Il loop si ripete per ogni gruppo di soggetti
    for jj=1:length(m_vettore)
      
        
        P_high = 0.8;           % Probabilità di successo per la scelta corretta    standard 0.8
        
        Ns = 1;                 % Numero di neuroni in S
        Nc = 2;                 % Numero di neuroni in C
        % Definisco gli stimoli possibili
        S_high = 1.0;           
        S_small = 0.3;
    
        % Essendo lo stimolo rappresentabile con un solo neurone, in quanto la
        % scelta presentata al soggetto è sempre la stessa, la definisco come
        % uno scalare dal valore costante
        S1 = zeros(Ns,1);
        S1(1) = S_high;
    
        % Definisco le scelte da premiare, decise in base alla probabilità di
        % successo
        Correct_winner_1 = 1;
        Correct_winner_2 = 2;
    
        N_training = N_training_vett(kk);   % Numero di soggetti nel loop seguente
        N_epochs = 40;                      % Numero di iterazioni del loop, ovvero di tentativi di ogni soggetto per apprendere la regola
        Dop_tonic = Dop_tonic_vett(kk);     % Valore di dopamina nel loop seguente
        m_value=m_vettore(jj);


         % Check if tonic_value is an integer, if so, add '.0' %MS
        if floor(Dop_tonic) == Dop_tonic
            tonic_str = [num2str(Dop_tonic) '.0'];
        else
            tonic_str = num2str(Dop_tonic);
        end

        % Check if m_value is an integer, if so, add '.0'%MS
        if floor(m_value) == m_value
            m_str = [num2str(m_value) '.0'];
        else
            m_str = num2str(m_value);
        end

      
 %% Loop che ripete la basal task per ogni soggetto
        for ii = 1:N_training
            name = strcat('W_tot_post_adapt_',num2str(ii));
            rng(10+ii)
            noise1=zeros(Nc,N_epochs);
            noise1(1:Ns,1:N_epochs) =  noise1_value*randn(Ns,N_epochs)*1;    % Rumore in C per ogni soggetto 
            
            % Pesi iniziali delle sinapsi
            Wgc = 0.5*diag(ones(Nc,1));   %  weights from cortex to GO
            Wgs = 0.5*(ones(Nc,Ns));      %  weights from stimuli to GO
            Wnc = 0.5*diag(ones(Nc,1));   %  weights from cortex to NOGO
            Wns = 0.5*(ones(Nc,Ns));      %  weights from stimuli to NOGO  
    
            % Inizializzazione vettori
    
            Wgc_epocs = zeros(Nc,Nc,N_epochs);
            Wgs_epocs= zeros(Nc,Ns,N_epochs);
            Wnc_epocs = zeros(Nc,Nc,N_epochs);
            Wns_epocs = zeros(Nc,Ns,N_epochs);
            Wgc_epocs(:,:,1) = Wgc;
            Wgs_epocs(:,:,1) = Wgs;
            Wnc_epocs(:,:,1) = Wnc;
            Wns_epocs(:,:,1) = Wns;
    
            vett_correct_reward = zeros(N_epochs,1);
            vett_punishment = zeros(N_epochs,1);
            vett_wrong_reward = zeros(N_epochs,1);
            vett_no_response = zeros(N_epochs,1);
            vett_S = zeros(Nc,N_epochs);
            exitC_tot = zeros(Nc,N_epochs);
            exitGo_tot = zeros(Nc,N_epochs);
            exitNoGo_tot = zeros(Nc,N_epochs);
    
            for i = 1:N_epochs      % Loop che rappresenta la basal task
                exitC=zeros(Nc,1);
                Wgc = squeeze(Wgc_epocs(:,:,i));
                Wgs = squeeze(Wgs_epocs(:,:,i));
                Wnc = squeeze(Wnc_epocs(:,:,i));
                Wns = squeeze(Wns_epocs(:,:,i)); close all
    
                noiseC = noise1(:,i);
                S = S1;
    
                % Viene determinata la risposta corretta
                num_rand = rand(1,1);   
                if num_rand < P_high
                    Correct_winner = Correct_winner_1;
                else
                    Correct_winner = Correct_winner_2;
                end
                
                % Viene usato il modello
                Small_winner=[];
                [Uc,C,Ugo,Go,IGo_DA_Ach,Unogo,NoGo,INoGo_DA_Ach,Ugpe,Gpe,Ugpi,Gpi,Ut,T,Ustn,STN,E,t,Wgc_post,Wgs_post,Wnc_post,Wns_post,r,k_reward,ChI,sw] = BG_model_function_Ach(S,Wgc,Wgs,Wnc,Wns,Correct_winner,Small_winner,Dop_tonic,noiseC,Dop_Phasic,m_value);
    
                exitC(:) = C(:,end-1100);
                exitGo=Go(:,end);
                exitNoGo=NoGo(:,end);
    
                if r==1 && sw==1
                    vett_correct_reward(i) = 1;
                elseif r==-1 && sw==0
                    vett_punishment(i) = 1;
                elseif r==1 && sw==2
                    vett_wrong_reward(i) = 1;
                else
                    vett_no_response(i) = 1;
                end
    
                vett_S(:,i) = S;
    
                Wgc_epocs(:,:,i+1) = Wgc_post;
                Wgs_epocs(:,:,i+1) = Wgs_post;
                Wnc_epocs(:,:,i+1) = Wnc_post;
                Wns_epocs(:,:,i+1) = Wns_post;
    
                exitC_tot (:,i) = exitC;
                exitGo_tot (:,i) = exitGo;%1:size(exitGo,2)
                exitNoGo_tot (:,i) = exitNoGo;
                clear Wgc Wgs Wnc Wns
            end
    
            %%
            correct_reward_tot = sum(vett_correct_reward);       %rewards
            punishment_tot = sum(vett_punishment);           %punishments
            no_response_tot = sum(vett_no_response);           %no answers
            wrong_reward_tot = sum(vett_wrong_reward);       %small reward
    
            [~,amax]=max(exitC_tot);
            vett_reward=(amax==1);
    
            %% save the date to the file named W_tot_new. This name can be subsequently changed
            save (name, 'Wgc_epocs', 'Wgs_epocs', 'Wnc_epocs', 'Wns_epocs', 'vett_correct_reward', 'vett_punishment', 'vett_no_response', 'vett_wrong_reward', 'vett_S', 'Dop_tonic','Dop_Phasic','m_value','N_epochs','exitC_tot',"exitGo_tot","exitNoGo_tot","vett_reward")
        end
        clear("C","Go","NoGo","Correct_winner_1","Correct_winner_2","Correct_winner")
    
        % Inverto le ricompense
        Correct_winner_1 = 2;
        Correct_winner_2 = 1;
    
        %% Loop che ripete la reversal task per ogni soggetto
        for ii = 1:N_training
            name = strcat('W_tot_post_adapt_',num2str(ii));
            load(name,'Wgc_epocs', 'Wgs_epocs', 'Wnc_epocs', 'Wns_epocs')
            
            Epoch = 40;     % Definisco l'epoca iniziale del reversal
    
            % Inizializzo vettori e matrici
            Wgc = squeeze(Wgc_epocs(:,:,Epoch+1));
            Wgs = squeeze(Wgs_epocs(:,:,Epoch+1));
            Wnc = squeeze(Wnc_epocs(:,:,Epoch+1));
            Wns = squeeze(Wns_epocs(:,:,Epoch+1));
    
            Wgc_epocs = zeros(Nc,Nc,N_epochs);
            Wgs_epocs= zeros(Nc,Ns,N_epochs);
            Wnc_epocs = zeros(Nc,Nc,N_epochs);
            Wns_epocs = zeros(Nc,Ns,N_epochs);
            Wgc_epocs(:,:,1) = Wgc;
            Wgs_epocs(:,:,1) = Wgs;
            Wnc_epocs(:,:,1) = Wnc;
            Wns_epocs(:,:,1) = Wns;
    
            vett_correct_reward = zeros(N_epochs,1);
            vett_punishment = zeros(N_epochs,1);
            vett_wrong_reward = zeros(N_epochs,1);
            vett_no_response = zeros(N_epochs,1);
            vett_S = zeros(Nc,N_epochs);
    
            % Rumore in C
            rng(10+ii)
            noise1 = zeros(Nc,N_epochs);
            noise1(1:Ns,1:N_epochs) = noise1_value*randn(Ns,N_epochs);
            for i = 1:N_epochs          % Loop che rappresenta la basal task
    
                exitC=zeros(Nc,1);
                Wgc = squeeze(Wgc_epocs(:,:,i));
                Wgs = squeeze(Wgs_epocs(:,:,i));
                Wnc = squeeze(Wnc_epocs(:,:,i));
                Wns = squeeze(Wns_epocs(:,:,i));
                noiseC = noise1(:,i);
               
                % Viene decisa la scelta da ricompensare
                num_rand = rand(1,1);   
                S = S1;
                if num_rand < P_high
                    Correct_winner = Correct_winner_1;
                else
                    Correct_winner = Correct_winner_2;
                end
    
                % Viene eseguito il modello e raccolti i risultati
                Small_winner=[];
                [Uc,C,Ugo,Go,IGo_DA_Ach,Unogo,NoGo,INoGo_DA_Ach,Ugpe,Gpe,Ugpi,Gpi,Ut,T,Ustn,STN,E,t,Wgc_post,Wgs_post,Wnc_post,Wns_post,r,k_reward,ChI,sw] = BG_model_function_Ach(S,Wgc,Wgs,Wnc,Wns,Correct_winner,Small_winner,Dop_tonic,noiseC,Dop_Phasic,m_value);
                exitC(:) = C(:,end-1100);
                exitGo(:)=Go(:,end);
                exitNoGo(:)=NoGo(:,end);
    
                if r==1 && sw==1
                    vett_correct_reward(i) = 1;
                elseif r==-1 && sw==0
                    vett_punishment(i) = 1;
                elseif r==1 && sw==2
                    vett_wrong_reward(i) = 1;
                else
                    vett_no_response(i) = 1;
                end
    
                vett_S(:,i) = S(:);
                Wgc_epocs(:,:,i+1) = Wgc_post;
                Wgs_epocs(:,:,i+1) = Wgs_post;
                Wnc_epocs(:,:,i+1) = Wnc_post;
                Wns_epocs(:,:,i+1) = Wns_post;
            % end
    
            exitC_tot (:,i) = exitC;
            exitGo_tot (:,i) = exitGo;
            exitNoGo_tot (:,i) = exitNoGo;
            clear Wgc Wgs Wnc Wns
            end
    
            reward_tot = sum(vett_correct_reward);            %rewards
            punishment_tot = sum(vett_punishment);            %punishments
            no_answer_tot = sum(vett_no_response);            %no answers
            small_reward_tot = sum(vett_wrong_reward);        %small reward
    
            [~,amax]=max(exitC_tot);
            vett_reward=(amax==2);
    
            % save the date to the file named W_tot_new. This name can be subsequently changed
            %save (strcat('Rev_',name), 'Wgc_epocs', 'Wgs_epocs', 'Wnc_epocs', 'Wns_epocs', 'vett_correct_reward', 'vett_punishment', 'vett_no_response', 'vett_wrong_reward', 'vett_S', 'Dop_tonic','Dop_Phasic','m_value', 'N_epochs','exitC_tot',"exitGo_tot","exitNoGo_tot","vett_reward")
        end
    
        %% Parte del programma che organizza i risultati e grafica i valori di C per epoca
    
        % inizializzo vettori
        Nit=40;
        Correct_winner_M=zeros(Nit,N_training);
        max_cons=zeros(N_training,1);
        exitCtotmu=zeros(size(exitC_tot));
        exitC_M=zeros([size(exitC_tot),N_training]);
    
        v_rew_basal=zeros(Nit,N_training);
    
        for ii=1:N_training     % verifico se i soggetti hanno superato la basal task
            name = strcat('W_tot_post_adapt_',num2str(ii));
            load(name)
            Correct_winner_M(:,ii)=vett_correct_reward;
            [~,amax]=max(exitC_tot);
            I=find(amax==1);
            k=0;
            result=[];
            for i=2:length(I)
                if I(i)-I(i-1) ~=1
                    result = [result; k+1];
                    k=0;
                else
                    k=k+1;
                end
            end
    
            v_rew_basal(:,ii) = vett_reward*100;
            
            result = [result;k+1];
            max_cons(ii)=max(result);
            exitCtotmu=exitCtotmu+exitC_tot;
            exitC_M(:,:,ii) = exitC_tot;
    
        end
        disp(strcat(("Tonic DA: "),num2str(tonic_str)))
        disp(strcat(("m_value DA: "),num2str(m_str)))
        %        disp(strcat(num2str(N_training)," ",pazienti(kk)))
    
        s=size(find(max_cons>=8),1);        % Verifico se i soggetti hanno superato la basal task
        disp("Numero di pazienti che hanno superato la prima fase:")
        disp(s);
    
        figure      % Grafico i valori in uscita di C ad ogni epoca in media
        exitCtotmu=exitCtotmu/N_training;
        exitCtotstd = std(exitC_M,0,3);
        subplot(211)
        errorbar(exitCtotmu(1,:),exitCtotstd(1,:));
        subplot(212)
        errorbar(exitCtotmu(2,:),exitCtotstd(2,:));
    
        % Salvo il numero di soggetti che ha superato il basal task e la % di
        % soggetti che in media hanno passato la prova a una data epoca, e la
        % relativa deviazione standard
        win_first_phase=s;
        mu_v_rew_basal = mean(v_rew_basal,2);
        std_v_rew_basal = std(v_rew_basal,0,2);
    
        % Ripeto i passaggi per la reversal task
    
        Correct_winner_M_Rev=zeros(Nit,N_training);
        max_cons_Rev=zeros(N_training,1);
        exitCtotmu=zeros(size(exitC_tot));
        exitC_M=zeros([size(exitC_tot),N_training]);
        
        v_rew_rev=zeros(Nit,N_training);
    
        for ii=1:N_training     % verifico se i soggetti hanno superato la reversal task
            name = strcat('Rev_W_tot_post_adapt_',num2str(ii));
            load(name)
            Correct_winner_M_Rev(:,ii)=vett_correct_reward;
            [~,amax]=max(exitC_tot);
            I=find(amax==2);
            k=0;
            result=[];
            for i=2:length(I)
                if I(i)-I(i-1) ~=1
                    result = [result; k+1];
                    k=0;
                else
                    k=k+1;
                end
            end
    
            v_rew_rev(:,ii) = vett_reward*100;
    
            result = [result;k+1];
            max_cons_Rev(ii)=max(result);
            exitCtotmu=exitCtotmu+exitC_tot;
            exitC_M(:,:,ii) = exitC_tot;
            
        end
    
        s=size(find(max_cons_Rev>=8),1);            % Verifico se i soggetti hanno superato la reversal task
        win_second_phase = s;
        disp("Numero di pazienti che hanno superato la seconda fase:")
        disp(s);
    
        disp("Numero di pazienti che hanno superato entrambe le fasi:")
        s=size(find(max_cons_Rev>=8 & max_cons>=8),1);
        disp(s);
    
        hold on         % Grafico i valori in uscita di C ad ogni epoca in media
        exitCtotmu=exitCtotmu/N_training;
        exitCtotstd = std(exitC_M,0,3);
        tit=strcat('Evolution of C neural activity by epoch for Tonic DA: ',num2str(tonic_str),',m value: ',num2str(m_str));
        title(tit)
        subplot(211)
        hold on
        errorbar(exitCtotmu(1,:),exitCtotstd(1,:));
        ylim([0 1])
        xlabel('Epoch'), ylabel('Mean C neural activity'),title(strcat("Neuron coding choice 1 for Tonic DA: ",num2str(tonic_str),',m value: ',num2str(m_str)))
        legend('Basal Task','Reversal Task')
        subplot(212)
        hold on
        % plot(exitCtotmu(2,:));
        errorbar(exitCtotmu(2,:),exitCtotstd(2,:));
        ylim([0 1])
        xlabel('Epoch'), ylabel('Mean C neural activity'),title(strcat("Neuron coding choice 2 for Tonic DA: ",num2str(tonic_str),',m value: ',num2str(m_str)))
        % xlabel('Epoch'), ylabel('Mean C neural activity'),title(strcat("Neuron coding choice 2 for tonic dopamine ",num2str(Dop_tonic_vett(kk))))
        % xlabel('Epoch'), ylabel('Mean C neural activity'),title(strcat("Neuron coding choice 2 for tonic dopamine ",num2str((kk))))
        legend('Basal Task','Reversal Task')
       
        % Salvo il numero di soggetti che ha superato il reversal task e la % di
        % soggetti che in media hanno passato la prova a una data epoca, e la
        % relativa deviazione standard, partendo dal valore finale della basal
        % task
        win_task=s;
        mu_v_rew_rev = mean(v_rew_rev,2);
        mu_v_rew_rev = [(100-mu_v_rew_basal(end)); mu_v_rew_rev];
        std_v_rew_rev = std(v_rew_rev,0,2);
        std_v_rew_rev = [std_v_rew_basal(end);std_v_rew_rev];
         

        %save (strcat('Rresults_tonic',num2str(tonic_str),'_m',num2str(m_str),'.mat'),"win_first_phase","win_second_phase","win_task","N_training","mu_v_rew_basal","std_v_rew_basal","mu_v_rew_rev","std_v_rew_rev")
        
    end
end

%%
% % %% Grafico la % di vincitori e la relativa deviazione standard per epoca
% % 
% % N_types = length(Dop_tonic_vett);
% % winner_first_phase=zeros(N_types,1);
% % winner_task=zeros(N_types,1);
% % n_training=zeros(N_types,1);
% % 
% % linestyles = ["r>-.","bo--","k*--","r*--","b*--","ko--","go--","g*--","r>--","b>--","k>--","g>--"];
% % 
% % v_label = strings(N_types,1);
% % x_basal = 1:40;
% % x_rev = 0:40;
% % 
% % for i=1:N_types
% %     
% %     figure
% % 
% %     % carico i risultati salvati in precedenza
% %     name = strcat('winners_dop_tonic_',num2str(Dop_tonic_vett(i)),'.mat');
% %     load(name)
% %     
% %     winner_first_phase(i)=win_first_phase;
% %     winner_task(i)=win_task;
% %     n_training(i)=N_training;
% % 
% %     errorbar(x_basal,mu_v_rew_basal,std_v_rew_basal);
% %     hold on
% %     errorbar(x_rev,mu_v_rew_rev,std_v_rew_rev);
% % 
% %     ylim([0 100])
% %     title(strcat("% of correct answer for tonic dopamine ",num2str(Dop_tonic_vett(i))," by epoch in basal and reversal"))
% %     legend("Basal Task", "Reversal Task")
% %     xlabel('Epoch'), ylabel('% of correct answer')
% %     v_label(i) = strcat("Tonic Dopamine: ", num2str(Dop_tonic_vett(i)));
% % end
% % 
% % %% Grafico la % di soggetti che hanno passato basal e reversal task
% % 
% % figure
% % 
% % x=categorical(["acquisition" "reversal"]);
% % for i = 1:N_types
% % 
% %     hold on;
% % 
% %     plot(x,[winner_first_phase(i)/n_training(i)*100 winner_task(i)/n_training(i)*100],linestyles(i));
% % 
% % end
% % 
% % ylim([0 100]);
% % 
% % hold on
% % xlabel('Stage'), ylabel('Percentage of people passing (%)'), title('Probabilistic Reversal')
% % legend(v_label)
