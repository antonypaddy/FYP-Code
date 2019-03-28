% Script to add noise to a set of mitdb records

signals = {'100','101','102','103','104',...
    '105','106','107','108',...
    '109','111','112','113',...
    '114','115','116','117',...
    '118','119','121','122',...
    '123','124','200','201',...
    '202','203','205','207',...
    '208','209','210','212',...
    '213','214','215','217',...
    '219','220','221','222',...
    '223','228','230','231',...
    '232','233','234'};

% types of noise to add
noise(:,1:2) = rdsamp('nstdb/bw');
noise(:,3:4) = rdsamp('nstdb/ma');
noise(:,5:6) = rdsamp('nstdb/em');


% extensions to add to names of records based osnr and noise type 
for type = 1:2:5
    if type == 1
        add = {'e06-em','e00-em','e_6-em','e_12-em','e_18-em','e_24-em','e_30-em','e_36-em','e_42-em','e_48-em','e_54-em'};
    elseif type == 3
        add = {'e06-ma','e00-ma','e_6-ma','e_12-ma','e_18-ma','e_24-ma','e_30-ma','e_36-ma','e_42-ma','e_48-ma','e_54-ma'};
    elseif type == 5
        add = {'e06-bw','e00-bw','e_6-bw','e_12-bw','e_18-bw','e_24-bw','e_30-bw','e_36-bw','e_42-bw','e_48-bw','e_54-bw'};
    end
    % SNR needed 
    SNR = fliplr([-54:6:6]);
    
    % Add noise with function
    for snr = 1:11
        for counting = 1:48
            signal_name = strcat('mitdb/',signals{counting});
            record_name = strcat(signals{counting},add{snr});
            [signalorig,Fs]=rdsamp(signal_name);
            ann = rdann(signal_name,'atr');
            nst(signalorig,noise(:,type:type+1),ann,Fs,SNR(snr),record_name);
        end
    end
end
