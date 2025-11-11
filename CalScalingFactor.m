function factors = CalScalingFactor(channel_path, AllUsers, pg_dB)
    if numel(AllUsers) ~= numel(pg_dB)
        error('Check pg.');
    end
    
    nUE = numel(AllUsers);
    files = dir([channel_path, '*.mat']);
    nTTIs = numel(files);
    
    powers = zeros(nUE, 1);
    for i = 1:nTTIs
        filepath = [channel_path, files(1).name];
        load(filepath, 'H');
        H = H(:,:,:,AllUsers);
        
        powers = powers + 10*log10(squeeze(mean(mean(mean(abs(H(:,:,:,:)).^2)))));
    end
    powers = powers / nTTIs;
    pg_dB = reshape(pg_dB, [nUE, 1]);
    factors = 10 .^ ((pg_dB-powers) / 20);
    factors = reshape(factors, [1,1,1,nUE]);
end

