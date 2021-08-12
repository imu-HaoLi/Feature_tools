function cepDpDD = log_mel(Sequence)

samples = Sequence;
sr = 16000;
dorasta = 1;
modelorder = 12;

pspectrum = powspec(samples, sr, 0.02, 0.01); %winsize 320, stepsize 160

% next group to critical bands
aspectrum = audspec(pspectrum, sr , 40 , 'mel');
nbands = size(aspectrum,1);

% put in log domain
cep = log10(aspectrum);


    cepDpDD = cep;
end
