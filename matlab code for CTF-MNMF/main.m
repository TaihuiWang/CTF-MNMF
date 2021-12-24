clc;
clear all;
addpath('./bss_eval')

%% parameter of CTF-MNMF
type = 1;
refMic = 1;
nb=2;
it=300;
seed=1;
windowSize=128;
shiftCof=0.25;

%% STFT window length and shift size
fsResample = 16000;
fftSize = fsResample*windowSize/1000;   
shiftSize = shiftCof*fftSize; 

%% read mixtures and images
N=2;
NameMixture = ['data/1mixture.wav'];
mix = audioread(NameMixture);
NameImage = ['data/1image.wav'];
ImageAll = audioread(NameImage);
M=size(mix,2);
%% perform separation
RandStream.setGlobalStream(RandStream('mt19937ar','Seed',seed));
fprintf('CTF-MNMF1\n');
[sep, Q] = bss_CTFMNMF1(mix, N, nb, fftSize, shiftSize, it, refMic);

%% evaluation
[SDR0,SIR0,SAR0,perm0]=bss_eval_sources( [mix(:,refMic)';mix(:,refMic)'],[ImageAll(:,1)'; ImageAll(:,2)']);
[SDR,SIR,SAR,perm]=bss_eval_sources(sep',[ImageAll(:,1)';ImageAll(:,2)']);

outNameEst = ['output/sep.wav'];
audiowrite(outNameEst,[sep(:,perm(1)) sep(:,perm(2))],fsResample);