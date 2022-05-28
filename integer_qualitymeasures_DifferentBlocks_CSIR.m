
close all;
clear all; 

array1=[]
array1b=[]
array2=[]
array2b=[]
array3=[]
array3b=[]
array4=[]
array4b=[]
array5=[]
array5b=[]

%x=1536;              %why was this block chosen, explain (show other blocks)...recoring 4/28/2022
%y=1599;

% Load and process Spanish dataset of a ship at sea
load('DAP_2010-10-14_09-43-33_010_zayaan_inbound_singlebin_P455_G1_sb_HRR.mat');
sb_HRR.G1.HRR_NoMC_calib = transpose(sb_HRR.G1.HRR_NoMC_calib);

for x=1:64:4409
y=x+64;

HRR_profiles = sb_HRR.G1.HRR_NoMC_calib(x:y,:);
%HRR_profiles = abs(HRR_profiles)/max(max(abs(HRR_profiles)));   %normalised
%ProfileRepetitionFreq = Dataset.ProfileRepitionFreq;
Range_axis = sb_HRR.G1.xaxis_downrange_m; 

HRR_profiles_normalised = abs(HRR_profiles)/max(max(abs(HRR_profiles)));   %normalised

%contrast function
sum_unalign=sum(abs(HRR_profiles_normalised))/64;
Contrast_unalign = sum(sum_unalign.^2,2)/96;  %when plotting the values , plot in db (10*log10{arg})
%entropy
Entropy_unalign =-(sum(sum_unalign.*log(sum_unalign)))/96;


%%
%%% Align HRR profiles using envelope correlation (Single reference profile) %%%
Aligned_Profiles = [ ];
ShiftsRequiredArray = [ ];

for i=1:(length(HRR_profiles(:,1)))

    [a b] = xcorr( abs(HRR_profiles(1,:)), abs(HRR_profiles(i,:)) ); %choose reference profile e.g. 3 TEST TEST TEST
    [maxval maxidx] = max(a);
    ShiftReq = b(maxidx);
    Shifted_M = circshift(HRR_profiles(i,:),ShiftReq);
    Aligned_Profiles = [Aligned_Profiles; Shifted_M];
    ShiftsRequiredArray = [ShiftsRequiredArray; i , ShiftReq];  %for shift plot
end

Aligned_Profiles_normalised = abs(Aligned_Profiles)/max(max(abs(Aligned_Profiles)));   %normalised

%contrast and entropy function
sum_peak=sum(abs(Aligned_Profiles_normalised))/64;
Contrast_peak = sum(sum_peak.^2,2)/96;
%entropy
Entropy_peak =-(sum(sum_peak.*log(sum_peak)))/96;

%%
%%% Align HRR profiles using envelope correlation (Adjacent reference profile) %%%
Aligned_Profiles2 = [abs(HRR_profiles(1,:))];  %start array with first profile
ShiftsRequiredArray = [ ];

for i=2:(length(HRR_profiles(:,1)))
    
    [a b] = xcorr( abs(Aligned_Profiles2(i-1,:)), abs(HRR_profiles(i,:)) ); %each profile get xcorr with the previous updated profile
    [maxval maxidx] = max(a);
    ShiftReq = b(maxidx);
    Shifted_M = circshift(HRR_profiles(i,:),ShiftReq);
    Aligned_Profiles2 = [Aligned_Profiles2; Shifted_M];
    ShiftsRequiredArray = [ShiftsRequiredArray; i , ShiftReq];  %for shift plot
end

Aligned_Profiles_normalised2 = abs(Aligned_Profiles2)/max(max(abs(Aligned_Profiles2)));   %normalised


%contrast and entropy function
sum_adj=sum(abs(Aligned_Profiles_normalised2))/64;
Contrast_adj = sum(sum_adj.^2,2)/96;
%entropy
Entropy_adj =-(sum(sum_adj.*log(sum_adj)))/96;

%%
%%% Align HRR profiles using envelope correlation (continous/cumalative/non_coherent avg reference profile) %%%

Aligned_Profiles3 = [abs(HRR_profiles(1,:))];
ShiftsRequiredArray = [ ];
   
[a b] = xcorr( abs(HRR_profiles(1,:)), abs(HRR_profiles(2,:)) );   %align first 2 profiles because mean function in matlab broken
[maxval maxidx] = max(a);
ShiftReq = b(maxidx);
Shifted_M = circshift(HRR_profiles(2,:),ShiftReq);
Aligned_Profiles3 = [Aligned_Profiles3; Shifted_M];
ShiftsRequiredArray = [ShiftsRequiredArray; 2 , ShiftReq];  %for shift plot

for i=3:(length(HRR_profiles(:,1)))
     M= mean(Aligned_Profiles3(1:i-1,:));
    [a b] = xcorr( abs(M), abs(HRR_profiles(i,:)) ); %each profile get xcorr with the previous updated profile
    [maxval maxidx] = max(a);
    ShiftReq = b(maxidx);
    Shifted_M = circshift(HRR_profiles(i,:),ShiftReq);
    Aligned_Profiles3 = [Aligned_Profiles3; Shifted_M];
    ShiftsRequiredArray = [ShiftsRequiredArray; i , ShiftReq];  %for shift plot
end

Aligned_Profiles_normalised3 = abs(Aligned_Profiles3)/max(max(abs(Aligned_Profiles3)));   %normalised

%contrast and entropy function
sum_avg=sum(abs(Aligned_Profiles_normalised3))/64;
Contrast_avg = sum(sum_avg.^2,2)/96;
%entropy
Entropy_avg =-(sum(sum_avg.*log(sum_avg)))/96;

%%
%%% Align HRR profiles using envelope correlation (Sliding Window avg reference profile :moving avg) %%%

n= 22; %window_length  or 22 or 54-64
%Aligned_Profiles4 = [abs(HRR_profiles(1:n,:))];  before using non-coherentavg

%FIRST ALIGN THE WINDOW ITSELF USING NON CHERENT AVG
Aligned_Profiles4 = [abs(HRR_profiles(1,:))];
ShiftsRequiredArray = [1,0 ];
   
[a b] = xcorr( abs(HRR_profiles(1,:)), abs(HRR_profiles(2,:)) );   %align first 2 profiles because mean function in matlab broken
[maxval maxidx] = max(a);
ShiftReq = b(maxidx);
Shifted_M = circshift(HRR_profiles(2,:),ShiftReq);
Aligned_Profiles4 = [Aligned_Profiles4; Shifted_M];
ShiftsRequiredArray = [ShiftsRequiredArray; 2 , ShiftReq];  %for shift plot


for i=3:n
    M= mean(Aligned_Profiles4(1:i-1,:));
    [a b] = xcorr( abs(M), abs(HRR_profiles(i,:)) ); 
    [maxval maxidx] = max(a);
    ShiftReq = b(maxidx);
    Shifted_M = circshift(HRR_profiles(i,:),ShiftReq);
    Aligned_Profiles4 = [Aligned_Profiles4; Shifted_M];
    ShiftsRequiredArray = [ShiftsRequiredArray; i , ShiftReq];  %for shift plot
end


for i=n+1:(length(HRR_profiles(:,1)))
    M2= mean(Aligned_Profiles4(i-n:i-1,:));    %moving avg
    [a b] = xcorr( abs(M2), abs(HRR_profiles(i,:)) ); %each profile get xcorr with the previous updated profile
    [maxval maxidx] = max(a);
    ShiftReq = b(maxidx);
    Shifted_M = circshift(HRR_profiles(i,:),ShiftReq);
    Aligned_Profiles4 = [Aligned_Profiles4; Shifted_M];
    ShiftsRequiredArray = [ShiftsRequiredArray; i , ShiftReq];  %for shift plot
end

Aligned_Profiles_normalised4 = abs(Aligned_Profiles4)/max(max(abs(Aligned_Profiles4)));   %normalised

%contrast and entropy function
sum_sw_avg=sum(abs(Aligned_Profiles_normalised4))/64;
Contrast_sw_avg = sum(sum_sw_avg.^2,2)/96;
%entropy
Entropy_sw_avg =-(sum(sum_sw_avg.*log(sum_sw_avg)))/96;
%%
% the adj profiles dont work because adj profiles too similar


% Generate Range-Doppler map


% Plot ISAR image profiles

%%
%{
%%%Testing to see best window
%%% Align HRR profiles using envelope correlation (Sliding Window avg reference profile :moving avg) %%%
array_con=[];
array_ent= [];
array_ent_ni=[];
for n=1:64           %chnage window length

%n= 70; %window_length
%Aligned_Profiles4 = [abs(HRR_profiles(1:n,:))];  before using non-coherentavg

%FIRST ALIGN THE WINDOW ITSELF USING NON CHERENT AVG
Aligned_Profiles4 = [abs(HRR_profiles(1,:))];
ShiftsRequiredArray = [1,0 ];
   
[a b] = xcorr( abs(HRR_profiles(1,:)), abs(HRR_profiles(2,:)) );   %align first 2 profiles because mean function in matlab broken
[maxval maxidx] = max(a);
ShiftReq = b(maxidx);
Shifted_M = circshift(HRR_profiles(2,:),ShiftReq);
Aligned_Profiles4 = [Aligned_Profiles4; Shifted_M];
ShiftsRequiredArray = [ShiftsRequiredArray; 2 , ShiftReq];  %for shift plot


for i=3:n
    M= mean(Aligned_Profiles4(1:i-1,:));
    [a b] = xcorr( abs(M), abs(HRR_profiles(i,:)) ); 
    [maxval maxidx] = max(a);
    ShiftReq = b(maxidx);
    Shifted_M = circshift(HRR_profiles(i,:),ShiftReq);
    Aligned_Profiles4 = [Aligned_Profiles4; Shifted_M];
    ShiftsRequiredArray = [ShiftsRequiredArray; i , ShiftReq];  %for shift plot
end



for i=n+1:(length(HRR_profiles(:,1)))
    M2= mean(Aligned_Profiles4(i-n:i-1,:));    %moving avg
    [a b] = xcorr( abs(M2), abs(HRR_profiles(i,:)) ); %each profile get xcorr with the previous updated profile
    [maxval maxidx] = max(a);
    ShiftReq = b(maxidx);
    Shifted_M = circshift(HRR_profiles(i,:),ShiftReq);
    Aligned_Profiles4 = [Aligned_Profiles4; Shifted_M];
    ShiftsRequiredArray = [ShiftsRequiredArray; i , ShiftReq];  %for shift plot
end

Aligned_Profiles_normalised4 = abs(Aligned_Profiles4)/max(max(abs(Aligned_Profiles4)));   %normalised
%{
% Plot Aligned HRR profiles
NumProfiles_axis = 1:size(Aligned_Profiles4,1); 
figure; imagesc(Range_axis, NumProfiles_axis, 20*log10(abs(Aligned_Profiles4))); colorbar
xlabel('Range (m)');
ylabel('Profile number');
title('Aligned range profiles : Envelope corr moving avg ref profile');
axis xy;
colormap('jet');
%}

%{
%PLOTTING THE SHIFTS PER PROFILE
figure;
plot(ShiftsRequiredArray(:,1),ShiftsRequiredArray(:,2))
grid on
xlabel("Profile Number")
ylabel("Shifts required")
title('Number of shifts required per profile: S-W Avg')
%}

%contrast and entropy function
sum_sw_avg=sum(abs(Aligned_Profiles_normalised4));
Contrast_sw_avg = sum(sum_sw_avg.^2,2);
%entropy
Entropy_sw_avg =-(sum(sum_sw_avg.*log(sum_sw_avg)));

array_con= [array_con;n,Contrast_sw_avg];
array_ent= [array_ent;n,1/Entropy_sw_avg];    %invert entropy to get max instead of min
array_ent_ni= [array_ent_ni;n,Entropy_sw_avg]; %non inverted

end

disp(array_con);
array_con2=[array_con(:,2)];
[Max_contr,Index_contr] = max(array_con2) 

disp(array_ent_ni);
array_ent2_ni=[array_ent_ni(:,2)];
[Min_ent_ni,Index_ent_ni] = min(array_ent2_ni)    %because ni(non-inverted) we use min

disp(array_ent);
array_ent2=[array_ent(:,2)];
[Min_ent,Index_ent] = max(array_ent2)    %because we inverted we use max

%array_combined=[array_con2./Max_contr + array_ent2./Min_ent];    %normalised data so none gets more significance
array_combined=[normalize(array_con2,'range')+ normalize(array_ent2,'range')];
[Max_combined,Index_combined] = max(array_combined)

disp(array_ent_ni(Index_combined,2));          %non-inverted,non-normalised
disp(array_con(Index_combined,2));

%}

%%
%db values
Contrast_unalign_db=10*log10(Contrast_unalign);
Entropy_unalign_db =10*log10(Entropy_unalign);
array1=[array1;Contrast_unalign_db];
array1b=[array1b;Entropy_unalign_db];

Contrast_peak_db=10*log10(Contrast_peak);
Entropy_peak_db =10*log10(Entropy_peak);
array2=[array2;Contrast_peak_db];
array2b=[array2b;Entropy_peak_db];

Contrast_adj_db=10*log10(Contrast_adj);
Entropy_adj_db =10*log10(Entropy_adj);
array3=[array3;Contrast_adj_db];
array3b=[array3b;Entropy_adj_db];

Contrast_avg_db=10*log10(Contrast_avg);
Entropy_avg_db =10*log10(Entropy_avg);
array4=[array4;Contrast_avg_db];
array4b=[array4b;Entropy_avg_db];

Contrast_sw_avg_db=10*log10(Contrast_sw_avg);
Entropy_sw_avg_db =10*log10(Entropy_sw_avg);
array5=[array5;Contrast_sw_avg_db];
array5b=[array5b;Entropy_sw_avg_db];

end

figure;
plot(array1./array1)
hold on;
plot(array2./array1)
hold on;
plot(array3./array1)
hold on;
plot(array4./array1)
hold on;
plot(array5./array1)

legend('unaligned', 'single reference' , 'adj reference', 'n-c avg', 's-w avg')

xlabel('Blocks ');
ylabel('Contrast in dB');
title('Contrast per range alignment algorithm for every block (of 64)');

figure;
plot(array1b./array1b)
hold on;
plot(array2b./array1b)
hold on;
plot(array3b./array1b)
hold on;
plot(array4b./array1b)
hold on;
plot(array5b./array1b)

legend('unaligned', 'single reference' , 'adj reference', 'n-c avg', 's-w avg')

xlabel('Blocks ');
ylabel('Entropy in dB');
title('Entropy per range alignment algorithm for every block (of 64)');



