clear;
close all;
warning('off');
%% parameters
main_magnetic_field = 7; %[T]

n_TR =800;

T1=900e-3; T2 = 40e-3; % [s] : T1 and T2
PD = 1; % Proton Density
TR = 5e-3; TE = 2e-3; % Repetition time, Echo time
flip_angle=5;
inhomogeneity = 0; % [Hz]

phantom_T1 = load('phantom101').phantom_T1;
phantom_T2 = load('phantom101').phantom_T2;
phantom_PD = load('phantom101').phantom_PD;


FOVx_offset = size(phantom_T1,1); FOVy_offset = size(phantom_T1,2); FOVz_offset = 1; % # of spins

% FOVz = 900; % slice thickness: [mm]

%%
R1 = ones(FOVx_offset*FOVy_offset*FOVz_offset,1)*1/T1; % R1 values: FOVz x 1
R2 = ones(FOVx_offset*FOVy_offset*FOVz_offset,1)*1/T2; % R2 values: FOVz x 1
M0 = ones(FOVx_offset*FOVy_offset*FOVz_offset,1)*PD; % Proton Density values: FOVz x 1

R1 = reshape(1./phantom_T1,FOVx_offset*FOVy_offset*FOVz_offset,1); % R1 values: FOVz x 1
R2 = reshape(1./phantom_T2,FOVx_offset*FOVy_offset*FOVz_offset,1); % R2 values: FOVz x 1
M0 = reshape(phantom_PD,FOVx_offset*FOVy_offset*FOVz_offset,1); % Proton Density values: FOVz x 1


gam = 42.57747892e6; % Gyromagnetic ratio, hydrogen


%% RF pulse
figure;
RF_bandwidth = 500; num_lobes = 1;  phase = 0; sampling_rate = 100000;
[rf_pulse_res,t_rf] = rf_pulse_sinc(flip_angle,RF_bandwidth,num_lobes,sampling_rate,phase);
tstep = t_rf(2)-t_rf(1);
TP = [tstep:tstep:TR];

t_size = size(TP,2);

RF_phaselist= phase_shift_angle(10000,117);
RF_phaselist = RF_phaselist/pi*180;

rf_t_size = size(rf_pulse_res,2);
%% Magnetic field

B0 = zeros(FOVx_offset,FOVy_offset,FOVz_offset, 3, t_size); % main magnetic field initialization
B1 = zeros(FOVx_offset*FOVy_offset*FOVz_offset, 3, t_size); % RF pulse initialization
B_Gz = zeros(FOVx_offset,FOVy_offset,FOVz_offset, 3, t_size); % z-dir gradient initialization
B_Gro = zeros(FOVx_offset,FOVy_offset,FOVz_offset, 3, t_size);

B0(:,:,:,3,:)= main_magnetic_field-main_magnetic_field + inhomogeneity/gam*main_magnetic_field; % rotating frame of reference

% z-direction Gradient:
Gz = (2*pi*RF_bandwidth) / ( FOVz_offset * gam);  % [mT/m] : Rfpulse BW = gamma * Gz * slice thickness
Gz = repmat(reshape(Gz*([1:FOVz_offset] - ceil(FOVz_offset/2)),[1 1 FOVz_offset]),[FOVx_offset FOVy_offset 1]);
B_Gz(:,:,:,3,1:rf_t_size) = repmat(Gz,[1 1 1 1 rf_t_size]);
B_Gz(:,:,:,3,rf_t_size: rf_t_size + rf_t_size/2/4) = repmat((-Gz)*4,[1 1 1 1 rf_t_size/2/4+1]);

% readout Gradient: x
dwell_time = 1/50000;
Gro = (2*pi) / (dwell_time * gam * FOVx_offset); % [mT/m]
Gro = repmat(reshape(Gro*([1:FOVx_offset] - ceil(FOVx_offset/2)),[FOVx_offset 1 1]),[1 FOVy_offset FOVz_offset]);
B_Gro(:,:,:,3,TE*sampling_rate+1:TE*sampling_rate+100) = repmat(Gro,[1 1 1 1 100]);


B = B0 + B_Gz+B_Gro;
%%
B = reshape(B,[FOVx_offset*FOVy_offset*FOVz_offset 3 t_size]);

M_dynamic = [];
echo_value= [];
M_input = repmat([0;0;1],[1 FOVx_offset*FOVy_offset*FOVz_offset]);

for tr = 1:n_TR
    
    [rf_pulse_res,t_rf] = rf_pulse_sinc(flip_angle,RF_bandwidth,num_lobes,sampling_rate,RF_phaselist(tr));
    RFx = reshape(real(rf_pulse_res),[1 1 rf_t_size]);
    RFy = reshape(imag(rf_pulse_res),[1 1 rf_t_size]);
    B1(:,1,1:rf_t_size) = repmat(RFx,[FOVx_offset*FOVy_offset*FOVz_offset 1 1]);
    B1(:,2,1:rf_t_size) = repmat(RFy,[FOVx_offset*FOVy_offset*FOVz_offset 1 1]);
    B_st = B+B1;
    % if tr ==150 % delay
    %     for loc_x = 1:FOVz_offset
    %     ode = @(t,M) [-R2(loc_x) 0 0;
    %         0 -R2(loc_x) 0;
    %         0 0 -R1(loc_x)] * M+ [0; 0; M0(loc_x)*R1(loc_x)];
    %     [timeline, M_delay(loc_x,:,:)] = ode23s(ode, [tstep:tstep:10e-3], M_input(:,loc_x)');
    %     end
    %     M_input = [squeeze(M_delay(:,end,1))';squeeze(M_delay(:,end,2))';squeeze(M_delay(:,end,3))'];
    % end
    parfor loc_x = 1:FOVx_offset*FOVy_offset*FOVz_offset


        warning('off');
        B_temp = squeeze(B_st(loc_x,:,:));
        ode = @(t,M) [-R2(loc_x) gam*B_temp(3,floor(t/tstep)) -gam*B_temp(2,floor(t/tstep));
            -gam*B_temp(3,floor(t/tstep)) -R2(loc_x) gam*B_temp(1,floor(t/tstep));
            gam*B_temp(2,floor(t/tstep)) -gam*B_temp(1,floor(t/tstep)) -R1(loc_x)] * M+ [0; 0; M0(loc_x)*R1(loc_x)];

        [~, sol] = ode23s(ode, TP, M_input(:,loc_x)');
        M_dynamic(loc_x, :, :) = sol;
        te_series(loc_x, tr, :) = sol(TE*sampling_rate, :);

    end

    display([num2str(tr),' TR: ',num2str(toc)])
    tic;
    M_input = [squeeze(M_dynamic(:,end,1))';squeeze(M_dynamic(:,end,2))';squeeze(M_dynamic(:,end,3))'];
end

te_series = reshape(te_series,[FOVx_offset FOVy_offset FOVz_offset size(te_series,2) 3]);
Echo_Mxy = abs(te_series(:,:,:,:,1) + 1i*te_series(:,:,:,:,2));
Echo_Mz = te_series(:,:,:,:,3);
save()
% Echo_Mxy = squeeze(abs(te_series(:,:,1) + 1i*te_series(:,:,2)));
% Echo_Mz = squeeze(te_series(:,:,3));
% figure; hold on;
% plot(squeeze(mean(Echo_Mxy,2)));
% plot(squeeze(mean(Echo_Mz,2)));
% xlabel('# of TR');

figure; hold on;
plot(squeeze(nanmean(Echo_Mxy,[1 2 3])));
plot(squeeze(nanmean(Echo_Mz,[1 2 3])));
xlabel('# of TR');

figure; hold on;
plot(squeeze(mean(Echo_Mxy(1,:,:,200:280),[1 2 3])));
plot(squeeze(mean(Echo_Mxy(11,:,:,200:280),[1 2 3])));
plot(squeeze(mean(Echo_Mxy(5,:,:,200:280),[1 2 3])));
xlabel('# of TR');


figure;
plot(Echo_Mxy(1,:));

% figure; hold on;
% plot(Echo_Mxy(149:end));
%%

figure; 
imagesc(squeeze(mean(Echo_Mxy,4)));

%%
signal = squeeze(mean(Echo_Mxy(:,:,:,:),[2 3]));
for i = 1: size(signal,1)

s = signal(i,:);
N = size(s,2);
ts = 5e-3;
f = [-N/2:N/2-1]/ts/N;
a = fftshift(fft(s));
phase_5hz(i) = angle(a(157));
amp_5hz(i) = abs(a(157));

end
figure; plot(unwrap(phase_5hz));
figure; plot(amp_5hz);


%%
figure;
imagesc(Echo_Mxy(:,:,1,100));
figure; hold on;
plot(squeeze(mean(Echo_Mxy(50-10,:,:,201:241),[2 3]))');
plot(squeeze(mean(Echo_Mxy(22,:,:,201:241),[2 3]))');


a = squeeze(mean(Echo_Mxy(:,:,:,201:240),[2 3]));
a = smoothdata(a,1,"gaussian",5);
[~,I] = max(a,[],2);
figure;plot(-50:50,(I-10)*5);
figure;plot((I-10)*5);

%%
