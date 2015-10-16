clear;
data_path='C:\Users\Peter Hu\Desktop\twolevel_stepperdopr853\res';
parameter_path=strcat(data_path,'\relative_parameters.dat');
collection_effect=input('No.1 considers collection effect, No.0 not:\n');

fid=fopen(parameter_path,'r');
omega0=fread(fid,1,'double');
omegaL=fread(fid,1,'double');
Omega0=fread(fid,1,'double');
mu=fread(fid,1,'double');
len=fread(fid,1,'unsigned long');
dt=fread(fid,1,'double'); %a.u.
T=fread(fid,1,'double');
fclose(fid);
wmg=2*pi/dt;

%plot laser field && fft 
formatDouble='%f64';
source_path=strcat(data_path,'\source.txt');
source_fread=fopen(source_path,'r');
source=textscan(source_fread,formatDouble);
source=cell2mat(source);

dipole_path=strcat(data_path,'\dipole.txt');
dipole_fread=fopen(dipole_path,'r');
dipole=textscan(dipole_fread,'%f64 %f64 %f64');
dipole=cell2mat(dipole);
t=dipole(:,1);  %unit in T
dipole=dipole(:,2);
fclose(source_fread);fclose(dipole_fread);

figure
plot(t,source,'b-','linewidth',2);
title('Laser field');
xlabel('Time(In Optical Period )','fontsize',14);
ylabel('\it{E(t)}(a.u.)','fontsize',14);

%fft && plot
if collection_effect==1
    dipole=dipole*7.5E24*(5.29E-11)^3;
end
len2=length(dipole);
figure
plot(t,dipole,'r-','linewidth',2)
title('dipole');
xlabel('Time(In Optical Period)','fontsize',14);
ylabel('Dipole Intensity (a.u.)','fontsize',14);

fre=(0:round(len2/2)-1)/len2*wmg/omegaL;
FFA=fft(dipole);
hff=FFA;
FFA=abs(FFA(1:round(len2/2))/length(dipole));%
FFA=2*log10(abs(FFA)+eps);
figure
plot(fre,FFA,'r-','linewidth',2);
title('HHG');
xlabel('Harmonic Order(\omega/\omega_L)','fontsize',14);
ylabel('Harmonic Intensity(arb.unit)','fontsize',14);

%wavelet transform
%考虑尺度如何采样的问题，尺度平均的话，频率就不平均，反之亦然
%%resample first and freq units transform
index=1:2^2:len2;
tSample=t(index);   %单位，T
tSample=tSample*T*2.418884326505E-17;   %单位，s
dipoleSample=dipole(index);%a.u. 
Fs=2^13/2^2*(3/8*1.0E15);    %采样频率,考虑重采样
fc=centfrq('cmor1-1');   %Hz
freqL=input('Input the Lower limit of freqrange (unit in order): ');
freqU=input('Input the Upper limit of freqrange (unit in order): ');
freqL=freqL*omegaL*(1.0/2.418884326505E-17)/(2*pi);  %a.u. to Hz
freqU=freqU*omegaL*(1.0/2.418884326505E-17)/(2*pi);  %a.u. to Hz
freqrange=[freqL freqU];
scalerange=fc./(freqrange*(1/Fs));
scalesNum=1024; %尺度采样个数
freqs=linspace(freqL,freqU,scalesNum);
scales=fc./(freqs*(1/Fs));
%scales=linspace(scalerange(end),scalerange(1),scalesNum);
Coeffs=cwt(dipoleSample,scales,'cmor1-1');
figure;
SCImg = wscalogram('image',abs(Coeffs),'scales',scales,'ydata',dipoleSample,'xdata',tSample);

%freqs=scal2frq(scales,'cmor1-1',1/Fs)*(2*pi)*2.418884326505E-17/omegaL;
freqs=freqs*(2*pi)*2.418884326505E-17/omegaL;
[Freqs,Tcenter]=meshgrid(freqs,tSample);
figure;
surf(Tcenter,Freqs,log(abs(Coeffs)'));shading('interp');view(0,90);