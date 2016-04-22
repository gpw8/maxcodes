%deconvolve all of the data in the Pacaya_2015_mat folder for the statiosn and events given by stn and dates cells
%welcome to the restart
% march1st_2016_03_02
stn={'PS01','PS02','PS03','PS04','PS05','PS06','PS07','PS08','PS09','PS10','PS12','PS14','PS15','PS16','PS17','PS18','PS19'};
dates={'150110','150111','150112','150113','150114','150115','150116','150117','150118','150119','150120','150121'};
stacompvec={'PS01EHZ','PS01EHN','PS01EHE','PS02EHZ','PS02EHN','PS02EHE','PS03EHZ','PS03EHN','PS03EHE','PS04EHZ','PS04EHN','PS04EHE','PS05EHZ','PS05EHN','PS05EHE','PS06EHZ','PS06EHN','PS06EHE','PS07EHZ','PS07EHN','PS07EHE','PS08EHZ','PS08EHN','PS08EHE','PS09EHZ','PS09EHN','PS09EHE','PS10EHZ','PS10EHN','PS10EHE','PS11EHZ','PS11EHN','PS11EHE','PS12EHZ','PS12EHN','PS12EHE','PS13EHZ','PS13EHN','PS13EHE','PS14EHZ','PS14EHN','PS14EHE','PS15EHZ','PS15EHN','PS15EHE','PS16EHZ','PS16EHN','PS16EHE','PS17EHZ','PS17EHN','PS17EHE','PS18EHZ','PS18EHN','PS18EHE','PS19EHZ','PS19EHN','PS19EHE'};
load('p.mat')
load('s.mat')
load('z.mat')
corstn={'PS11','PS13'};
for luv=1:2
    if luv==1
        e=1:17;
        er=1:12;
    else
        e=1:2;
        er=1:12;
    end
    for i=e
        for j=er
            if luv==1
                filnm=strcat('/local/gpwaite_grp/Pacaya_2015_mat/',char(stn(i)),'/',char(stn(i)),'EHZ20',char(dates(j)),'.mat');
            end
            if luv==2
                filnm=strcat('/local/gpwaite_grp/Pacaya_2015_mat_corr/',char(corstn(i)),'/',char(corstn(i)),'EHZ20',char(dates(j)),'.mat');
            end
            if exist(filnm)>0 % exist function usage found on the MATLAB online help documentation.  http://www.mathworks.com/help/matlab/ref/exist.html?s_tid=gn_loc_drop
                load(filnm);
                compnm='EHZ20';
                if luv==1
                    ldnm=strcat(char(stn(i)),compnm,char(dates(j)));
                end
                if luv==2
                    ldnm=strcat(char(corstn(i)),compnm,char(dates(j)));
                end
                namd=ldnm(1:7);
                nmind=strcmp(stacompvec,namd);  %reference the name with the properly ordered (with respect to the Poles Zeros, and Sensitivity arrays) stacompvec to obtain an the array location which corresponds to that of the correct poles, zeros, and sensitivity
                nmind=find(nmind==1);
                d=wo.data;
                fs=wo.Fs;
                d=d-mean(d); %subtract the data set's mean from the values of the dataset
                win=tukeywin(length(d),.25); %perform a tukey window on the data-the tukey window produces a taper based of an equation utilizing a cosine-tukey win help found on matlab help documetation that comes with the program
                d=d.*win; %multiplying each element in the data set by its tukey correspondent
                nyquist=(fs/2); % finding the nyquist of the data for a later line that avoids aliasing
                lz=length(d); %finding the length of the tukeyed data-same as tukey
                nn=2^nextpow2(lz); %finding the value that would be the next power of 2 greater than the length of the data for the fft (allows fft to operate faster)
                %Load the poles and zeros from the imported .mat files provided by Federica Lanza
                poles=p(:,nmind); % load the poles and construct a matrix with columns corresponding to picks/components
                zeros=z(:,nmind); %load the zeros and construct a matrix with columns corresponding to picks/components
                normalization = 1/abs(polyval(poly(zeros),2*pi*1i)/polyval(poly(poles),2*pi*1i)); % by greg waite, poly creates a polynomial with the poles and zeros as roots of the coeficients and then polyval evaluates it over 2*pi*li, this gives the inverse of the magnitude of the ratio of the zero polynomial to the pole polynomial
                sensitivity=Sensitivity(nmind); %obtain a sensitivity value for the pick/component
                [B,A]=zp2tf(zeros,poles,normalization); %this takes the zeros, poles, and normalization value, and then find a numerator polynomial (B) and a denominator polynomial (A) which are the product of the foil multiplication method of the independent variable and the poles/zeros values
                ff=linspace(0,fs,nn); %changed from nn to nn %creates a vector of freqeuncies from 0 to the sample rate
                ww=ff*2*pi; %transforms the values of the ff frequency vector into angular frequencies
                hh=freqs(B,A,ww)*sensitivity; %this is used to compute the frequency response of the system with the knowledge of how the poles and zeros would affect the frequencies generated for the ww vector                            %check this max check complex
                %deconvolve
                %create a waterlevel to smooth
                waterlevel=1*10^-6; % waterlevel constant, used to fill in zones of low amplitude in the frequency spectrum
                zft=fft(d,nn); %changed from nn to nn %performing the fourier transform allows for the spectrum to be generated so that the amplitude of various frequencies can be observed
                if ~iscolumn(zft), zft=zft'; end %if zft is not arranged in columns then do so
                if ~iscolumn(hh), hh=hh'; end %if the frequency response is not arranged in columns then do so
                tmp=hh.*conj(hh); %multiply each element in the frequency response by its conjugate to remove the complex portions?
                gamma=max(tmp)*waterlevel; %taking the maximum of tmp (once the complex aspect has been removed) and then multiplying it by the waterlevel constant
                newzft=(zft.*conj(hh))./(tmp+gamma); %taking the conjugate of the frequency response and multiplying it by the fourier transform and then dividing it by the maximum of the tmp*waterlevel this gives us the amplitudes
                newz=ifft(newzft,'symmetric');
                deconvolvedd(:,1)=newz(1:lz,1)';
                filenm=char(strcat('dcnvlv',namd,'day',char(dates(j)),'.mat'));
                save(filenm,'deconvolvedd')
            end
            clearvars -except stn dates i corstn j p z Sensitivity stacompvec luv e er
        end
    end
end


load('handel.mat')
sound(y)