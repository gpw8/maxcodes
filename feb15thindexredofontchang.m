% feb15thindexredofontchangmarch3rdcorrectedapril5th_2016_04_05
%modifier good this this%%
%changes----on 1/2/2016 were only textual
%on 1/3/2016 should be mostly textual changes
%The purpose of this code is to find the best fit for coda spectral energy to
%obtain coda Q, we will be using a coda lapse time of the square root of 3
%multiplied by the P wave arrival time for an event (to estimate the S wave arrival time) then multiplied by a scalar. 

%------------------check comments to make sure that they reflect what is
%actually being coded
%dont renable search unless looking for suitable coda or start time windows
 tip=1; %start an overall counter 
% %create two for loops to check for a best fit of the start time and coda
% %window 

%on nov 24 the Q plot was disabled and the gridsearch for the lapse time
%and coda window renabled
% renable for search-----------------------------------------
%different events and different stations
eventdate={'150115','150116','150117','150117','150118','150118'};
timeat={'2045','1926','0350','0505','0131','1312'};
stnindx=[1,3,4,10,14,17;3,8,9,12,13,15;4,5,6,8,11,13;2,5,12,14,15,18;3,11,13,14,15,18;2,12,13,14,15,17];
%different events and same stations
% eventdate={'150116','150115','150116','150118','150119','150121'};
% timeat={'1925','2220','0119','2058','1845','1026'};
% stnindx=[2,4,5,6,7,9;2,3,8,10,11,15;1,3,4,5,6,9;1,4,6,7,8,10;1,3,5,11,12,13;2,3,4,6,7,9];
% %original
% eventdate={'150116', '150117','150116','150116','150117','150117','150118','150119','150121','150121'} %create a list of event days
% timeat={'1424','0708','0538','0653','1951','2102','2227','0127','1254','1647'} %create a list of event times
% stnindx=[1,2,4,5,6,8;1,4,5,8,9,12;1,2,3,4,5,8;1,2,4,5,9,11;2,6,12,13,14,17;3,4,7,8,10,12;1,3,4,7,10,12;1,2,4,6,8,11;1,2,3,4,6,7;2,3,4,6,7,8]
st=1;
Sadd=0;
nancount=1;
Cdwndadd=2;
for start=1:6 %start the eventday search
slapse=1;    %start the coda start time multiplier count
for S=1:10 %set window start length multplier; to be multipled by the total S wave arrival time assuming S waves travel at 3^(-1/2) the speed of P waves  
    swindow=1; %start the codawindow count
    for codawindow=3:1:20 %set the coda window length
%         S=4
%         codawindow=10
stnind=stnindx(start,:);%use the start variable to take the indexes of the target stations for the chosen event        
%everythings good just the load and fid path change line 4 5 7 88 these
 eventday=char(eventdate(st));%date of the pick-the format is year/month/day so 150114 would be January 14th, 2015. Get the date from the eventdate array 
 timea=char(timeat(st));%time of the pick hours and minutes-the format is hours/minutes so 1126 would be 11 hours and 26 minutes. Get the time from the timeat array 
smoothingwindow=4; %define the smoothing window for the data Calvert et al. 2013 sets theirs to be around 16 cycles for a moving smoothing window-we do not as the hanning smooth that we use will not accomidate such a large amount of cyles with the lower center frequencies
Snumber=S; %S wave arrival time multiplier for coda window start
% codawindow=10;
%event picking from the file, The following extraction code to index data
%from the cnv was inspired/aided from Dr. Gregory Waites extractpickfiles.m
%code
%change the line below to match the computer data repository location
%---------------------------------------------------------------------------
fid=fopen('/home/campus29/mwguetti/MATLAB/mATLAB/f24','r') %use fopen to open the file/produce a file identifier-special thanks to Federica Lanza for providing the f4 file and picked data.-method can be found from from mathworks help documentaion textscan fopen online and off
 C=textscan(fid,'%s %s %s %s %s %s %s %s %*[^\n]')% scan the text for the data to be placed in a cell array %-this and the next 2 lines methods obtained from the MATLAB R2015a program's internal documentation on textscan
 frewind(fid)
 Cd=textscan(fid,'%s','Delimiter','/n')
 fclose(fid)
 finddate=find(cellfun('isempty',regexp(C{1,1},'1501[12]'))==0);%this should find the dates assuming that nothing else contains the '15011' or '15012' sequences.-special thanks to 'YYC' answering 'bit-questions' question regarding cells in Stackoverflow that presented regexp as an option to find the cell values with the [12] in the string http://stackoverflow.com/questions/8056131/strfind-for-string-array-in-matlab. Also special thanks to the internal MATLAB help for the R2015a program used at school for helping with its detail on regexp and how to use it to find patterns in strings,  %special thanks to Jonas answering N.C.Rolly's question about finding empty cell arrays on stackoverflow http://stackoverflow.com/questions/3400515/how-do-i-detect-empty-cells-in-a-cell-array   
 for ij=1:length(finddate)% start a for loop to correct for incomplete hour/minute times and seconds%-the next 10 lines i created
 if isequal(Cd{1,1}{finddate(ij)}(9),char(32))==1 | isequal(Cd{1,1}{finddate(ij)}(10),char(32))==1 %find if there are blank values between the HHMM value digits-citaion?
     C{1,4}{finddate(ij)}=C{1,5}{finddate(ij)}; %reorder the cell to accommidate the unecessary cell offset that blanks in the HHMM cell would cause, as the program figures that values seperated by blanks to be seperate values and so assigns them extra cells 
     C{1,5}{finddate(ij)}=C{1,6}{finddate(ij)};
     C{1,6}{finddate(ij)}=C{1,7}{finddate(ij)};
     C{1,7}{finddate(ij)}=C{1,8}{finddate(ij)};
     C{1,8}{finddate(ij)}=[];
 end 
     C{1,2}{finddate(ij)}=strrep(Cd{1,1}{finddate(ij)}(8:11),char(32),'0'); %-concept of finding blanks in matlab, does this need citation
  C{1,3}{finddate(ij)}=strrep(Cd{1,1}{finddate(ij)}(13:17),char(32),'0');
 end        
for I=1:8 %create a for loop that will be used to reorder the components of C into a more useful cell array-method from Stackoverflow, or myself with cell array indexing, matlab documentaion on accessing cell array. A method can be found as asked by reve_etrange and answered by gnovice on Stackoverflow http://stackoverflow.com/questions/5349470/matlab-index-a-cell-array-with-cell-array-of-arrays-and-return-a-cell-array 
Aw(I,:)=C{1,I}; %this reorders the rows of the f24 into columns of cells
end
Aw(:,finddate(2:end)-1)=[];%this removes the zero columns (columns with a single zero value at their begining) that acted as buffers at the end of events-method from the internal program documetation for MATLAB R2015a titled Deleting Data from a Cell Array
Aw(:,end)=[];
en=strcmp(Aw(1,:),eventday);    %find all events with that specific day that occurred in the f24 file. Use strcmp to create a logical array with ones in the cells where the proper values is present (where the Aw value matches the eventday value) %concept of indexing used in this and other lines (i.e. using strncmp and find) garnered from Andrey Rubshtein's answer to Benjamin/Dennis Jaheruddin from Stackoverflow http://stackoverflow.com/questions/8061344/how-to-search-for-a-string-in-cell-array-in-matlab
Apt=find(en(1,:)==1);        %Use find to get the indicies.
ep=strcmp(Aw(2,Apt),timea);  %find all events with that specific hour and minute that occurred in the list of events that were sorted by day.
Aqt=find(ep(1,:)==1);        %find the index of the event with the right time from the list of indicies of events in Aw which have the right year,month,day.
Aqt=Apt(Aqt);                %get the index of the Aw column with the right year/month/days string and time string.
eventdaychk=char(Aw(1,Aqt));        %check that the day number for the event is ok-char method can be found in matlab documentaion
yearst=eventdaychk(1:2);      %year of the event 
monthst=eventdaychk(3:4);     %month of the event
dayst=eventdaychk(5:6);       %day of the event
timeachk=char(Aw(2,Aqt)); %find the time of the event
strhr=timeachk(1:2);              %start time hour for event
strmin=timeachk(3:4);             %start time minutes for each event
timeachksec=char(Aw(3,Aqt));            %start time seconds for each event
strtstr=strcat('Start time,','Year:',yearst,', ','Month:',monthst,', ','Day:',dayst,', ','Hour:',strhr,', ','Minutes:',strmin,', ','Seconds:',timeachksec);%strcat can be found in matlab documentation
disp(strtstr);
eventbegin=(eval(strhr)*3600)+(eval(strmin)*60)+eval(timeachksec);       %calculate the total seconds until the event start (origin time) 
eventbeginsample=eventbegin*125;  %define samples for start time---the following text is not necessary? this assumes that all loaded data will be in phase with the same start time for the data recording 
eventlat=char(Aw(4,Aqt));  %find the latitude of the event
if length(eventlat)>8 | length(eventlat)<8 %error prompts documetation can be found in the mathworks help documetation, greg mentioned to use error thresholds to keep lat long values in check
    error('latitude reading scheme malfunctioned')
end
eventlat=eval(eventlat(1:2))+eval(eventlat(4:end))/60; %convert it to decimal degrees while removing the letter N component-method of eval can be found on matlab documetation and can be found on Stackoverflow http://stackoverflow.com/questions/15050437/eval-command-in-matlab this is a site but I am not sure if I saw this to base my code off of?, greg proposed the method of conversion
eventlong=char(Aw(5,Aqt));    %find the longitude of the event 
if length(eventlong)>8 | length(eventlong)<8 
    error('longitude reading scheme malfunctioned')
end
eventlong=-(eval(eventlong(1:2))+eval(eventlong(4:end))/60); %convert it to numerical while removing the letter W component and changing the value to negative since we are in the western hemisphere
if eventlat>14.5 | eventlat<14 %create a check to make sure that the latitude is within reasonable bounds-this and the boundary conditions for the latitude were proposed by Dr. Gregory Waite
    error('Latitude is out of reasonable bounds')
end
if eventlong<-91 | eventlong>-90 %create a check to make sure that the longitude is within reasonable bounds
    error('Longitude is out of reasonable bounds')
end 
%special thanks to Lane community college's web pdf for the refresher on
%DMS do decimal degrees conversion. The method was used to check the
%conversions here.
%http://gis.lanecc.edu/gtft/gtft_readings/gtft_reading_wk2/Working_with_Geographic_Coordinates.pdf.
%Special thanks to google maps for confirming the general location with
%gps values
for Rt=1:12 %length of columns needed to get to the next event in the cnv file          
    if Aqt+Rt<length(Aw);
    gh(Rt)=strncmp(Aw(1,Aqt+Rt),Aw(1,Aqt),3); %compare the first three string characters of the first cell in each column for Rt colunms for similarity to find the next event time
    elseif Aqt+Rt==length(Aw);
    gh(Rt)=1    
    end 
end 
R=find(gh~=0); %find the cells in gh that have an instance of 1
R=R(1,1); %use the first instance the find results to index the next event column 
count=1; %start the count
countt=1;
% sta=1;   %start the station count
%load the poles
load('p.mat')
%load the zeros
load('z.mat')
%load the sensitivity values %choose all of the data that you want to load
%for an event 
load('s.mat')
stacompvec={'PS01EHZ','PS01EHN','PS01EHE','PS02EHZ','PS02EHN','PS02EHE','PS03EHZ','PS03EHN','PS03EHE','PS04EHZ','PS04EHN','PS04EHE','PS05EHZ','PS05EHN','PS05EHE','PS06EHZ','PS06EHN','PS06EHE','PS07EHZ','PS07EHN','PS07EHE','PS08EHZ','PS08EHN','PS08EHE','PS09EHZ','PS09EHN','PS09EHE','PS10EHZ','PS10EHN','PS10EHE','PS11EHZ','PS11EHN','PS11EHE','PS12EHZ','PS12EHN','PS12EHE','PS13EHZ','PS13EHN','PS13EHE','PS14EHZ','PS14EHN','PS14EHE','PS15EHZ','PS15EHN','PS15EHE','PS16EHZ','PS16EHN','PS16EHE','PS17EHZ','PS17EHN','PS17EHE','PS18EHZ','PS18EHN','PS18EHE','PS19EHZ','PS19EHN','PS19EHE'};
for Ld=Aqt+1:Aqt+R-1 %create a for loop that runs from the column after the chosen event name/time column to one column before the next event name/time column in the cnv derived array 
    wt=find(strncmp(Aw(:,Ld),'',2)); %find the index of each empty cell in the analysed column
    qs=char(Aw(1,Ld)); %make sure the first cell in the column contains a string and keep the string as the first name, the travel time for this station is located in the cell below it (to be used later) before another station name which will be used later in the code. This cell should contain only staion names
    qs=qs(1:4);%when using real cnv change from 1:4 to 1:5 %cut the string name to give just the staion and not he uncertianty
    for Ki=2:wt(1)-1 %create a for loop that runs from the 2nd cell of each station name column (the first cell being occupied by the string with only the station name and uncertainty, already accounted for) to the charachter of numeric string 
       qsta=char(Aw(Ki,Ld)); %make sure that that the cell contents are in string form
       if isequal(length(qsta),12)==0 & isequal(length(qsta),4)==0 %isequal is found on matlab documentation 
           error('Travel time scheme malfunctions')
       end
       qst=eval(qsta(1:4)); %choose the time that is added to the event start time by first cutting and then evaluating the previously constructed string
       if qst>9 | qst<0 % make sure that the time jump is not over 3 digits and a decimal place, also that the qst value is not to small either 
           error('Arrival time too large or too small')
       end 
%        staloc(sta)=cellstr(qs); %creating a cell array of the station names used in the event for the purpose of locating them-cellstring can be found on matlab documetation
%        sta=sta+1; %update the station 
       %start day accommidations
%create a series of if statements to modify the eventbeginsample so that
%the start of the recording is accomidated and the data is properly indexd-
%time values provided by Ms. Federica lanza
if qs=='PS01' & eval(dayst)==13
    eventbeginsample=eventbeginsample-((22*3600)+(55*60)+23)*125;
end
if qs=='PS02' & eval(dayst)==10
      eventbeginsample=eventbeginsample-((18*3600)+(49*60)+35)*125;
end
if qs=='PS03' & eval(dayst)==10
       eventbeginsample=eventbeginsample-((21*3600)+(6*60)+35)*125;
end
if qs=='PS04' & eval(dayst)==13
         eventbeginsample=eventbeginsample-((20*3600)+(59*60)+30)*125;
end
if qs=='PS05' & eval(dayst)==13
         eventbeginsample=eventbeginsample-((19*3600)+(43*60)+57)*125;
end
if qs=='PS06' & eval(dayst)==13
         eventbeginsample=eventbeginsample-((15*3600)+(13*60)+52)*125;
end
if qs=='PS07' & eval(dayst)==13
         eventbeginsample=eventbeginsample-((17*3600)+(55*60)+19)*125;
end
if qs=='PS08' & eval(dayst)==15
         eventbeginsample=eventbeginsample-((22*3600)+(2*60)+22)*125;
end
if qs=='PS09' & eval(dayst)==13
         eventbeginsample=eventbeginsample-((20*3600)+(5*60)+32)*125;
end
if qs=='PS10' & eval(dayst)==13
         eventbeginsample=eventbeginsample-((17*3600)+(37*60)+26)*125;
end
if qs=='PS11' & eval(dayst)==14
         eventbeginsample=eventbeginsample-((20*3600)+(44*60)+2)*125;
end
if qs=='PS12' & eval(dayst)==14
         eventbeginsample=eventbeginsample-((21*3600)+(30*60)+4)*125;
end 
if qs=='PS13' & eval(dayst)==14
        eventbeginsample=eventbeginsample-((18*3600)+(54*60)+11)*125;
end
if qs=='PS14' & eval(dayst)==16
    eventbeginsample=eventbeginsample-((18*3600)+(22*60)+10)*125;
end
if qs=='PS15' & eval(dayst)==14
    eventbeginsample=eventbeginsample-((23*3600)+(27*60)+33)*125;
end
if qs=='PS16' & eval(dayst)==15
    eventbeginsample=eventbeginsample-((20*3600)+(16*60)+25)*125;
end
if qs=='PS17' & eval(dayst)==15
    eventbeginsample=eventbeginsample-((19*3600)+(16*60)+48)*125;
end
if qs=='PS18' & eval(dayst)==14
    eventbeginsample=eventbeginsample-((21*3600)+(5*60)+5)*125;
end
if qs=='PS19' & eval(dayst)==14
    eventbeginsample=eventbeginsample-((22*3600)+(24*60)+22)*125;
end
stationstartevent=eventbeginsample+125*qst; %add the station time to the event start time-P wave arrival time
codastart=(((stationstartevent-eventbeginsample)*sqrt(3))*Snumber)/125; %S wave arrival time, assuming that P waves are 3^(1/2) times faster than S waves * the multiplier for the coda window start time (lapse time)then divided by 125 to convert back into seconds-proposed by Dr. Gregory Waite? Calvert et al. 2013 mentions this method.
codastartarch(countt)=codastart; 
%The following deconvolution code was taken and adapted from Dr. Gregory Waites earthquake seismology course 
  % we are only using the vertical components so we will not even bother to
  % deconvolve the horizontal components-the idea of deconvolving first and
  % then filtering was proposed by Dr. Greg Waite
  compnm='EHZ20'; %since we are using only the vertical components this name addition is used for all stations. In this we are following De Sienna et al. 2014 (Attentuation and scattering tomography of the deep plumbing system of Mount Saitn Helens).
    ldnm=strcat(qs,compnm,eventday); %use the current station name +the vertical designation (compnm) and the event dat to create a load name 
    ldnnm=strcat('/local/gpwaite_grp/Pacaya_2015_mat/',qs,'/',ldnm,'.mat'); %create an overall path-loadname for the data-this will differ based on where the data is stored
    if qs=='PS13' | qs=='PS11' %reroute to corrected data provided by Ms. Lanza
        ldnnm=strcat('/local/gpwaite_grp/Pacaya_2015_mat_corr/',qs,'/',ldnm,'.mat')
    end  
    load(ldnnm); %load the data
    namd=strcat(qs,compnm); %make a name with just the currentr station name and the vertical designation
    namd=namd(1:7); %cut the "20" out of the namd name
    nmind=strcmp(stacompvec,namd);  %reference the name with the properly ordered (with respect to the Poles Zeros, and Sensitivity arrays) stacompvec to obtain an the array location which corresponds to that of the correct poles, zeros, and sensitivity
    nmind=find(nmind==1);  %find the index value of the previously obtained array location found in the last line
    data=wo.data; %access the current loaded data structure and rename the data into a new variable, knowledge of how to access data structures from Greg Waite Earthquake seismology course and Mathworks MATLAB help doucmentation
    fs=wo.Fs; %obtain the sample rate from the currently loaded data structure
    namdd(:,countt)=namd; 
    %loop to test all of the components
    %create new variables 
    %--------------------------
%     
%     figure
%     plot(eventbeginsample:eventbeginsample+(codastart+codawindow)*fs,data(eventbeginsample:eventbeginsample+(codastart+codawindow)*fs))
%      saveas(gcf,'tr.jpeg','jpeg') %saveas methods can be found in
%      mathworks matlab help documentation
%      close Figure 1
    %-----------------------
    %deconvolving the data
    d=data(eventbeginsample-1250:eventbeginsample+((codastart+codawindow)*fs)+1250); %rename the data and add 10 seconds (should be) worth of data before and after the event begins that will later be cut to act as a buffer against edge effects-greg
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
%     if iscolumn(zeros); zero=zeros'; end %if zeros are in columns make them in rows
%     if iscolumn(poles); pole=poles'; end  %if poles are in columns make them in rows
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
    clear newzft newz gamma tmp zft hh ww ff B A normalization
    clear d FA FB nn lz win data sensitivity 
 for filt=1:6 %this creates a filter loop which is used to filter the data-deconvolve it and then pick the data according tot he coda window start time (laspe time) and coda window length. It also plots the deconvolved and the picked data as well ast the sqaure of the picked data-calvert?  
     wn1a=[1,2,4,8,16,20]; %lower filter ranges  filters ranges from calvert
     wn2a=[2,4,8,16,32,60]; %upper filter ranges
     centerfreq=[1.5,3,6,12,24,40];
     deconvolved=deconvolvedd;
     wn1=wn1a(filt); %select the lower filter according to the loop progression
     wn2=wn2a(filt); % select the upper filter according to the loop progression 
    [FB,FA]=butter(4,[(wn1/nyquist) (wn2/nyquist)],'bandpass');  % the concept and implimentation of the butterworth filters was derived from Dr. Gregory Waite's filt_traces.m function which was based off a function written by Derek Schutt. -filter help by mathworks matalb documetation- greg informed me to set the order to 4 but not really higher to avoid edge effects
    deconvolved(:,1)=filter(FB,FA,deconvolved); % we use a bandpass filter that has variable frequency ranges as dictated by the for loop following the method by Calvert et al. 2013
    dh=strcat('d',ldnm,'frequencies',sprintf('%d',wn1),'-',sprintf('%d',wn2));
    clear FA FB 
    %get the vertical data for each station component
    %--------------------try to find a source term value of each station
    %per event by finding the maximum between the start of the origin time of the event and the start of the coda window. 
    %Check the standard deviation of these maximums to determine if they are relevant (ie are not affected by source radiation paterns and so are similar)-is so then the average S maximum for an event can be used as the average for the event.-This method was proposed by Dr. Greg Waite check this for many stations below 1 km depth   
       vertdat=deconvolved(codastart*fs+1251:1251+(codastart+codawindow)*fs,1).^2; 
       sourcemax(count)=max(deconvolved(1251:1251+codastart*fs,1).^2);%-created by Dr. Greg Waite
       %also it should be noted that all of the data is derived from the deconvolved data which is taken from the start of the event origin time to the end of the coda window plus a buffer zone after the coda window---------------------------------------------------------------------------------------------------------------
%        figure %create a figure to plot the vertical component-figures proposed by Dr. Gregory Waite
%        subplot(3,1,1)
%        plot(eventbeginsample:eventbeginsample+(codastart+codawindow)*fs,deconvolved(1251:1251+(codastart+codawindow)*fs,1)) %plot the data against the number of samples 
%        xlabel('Time (1/125) seconds') %label the time in 1/fs of a second i.e. samples
%        ylabel('Amplitude')
       titlev=dh; %use the name and the event day to make the title for the first plot which should be the vertical
       titleall{1,count}=titlev; %store the title in a cell array of titles that increase with each count   
       %        title(titlev) %use the title in the vertical plot
%        subplot(3,1,2)
%        plot(eventbeginsample+(codastart*fs):eventbeginsample+(codastart+codawindow)*fs,deconvolved(1251+(codastart*fs):1251+(codastart+codawindow)*fs,1))
%        xlabel('Time (1/125) s') %label the time in 1/fs of a second i.e. samples
%        ylabel('Amplitude')
%        title('Vertdat Prior to being squared')
%        subplot(3,1,3)
%        plot(eventbeginsample+(fs*codastart):eventbeginsample+(fs*(codastart+codawindow)),vertdat) %plot the data against the number of samples
%        xlabel('Time (1/fs) s') %label the time in 1/fs of a second i.e. samples
%        ylabel('Amplitude^2')
%        title('Vertdat')
%        set(gcf,'Position',[200,200,2000,1000]) %set described by Dr. Greg
%        saveas(gcf,dh,'jpeg') %save as with gcf can be found in the matlab help documetaion that comes with the program
       clear titlev dh
       dhc=round(smoothingwindow*fs*1/centerfreq(filt));
       if isequal((-1)^round(smoothingwindow*fs*1/centerfreq(filt)),1^round(smoothingwindow*fs*1/centerfreq(filt)))==0 %we can find out if a number is odd by using it as an exponent of 1 and fininding out if -1 raised to the number is equivalent to +1 raised to the same number, if it is then the number is even.
         dhc=round(smoothingwindow*fs*1/centerfreq(filt))+1; %if the smoothing window is odd then add 1 
       end 
%        close Figure 1 %method seen on matlab help documetaion
%------------------------------should the smooth be applied to to the
%full-no
%traces or just the windowed coda? vertdat=vertdat(fs*codastart:fs*(codastart+codawindow))
       vertsave(1:length(vertdat),count)=vertdat; %save the unsmoothed squared data for use later in the subplots 
       vertdat=hanningsmooth(vertdat,dhc); %uses the function hanning smooth by Dr. Gregory Waite and uses the smoothing window multiplied by the center frequency of the filtered range and 125 sample per second-suggested by Dr. Waite
       pA=vertdat(1:end); 
       pickevent(1:length(pA(:)),count)=pA;%add the picked event data to the pickevent vector %-creating multidimensional vectors method found in matlab help. More can be found on Stackoverflow asked by Theodoros Theodoridis http://stackoverflow.com/questions/23376111/multidimensional-arrays-multiplication-in-matlab 
%        pN=vertdat((floor(xn(1)-(stationstartevent-100))+1):(floor(xn(2)-(stationstartevent-100))+1)); %pick from the vertical data
%        picknoise(1:length(pN(:)),count)=pN;%add the picked noise data to the picknoise vector
       count=count+1; %update the count
       clear deconvolved 
 end
 clear deconvolvedd ldnm
  eventbeginsample=eventbegin*fs;
       clear vertdat wo.data cell
       countt=countt+1;
       if Ki~=(wt(1)-1) %use this if statement to avoid using the last string cell of the column as a station names since it is the travel time of the station at wt-2
       qs=qsta(5:8); %obtain the station name for the next iteration (for the first iteration this replaces the name of the station found in the top cell of the column with the name in the current cell, to be use with the travel time in the next cell lower of the column)
       end 
    end
clear qs %clear qs upon the end of the pick loop
end
count=count-1;


for m=1:6 %-coding lines dealing with the possible S wave maximum (or maximum value between the event origin time and the start of the coda window) were proposed by Dr. Greg Waite.
pickeventa(:,:,m)=pickevent(:,m:6:end); %-creating multidimensional vectors method found in matlab help. More can be found on Stackoverflow asked by Theodoros Theodoridis http://stackoverflow.com/questions/23376111/multidimensional-arrays-multiplication-in-matlab
vertsavea(:,:,m)=vertsave(:,m:6:end);%reorganize the saved unsmoothed vertdat arrays in a similar fashion to the pickeventa
Sourcemax(:,m)=sourcemax(m:6:end);
titlealla(1,:,m)=titleall(1,m:6:end);
stdsourcemax(m)=std(Sourcemax(:,m));
averagesourcemax(m)=mean(Sourcemax(:,m));
end
nstns=length(pickeventa(1,:,1));
Sourcemaxtable=table(titlealla,stdsourcemax,averagesourcemax);
kont1=1;
%%
%grid search method proposed and mostly designed by Dr. Gregory Waite, grid search methods can be found on Stackoverflow
 for freqsza=1:6
     centerfreq=[1.5,3,6,12,24,40];
     kount1=1;
     
     for numstations=1:nstns;
             K=1:length(pickeventa(:,numstations,freqsza));
             K2=(K'+(codastartarch(numstations)*fs)-1)/fs; %start the model with the proper lapse time as measures from the origin time-De. Sienna et al. 2014 (MSH) page 8227-Coding lines (this and next) created by Dr, Greg Waite
             Et=log(pickeventa(:,numstations,freqsza).*K2.^(3/2));%linear regression method from Calvert et al 2013., Wolfram Mathworld helped with the concepts of regression as well as the correlation coeficients others http://mathworld.wolfram.com/CorrelationCoefficient.html? etc., matlab in program help helped with the related code lines.-Greg suggested the use of the linear regression special thanks to user4402918 for his post and kkuilla's answer regarding MATLAB Correlation coefficiants http://stackoverflow.com/questions/28995650/correlation-coefficients-in-matlab. I looked that the MATLAB online documentation as well for reference.  
             polye=polyfit(K2,Et,1);%-matlab in program help, use polyfit to perform a linear regression which will give us the coefficiant of the slope that will be used to find Q
             if polye(1)<0
             Qc=(polye(1).^-1)*(-2*pi*centerfreq(freqsza)); %extract Q from the slope coefficient. 
             linfit=polyval(polye,K2); %evaluate the fit polynomial over the time K2-this and the next three lines were taken from the internal matlab help on linear regression  
             sr=sum((Et-linfit).^2); %sum the square of the residual values between the fit and the data
             stot=(length(Et)-1)*var(Et); %multiply the length of the actual data-1 times the variance of the actual data
             R2=1-sr/stot; %find the coeficient of determination by subtracting the ratio of the sum of the squared residual values over the sum of the differences (squared) between the data and its mean -MATLAB internal help documentaion R2015a "Linear Regression", from 1   
%         frequencyusedstation(kount1)=freqsza;
%         StationQ(kount1)=tzQ(ind);
          Qtocheck(kount1)=Qc;
          pickcheck(:,kount1)=vertsavea(:,numstations,freqsza); %use the pickcheck variable to store the vertsavea data for a particular station and frequency
           %for the particular station and frequency find the NaN values using isnan on pickcheck(:,kount1) and then using find(isnan(pickcheck(:,kount1))==1) to find the indicies of the NaN values and store them in a variable special thanks to Marc for
        %answering Graviton's question on Stackoverflow regarding finding NaN values for it
        %demonstrated how to perform this action with Graviton's method used in
        %this code http://stackoverflow.com/questions/1713724/find-all-nan-elements-inside-an-array
         if sum(isnan(pickcheck(:,kount1)))>0
                     error('NaN')
           end 
          logpick(:,kount1)=Et;
          if sum(isnan(logpick(:,kount1)))>0
                     error('NaN')
           end 
          linfittocheck(:,kount1)=linfit;
          if sum(isnan(linfittocheck(:,kount1)))>0
                     error('NaN')
           end 
%         Station(kount1)=numstations;
%         titleq(kount1)=titlealla(1,numstations,freqsza); 
%         Ecodat(:,kount1)=pickeventa(:,numstations,freqsza);
%         modlpsdq(:,kount1)=modlpsd(:,ind);
          rsend(kount1)=R2;%save the R^2 value in an array which updates with each passing kount1, with the data ultimately being used in the subplots
          titlend(kount1)=titlealla(1,numstations,freqsza); %save the title for the particular station and frequency 
          freqend(kount1)=freqsza; %save the frequency index
          stationend(:,kount1)=namdd(:,numstations); %for a particular station and frequncy save the stationname which will be used later in the subplots, since the namdd variable contains the names of the stations and components in the order in which they wer deconvolved they should be properly indexed by numstations per event 
          windw(kount1)=codawindow; %save the coda window in an updating array that updates with each passing kount1
          laps(kount1)=S; %save the coda laspe time multiplier in an array that updates with each passing kount1
%         figure
        %renable------------------------------------------------------------------------------------
%         hold on
%         plot(codastartarch(numstations)*fs:(codastartarch(numstations)+codawindow)*fs,pickeventa(:,numstations,freqsza)/max(pickeventa(:,numstations,freqsza)),'r')
%         xlabel('Time (1/fs seconds) from event origin')
%         ylabel('Power spectral density (kg meters^2/seconds^2)')
%         tlt=strcat('Minimum Q for station',titlealla(1,numstations,freqsza),'--','rmsq',sprintf('%d',bvsta(kount1)),'--','Q coda:',sprintf('%d',StationQ(kount1)),'--','Center freqency:',sprintf('%d',centerfreq(freqsza)));
%         title(tlt)
%         hbp=char(titlealla(1,numstations,freqsza));
%         hbp=hbp(1,2:8);
%         hqp=char(titlealla(1,numstations,freqsza));
%         hqp=hqp(1,15:16);
%         tltsave=strcat('StnminQ',hbp,'Dy',hqp,'otime',timea,'Q',sprintf('%.3f',tzQ(ind)),'Cf',sprintf('%.1f',centerfreq(freqsza)),'.jpeg');        
%         plot((codastartarch(numstations))*fs:(codastartarch(numstations)+codawindow)*fs,modlpsd(:,ind)/max(modlpsd(:,kount)),'b')
%         hold off
%         set(gcf,'Position',[200,200,2000,1000])
%         saveas(gcf,tltsave,'jpeg')
        kount1=kount1+1; %update the kount1 variable
             else 
                 disp('positive slope for fit, data modified for station iteration')
                 nname=strcat('S: ',sprintf('%d',S),', codawindow: ',sprintf('%d',codawindow),', Freq: ',sprintf('%d',freqsza),' ,stn: ',namdd(:,numstations)');
                 nanproblem{nancount}=nname;
                 nancount=nancount+1;
                 Qtocheck(kount1)=123456789;
                 pickcheck(:,kount1)=vertsavea(:,numstations,freqsza); %use the pickcheck variable to store the vertsavea data for a particular station and frequency
                 if sum(isnan(pickcheck(:,kount1)))>0
                     error('NaN')
                 end 
                 logpick(:,kount1)=repmat(123456789,[length(Et),1]);
                 linfittocheck(:,kount1)=repmat(123456789,[length(K2),1]);
                 rsend(kount1)=123456789;
                 titlend(kount1)=titlealla(1,numstations,freqsza); %save the title for the particular station and frequency 
                 freqend(kount1)=freqsza; %save the frequency index
                 stationend(:,kount1)=namdd(:,numstations); %for a particular station and frequncy save the stationname which will be used later in the subplots, since the namdd variable contains the names of the stations and components in the order in which they wer deconvolved they should be properly indexed by numstations per event 
                 windw(kount1)=codawindow; %save the coda window in an updating array that updates with each passing kount1
                 laps(kount1)=S; %save the coda laspe time multiplier in an array that updates with each passing kount1
             kount1=kount1+1;
             end 
        clear R2 ploye Et Qc linfit sr stot K K2 lonan linan pnan nname% clear the variabels that are no loger necessary after this iteration.
%          close Figure 1
     end
     %this is used to find the overall minimum rmsq for a frequency and will
     %be removed in the final draft of the code
     % for each frequency find the total minimum of Q 
%    bvfreq(freqsza)=min(bvsta);
%    ref(1)=find(bvsta==bvfreq(freqsza));
%    frequencyfinal(freqsza)=centerfreq(freqsza);
%    FrequencyQ(freqsza)=StationQ(ref);
%    Frequncystation(freqsza)=Station(ref);
%    figure
%    renable-------------------------------------------------------------------------------------------------------------
%    hold on
%         plot(1:length(Ecodat(:,ref)),Ecodat(:,ref)/Ecodat(1,ref),'r')
%         xlabel('Time (1/fs seconds) from start time lapse')
%         ylabel('Power spectral density (kg meters^2/seconds^2)')
%         tlt2=strcat('Total minimum Q for frequency',titleq(ref),'--','rmsq:',sprintf('%d',bvfreq(freqsza)),'--','Q coda:',sprintf('%d',FrequencyQ(freqsza)),'--','Center freqency:',sprintf('%d',frequencyfinal(freqsza)));
%         hop=char(titleq(ref));
%         hop=hop(1,1:16);
%         tlt2save=strcat('FreqminQ',hop,'origintime',timea,'Q',FrequencyQ(freqsza),'Centerfreq',sprintf('%d', frequencyfinal(freqsza)),'.jpeg');
%         title(tlt2)
%         plot(1:length(Ecodat(:,ref)),modlpsdq(:,ref),'b')
%         hold off
%         set(gcf,'Position',[200,200,2000,1000])
%         saveas(gcf,tlt2save,'jpeg')
        Qtochecka(:,kont1)=Qtocheck;
        if sum(isnan(Qtocheck))>0
                     error('NaN')
           end 
        linfittochecka(:,:,kont1)=linfittocheck;
        if sum(sum(isnan(linfittocheck)))>0
                     error('NaN')
           end 
        pickchecka(:,:,kont1)=pickcheck;
         if sum(sum(isnan(pickcheck)))>0
                     error('NaN')
           end 
        loggpick(:,:,kont1)=logpick;
         if sum(sum(isnan(logpick)))>0
                     error('NaN')
           end 
        rsenda(:,kont1)=rsend;
        if sum(isnan(rsend))>0
                     error('NaN')
           end 
        titlenda(:,kont1)=titlend;
        if sum(cellfun('isempty',titlend))>0 %guetmade
                     error('No name for title')
           end 
        freqenda(:,kont1)=freqend;
        if sum(isnan(freqend))>0
                     error('NaN')
           end 
        stationenda(:,:,kont1)=stationend;
           if sum(sum(isnan(stationend)))>0
                     error('NaN')
           end
        lapps(:,kont1)=laps;
        if sum(isnan(laps))>0
                     error('NaN')
           end 
        window(:,kont1)=windw;
        if sum(isnan(windw))>0
                     error('NaN')
           end 
clear kount1 Qtocheck rsend titlend freqend stationend laps windw logpick pickcheck linfittocheck
         kont1=kont1+1;
%         close all %remmeber to disable the close alls when plotting
 end
% 
% for hap=1:6 
% if abs(Qtocheck(nstns*(hap-1)+1:nstns*hap)-mean(Qtocheck(nstns*(hap-1)+1:nstns*hap)))<20 & rmsqend(nstns*(hap-1)+1:nstns*hap)<0.2
%     lapse(tip)=S;
%     windowcda(tip)=codawindow;
%     Qfin(tip)=Qtocheck(nstns*(hap-1)+1);
%     rmsqfin(tip)=rmsqend(nstns*(hap-1)+1);
%     titlfin(:,tip)=titlend(nstns*(hap-1)+1:nstns*hap);
%     freqfin(tip)=freqend(nstns*(hap-1)+1);
%     tip=tip+1;
% end %-------------check the and titlefin

Qsc(:,:,swindow)=Qtochecka;
rssc(:,:,swindow)=rsenda;
freqsc(:,:,swindow)=freqenda;
stationsc(:,:,:,swindow)=stationenda;
picksc{swindow}= pickchecka;
logpi{swindow}=loggpick;
linfitsc{swindow}=linfittochecka;
windowsc(:,:,swindow)=window;
lapsesc(:,:,swindow)=lapps;
% daysc(:,swindow)=eventday;
% timesc(swindow)=timea;
swindow=swindow+1;
clearvars -except tip S Sadd Cdwndadd nanproblem nancount slapseroww wccoll wn1a wn2a logpi logpic logpickscan codawindow start swindow st slapse Qsc picksc linfitsc picksca linfitsca pickscan linfitscan Qsca rssca freqsca stationsca lapsesca windowsca daysca timesca rssc freqsc stationsc lapsesc windowsc daysc timesc Qscan rsscan freqscan stationscan lapsescan windowscan dayscan timescan stnindx eventdate eventday timea timeat stnind

    end
    Qsca(:,:,:,slapse)=Qsc;
    rssca(:,:,:,slapse)=rssc;
    freqsca(:,:,:,slapse)=freqsc;
    stationsca(:,:,:,:,slapse)=stationsc;
    picksca{:,slapse}=picksc;
    logpic{:,slapse}=logpi;
    linfitsca{:,slapse}=linfitsc;
    lapsesca(:,:,:,slapse)=lapsesc;
    windowsca(:,:,:,slapse)=windowsc;
%     daysca(:,slapse)=daysc;
%     timesca(:,slapse)=timesc;
 slapse=slapse+1;
 clearvars -except tip S Sadd nanproblem nancount Cdwndadd slapseroww wccoll wn1a logpic logpickscan wn2a start st slapse picksca linfitsca pickscan linfitscan Qsca rssca freqsca stationsca lapsesca windowsca daysca timesca Qscan rsscan freqscan stationscan lapsescan windowscan dayscan timescan stnindx stnind eventdate eventday timea timeat
end

%matlab online documentation helped to use size to troubelshoot.
Qscan{st}=Qsca;%special thanks to Dan for answering potAito's questions regarding creating variables with names from strings in that it inspqired me to use cell arrays for Qscan and the other final scan variables http://stackoverflow.com/questions/16099398/create-variables-with-names-from-strings
rsscan{st}=rssca;
freqscan{st}=freqsca;
stationscan{st}=stationsca;
pickscan{:,:,st}=picksca;
logpickscan{:,:,st}=logpic;
linfitscan{:,:,st}=linfitsca;
lapsescan{st}=lapsesca;
windowscan{st}=windowsca;
     %find the NANs in the Qscan matrix special thanks to Marc for
        %answering Graviton's question on Stackoverflow regarding finding NaN values for it
        %demonstrated how to perform this action with Graviton's method used in
        %this code http://stackoverflow.com/questions/1713724/find-all-nan-elements-inside-an-array
        Qnan=find(isnan(Qscan{st})==1);
        lnan=find(isnan(lapsescan{st})==1);
        wnan=find(isnan(windowscan{st})==1);
        rnan=find(isnan(rsscan{st})==1);
        stnan=find(isnan(stationscan{st})==1);
        fnan=find(isnan(freqscan{st})==1);
        lognan=find(isnan(logpic)==1)
        linftnan=find(isnan(linfitsca)==1)
        pNan=find(isnan(picksca)==1)
        %--------------------------------------------------------------shoudl I add a doubel layer of nan protection? 
        if length(pNan)>0 | length(lognan)>0 | length(linftnan)>0 | length(Qnan)>0 | length(lnan)>0 | length(wnan)>0 | length(rnan)>0 | length(stnan)>0 | length(fnan)>0%if the reference matrices have NaNs then something is wrong and produce error prompt
        error('NaN values in the reference matrices')
        end 
% dayscan(:,:,st)=daysca;
% timescan(:,:,st)=timesca;
%%
for df=1:6
    for gf=1:6
        for wf=1:length(lapsescan{st}(stnind(gf),df,1,:))
            for cf=1:length(windowscan{st}(stnind(gf),df,:,1))
                Qsn(cf)=Qscan{st}(stnind(gf),df,cf,wf);
                Sscann(cf)=lapsescan{st}(stnind(gf),df,cf,wf);
                wndscann(cf)=windowscan{st}(stnind(gf),df,cf,wf);
                stnscann(:,cf)=stationscan{st}(:,stnind(gf),1,cf,wf);
                pickscann{cf}=pickscan{st}{wf}{cf}(:,stnind(gf),df);
                logpickscann{cf}=logpickscan{st}{wf}{cf}(:,stnind(gf),df);
                linfitscann{cf}=linfitscan{st}{wf}{cf}(:,stnind(gf),df);
                frqscann(cf)=freqscan{st}(stnind(gf),df,cf,wf);
                rsscann(cf)=rsscan{st}(stnind(gf),df,cf,wf);
            end
            Q(:,wf)=Qsn;
            Sscan(:,wf)=Sscann;
            wndscan(:,wf)=wndscann;
            stnscan(:,:,wf)=stnscann;
            lpickscan{:,wf}=logpickscann;
            pckscan{:,wf}=pickscann;
            linscan{:,wf}=linfitscann;
            frqscan(:,wf)=frqscann;
            rscan(:,wf)=rsscann;
        end
        figure
        %using a modification of the method proposed by pipo answering Yves Gaudemer's
        %question regarding "imagesc,pcolor,and NaN" in which pipo detailed
        %how to plot NaNs with a different color scheme than their
        %surrounding data in imagesc. Pipo's method is used here until the
        %line [rwind,rlaps]=find(rscan==max(max(rscan))), then from
        %imagesc(x,y,Q,clims) to the following xlabel, and then from
        %imagesc(x,y,rscan) to its following xlabel.
        %http://www.mathworks.com/matlabcentral/newsreader/view_thread/140607.
        %The matlab interal help documentation helped me learn about color
        %map and how to choose black as the color to replace the NaN
        %placeholder value sections. Bill Cheatham's anwers to Mercury's
        %question regarding "Matlab imagesc plotting nan values" in
        %Stackoverflow gave me the idea to use black as a color and how to
        %pursue doing it
        %http://stackoverflow.com/questions/14933039/matlab-imagesc-plotting-nan-values.
        %On february 24th I noticed that clims might cause the use of the
        %largest value set at >100000000000 for the NaNs or positive slope
        %attributed values for Q and rscan (not NaNs for rscan since there would be already an error prompt) (i.e. values given 123456789)
        %not to fuction using Pipo's method. I talked to Greg and he
        %informed me that this was the case and that a good alternative
        %would be to set the nans to zero and use zero as the smallest
        %value and to be colored black. I set all the 123456789 (for NaN
        %and positive slope attributed values, hopefully not others) to
        %zeros for Q and changed the colormaps so that the first row would be
        %black which corresponds to the smallest values (aka. zero). For
        %rscan I first changed the 123456789 values to -10000 for the purposes of the second subplot but then I changed the -10000 to zeros for the third subplot to use the method Greg and Pipo detailed. 
        [Qw,Ql]=find(Q==123456789);
        Q(Qw,Ql)=0;
        [rw,rl]=find(rscan==123456789);
        rscan(rw,rl)=-10000;
        [rwind,rlaps]=find(rscan==max(max(rscan)));%-from matlab help documentation
        subplot(5,1,1)
        y=[wndscan(1,wf),wndscan(end,wf)];%start the x axis at the first value of the coda window length (wndscan(1,wf) which should corrspond to the first value of the codawindow array).
        x=[Sscan(cf,1),Sscan(cf,end)];%start the y axis at the first value of the of the coda lapse time multiplier array (Sscan(cf,1) which should correspond to the fist value of the S array). 
        clims=[0 1000];
        imagesc(x,y,Q,clims);
        qd=colormap;
        qd(1,:)=[0,0,0];
        colormap(qd);
        ylabel('Coda window length (seconds)')
        xlabel('Coda lapse time mulitplier')
        titletop=(strcat('date',eventday,'time',timea,'frequencies',sprintf('%d',wn1a(freqscan{st}(stnind(gf),df,1,1))),'-',sprintf('%d',wn2a(freqscan{st}(stnind(gf),df,1,1)))));
        titlebot=(strcat('station',stationscan{st}(:,stnind(gf),df,1,1)'));
        titleg=strcat(titletop,titlebot); %create a title variable that is composed of the concatenated parts of the title which will be used in the save name of the image
        title({titletop;titlebot},'FontSize',8) %format the title to be two lines with a font size of 8 
        colorbar('eastoutside') %create a color bar on the right hand side of the image
        subplot(5,1,2)%plot the log(squared traces*time^(3/2))for the smallest coda lapse time multiplier and the longest coda window length then plot the corresponding fit on the same graph
        if isequal(max(max(rscan)),-10000)==0
        hold on
        plot((1:length(lpickscan{rlaps}{rwind}))/125,lpickscan{rlaps}{rwind},'b')%plot the log(squared traces*time^(3/2))for the smallest coda lapse time multiplier and the longest coda window length in blue
        plot((1:length(linscan{rlaps}{rwind}))/125,linscan{rlaps}{rwind},'r')%plot the corresponding fit on the same graph in red
        xlabel('Time from coda window start(seconds)')
        ylabel('Log of the product of the trace squared amplitudes multiplied by time^(3/2)')
        title(strcat('Fit and log(smoothed squared vertical component traces*time^(3/2)) for the best R^2 fit with Q: ', sprintf('%d',Q(rwind,rlaps)),', R^2',sprintf('%d',rscan(rwind,rlaps))),'FontSize',6)
        %special thanks to Framerius for answering Devon's question
        %on Stackoverflow regarding "Custom x-axis values in a matlab plot" for describing
        %how to change the x axis values using set-this method was used in
        %this code.
        %http://stackoverflow.com/questions/13568190/custom-x-axis-values-in-a-matlab-plot.
        %This method was used for subplots 2,4,and 5. I also had help from the MATLAB R2015a internal documentation titled "Changing Tick Marks and Tick Labels of Graph". 
        set(gca,'XTick',1:length(lpickscan{rlaps}{rwind})/125)
        hold off
        end
%         subplot(5,1,3)%plot the log(squared traces*time^(3/2)) for the largest coda lapse time multiplier and the largest coda window length then plot the corresponding fit on the same graph 
%         hold on
%         plot(1:length(lpickscan{end}{end}),lpickscan{end}{end},'b')%plot the log(squared traces*time^(3/2)) for the largest coda lapse time multiplier and the largest coda window length in blue
%         plot(1:length(linscan{end}{end}),linscan{end}{end},'r')%plot the corresponding fit on the same graph in red
%         xlabel('Time from coda window start(1/fs seconds)')
%         ylabel('Log of the product of the trace squared amplitudes multiplied by time^(3/2)')
%         title ('Fit and log(smoothed squared vertical component traces*time^(3/2)) for the largest coda laspe time and coda window length')
%         hold off
%         [sQ,Qq]=size(Q);
%         [sW,wq]=size(wndscan);
%         [sSc,Scq]=size(Sscan);
%         Q1=reshape(Q,[1,sQ*Qq]);
%         wndscan1=reshape(wndscan,[1,sW*wq]);
%         Sscan1=reshape(Sscan,[1,sSc*Scq]);
%         scatter3(wndscan1,Sscan1,Q1); %special thanks to Walter Roberson for answering tehslax question regarding 2d intensity maps which inspired the use of scatter http://www.mathworks.com/matlabcentral/newsreader/view_thread/289546 as well as the mathworks help documentation which detailed scatter3
%         xlabel('Coda window length (seconds)')
%         ylabel('Coda lapse time multiplier')
%         zlabel('Best fit Q value')
        subplot(5,1,3) %plot the R^2 values for eth particular station and filter 
        y=[wndscan(1,wf),wndscan(end,wf)];%start the x axis at the first value of the coda window length (wndscan(1,wf) which should corrspond to the first value of the codawindow array).
        x=[Sscan(cf,1),Sscan(cf,end)];%start the y axis at the first value of the of the coda lapse time multiplier array (Sscan(cf,1) which should correspond to the fist value of the S array). 
        [rtw,rtl]=find(rscan==-10000);
        rscan(rtw,rtl)=0;
        clims=[0 1];
        imagesc(x,y,rscan,clims) %the use of imagesc was described in the MATLAB help documentation site from mathworks.com titled "imagesc" as well as in the mathworks internal documentation for R2015a
        ylabel('Coda window length (seconds)')
        xlabel('Coda lapse time multiplier')
%         zlabel('Coefficiant of determination R^2')
        rd=colormap;
        rd(1,:)=[0,0,0];
        colormap(rd);
        colorbar('eastoutside')
        subplot(5,1,4)%plot the unsmoothed squared vertical component traces for the smallest coda laspe time and largest coda window length
        plot((1:length(pckscan{1}{end}))/125,pckscan{1}{end},'b')
        xlabel('Time from coda window start(seconds)')
        ylabel('Unsmoothed squared trace amplitudes')
        title ('Unsmoothed squared vertical component traces for the smallest coda lapse time and largest coda window length')
        set(gca,'XTick',1:length(pckscan{1}{end})/125)
        subplot(5,1,5)%plot the unsmoothed squared vertical component traces for the largest coda laspe time and largest coda window length
        plot((1:length(pckscan{end}{end}))/125,pckscan{end}{end},'b')
        xlabel('Time from coda window start(seconds)')
        ylabel('Unsmoothed squared trace amplitudes')
        title ('Unsmoothed squared vertical component traces for the largest coda lapse time and largest coda window length')
        set(gca,'XTick',1:length(pckscan{end}{end})/125)
        set(gcf,'Position',[200,200,2000,1000]) %set described by Dr. Greg
        saveas(gcf,titleg,'jpeg')
        close Figure 1
        if tip==1 % use an if statement to create the first instance of slapseroww and wccoll which will be the storage arrays for the indecies of the coda laspe time multiplier and the coda window length, repectively, that correspond to R^2>=0.7
        slapseroww=tip; %since tip=1 only occurs once the slapseroww and wccoll arrays will never be reset to one. This allows us to establish the variables with an intitial value of 1 which will later be removed at the end.
        wccoll=tip; 
        end 
        [wccol,slapserow]=find(rscan>=0.7); %find all of the R^2 values in rscan array which are greater than or equal to 0.7 and save the indices of the S and coda window length which correspond to those values, the value of 0.7 was from Calvert et al. 2013  
        slapseroww=cat(1,slapseroww,slapserow);  %update the slapseroww array which stores the coda laspe time multiplier index arrays for each iteration
        wccoll=cat(1,wccoll,wccol);%update the wcoll array which stores the coda window length index arrays for each iteration
        tip=tip+1; %update the tip count
    end 
end 
clearvars -except tip start Sadd Cdwndadd nanproblem nancount slapseroww wccoll st wn1a wn2a pickscan logpickscan linfitscan Qscan rsscan freqscan stationscan lapsescan windowscan dayscan timescan stnindx eventdate timeat
st=st+1;
end
% http://www.mathworks.com/matlabcentral/newsreader/view_thread/289546
% http://stackoverflow.com/questions/2938775/extracting-data-points-from-a-matrix-and-saving-them-in-different-matrixes-in-ma
 slapseroww=slapseroww(2:end)+Sadd;
 wccoll=wccoll(2:end)+Cdwndadd;
 figure
 subplot(2,1,1)
 histogram(slapseroww)
 xlabel('Lapse time multiplier values')
 ylabel('Frequency counts')
 title('Frequency of lapse time multipliers that result in an R^2>=0.7')
 subplot(2,1,2)
 histogram(wccoll)
 xlabel('Coda window lengths')
 ylabel('Frequency Counts')
 title('Frequency of coda window lengths that result in an R^2>=0.7')
 set(gcf,'Position',[200,200,2000,1000]) %set described by Dr. Greg
 saveas(gcf,'laspemultadnwindowhistograms','jpeg')
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 %         for wind=1:(length(vertdat)/smoothingwindow)-1 %smoothing window as described in Calvet et al. %check to see if the 1 offset is truly necessary
%            vertdat(1:smoothingwindow)=hann(smoothingwindow).*(vertdat(1:smoothingwindow));
%            vertdat((wind*smoothingwindow)+1:(wind+1)*smoothingwindow)=hann(smoothingwindow).*vertdat((wind*smoothingwindow)+1:(wind+1)*smoothingwindow);
%        end
 
 
 
 
 
%%
%rotating the data- this code was inspired/adapted from Gregory Waites Volcano seismology class code: Lab 7 VLP particle motion 
% for k=1:(count/3) %since the North and east are to be rotated into radial and tangential compoenents this for loop allows "jumps of three regular indicies" to allow for loop to more easilty index the proper data
% %rotate the two horizontal components into radial and tangetial
% % use a previously constructed file that contains the latitude and
% % longitude coordinates for the stations
% stationloc=char(staloc(k)); %use the staloc cell (a cell list of participating stations in the order they were read from the array generated from the cnv file) to produce a string referenced to k which ranges from 1 to the number of stations that are involved with the event picks
% if length(staloc)~=count/3 %show an error message if the amount of recorded station names does not match with count/3; the amount of stations involved
%     error('number of station names does not correspond to count/3 i.e. number of stations')
% end
% stanumb=strcmp(stavec,stationloc); %search the stavec (the array of station names) for the name of the station corresponding to k
% stanum=find(stanumb==1); %once found obtain the index of the station in the array-this index corresponds to the index of the desired station locations in the imported station locations excel data (the stavec and the excel data correspond).
% %The imported station location fields will import under the vector names of
% %Lattitude and Longitude, use stanum to find the correct values for the
% %station currently under rotation
% stationlat=Latitude(stanum,1); 
% stationlong=Longitude(stanum,1);
% %This takes the azimuth from the source to the station and thus does not require the back azimuth equations
% angle=azimuth(eventlat,eventlong,stationlat,-stationlong,'radians'); % the longitudes are negative since we are in the western hemisphere-Greg, the eventlong is already negative
% %Rotating the angles into the radial and the transverse
% %pre rotation matrix
% Prmatsignal=[pickevent(:,3*k)';pickevent(:,(3*k)-1)'];  %if the components are in the right order the first row of this vector is the east component pick for the station in question and the second row ith the picked data of the north component
% Prmatnoise=[picknoise(:,3*k)';picknoise(:,(3*k)-1)']; %the same as said for the Prmatsignal goes for the Prmatnoise 
% % rotation matrix
% rot=[cos(angle),sin(angle); -sin(angle), cos(angle)]; % in this case the radial is pointing from the station to the source so is negative
% %rotate the matrix
% Aevent=rot*Prmatsignal;
% Bnoise=rot*Prmatnoise;
% %extract the radial and tangential components 
% radialsignal(k,:)=Aevent(1,:); %the first row of the Aevent  will be the radial component
% tangentialsignal(k,:)=Aevent(2,:); %the second row of the Aevent will be the tangential component
% radialnoise(k,:)=Bnoise(1,:); %the same appiles ot the radial and tangential noise
% tangentialnoise(k,:)=Bnoise(2,:);
% clear Aevent
% clear Bnoise
% 
% %generate the spectra for the radial and tangential components as well as
% %the vertical components
% 
%  spectrumsignal((k*3)-1,:)=fft(radialsignal(k,:),1024);  %spectrum of radial signal
%  spectrumsignal((k*3),:)=fft(tangentialsignal(k,:),1024); %spectrum of tangential signal
%  spectrumnoise((k*3)-1,:)=fft(radialnoise(k,:),1024);  %spectrum of radial noise
%  spectrumnoise((k*3),:)=fft(tangentialnoise(k,:),1024); %spectrum of tangential noise
%  %compute the spectrum of the vertical signal and noise
%  spectrumsignal((k*3)-2,:)=fft(pickevent(:,(3*k)-2),1024); %spectrum of vertical signal 
%  spectrumnoise((k*3)-2,:)=fft(picknoise(:,(3*k)-2),1024); %spectrum of vertical noise
%  end
% countorig=count; % create a scalar with the count number so that the variable count can be repurposed
% %%
% % grid search for tstar of the event
% for grid=1:countorig
%     %use the number of events to create a for loop that brings up the
%     %spectrums of the postrotated matrix so that each spectrum can be
%     %anaylsed using the grid search.
%     %remove the loaded and indexed noisespectum corresponding to the current event and component 
%     spectrums(:,grid)=spectrumsignal(grid,:)'-spectrumnoise(grid,:)';  %remove the noise from the data-post rotaion
%     %create a frequency vector from 1 to fs 
%     ff=linspace(0,fs,length(spectrumsignal(grid,:)'));
%     %use a counting system to generate an index for each individual
%     %combination of variables
%     clear count
%     count=1;  
%     %create a for loop to determine the best fit for the zero frequency amplitude
% for do=50000:5000:300000
%     %create a for loop to find the best fit for the corner frequency 
%         for fc=2:1:40
%             %create a for loop to find the best fit for the tstar
%             for tz=(1*10^-2):(1*10^-2):(40*10^-2)
%                     %av(count)=a; %tstarzero and alpha have been changed
%                     %for tstar since we are going with alpha=0
%                     tzv(count)=tz;  %use this structure to save variables using the counting scheme
%                     fcv(count)=fc;
%                     dov(count)=do;
%                     for fi=1:length(spectrums(:,grid)) %The length will reflect the length of the fft. If not then there is a problem 
%                         f=ff(fi); %create a frequency variable from the ff- frequency vector                    
%                         tstar=tz;
%                         H(fi)=(do*exp(-pi*f*tstar))/(1+(f/fc)^(2*2))^(1/2);%brune source model. This model may be innapropriate
%                         df(fi)=abs(H(fi))-abs(spectrums(fi,grid)); %subtracting the magnitude of the the model spectra from that of the observed spectra 
%                         
%                     end 
%                     Hmat(:,count)=H(41:end); %create a matrix of brune source model spectrums
%                     clear H
%                     % create a list of root mean squared results of the differences between the brune model and the observed spectrum
%                     rmsq(count)=sqrt((sum(df(41:end).^2))/length(spectrums(41:end,grid))); %cut anything below 5.0049 Hz out of the analysis  
%                     clear df
%                     %update the count 
%                     count=count+1;
%             end 
%         end
% end
% %find the minimum root mean square
% bv(grid)=min(rmsq);
% %find the index of the minimum root mean square 
% ref=find(rmsq==bv(grid));  
% %alpha(:,grid)=av(ref(1,grid));
% tstarzero(:,grid)=tzv(ref(:)); %since the rmsq has the same index numbers as the variables it can be used to find the apropriate values  
% freqenycorner(:,grid)=fcv(ref(:));
% zerofreqency(:,grid)=dov(ref(:));
% if length(tzv)~=length(rmsq) | length(tzv)~=length(fcv) | length(tzv)~=length(dov)
%     error('reference system fails')
% end 
% brune=Hmat(:,ref); %create a vector of amplitudes of the model spectrum
% if length(ref)>1
%     disp('mutliple brune models provide minimum root mean square differences')
% end
% brunea(:,:,grid)=brune; %save the vector in a column format with the number of columns equaling the number of components
% figure
%                     hold on
%                     %plot the model spectrum with the observed spectrum
%                     plot(ff(41:end),abs(brune),'r') %plot the magnitudes of the model data
%                     xlabel('Freqencies (Hertz)')
%                     ylabel('Amplitude')
%                     nmtc=char(titleall(1,grid)); %obtain the component name 
%                     for up=1:countorig/3 %create a for loop with the number of components/3 allows for the component reference as done with the pre rotation matricies
%                     if grid==(up*3)-2 %create the appropriate names-make sure these correspond to the events/components originally chosen
%                         nam=nmtc; %if vertical component
%                     end
%                         if grid==(up*3)-1 %if north component
%                             nam=strcat('radial event:',nmtc(6:end));
%                         end
%                             if grid==(up*3) %if east component
%                                nam=strcat('tangential event:',nmtc(5:end));
%                             end 
%                     end
%                             tl=strcat('Spectrum event (blue) with calculated (red)',nam,'rmsq:',sprintf('%d %s',bv(grid))); %create the title of the plot
%                     title(tl)
%                     plot(ff(41:end),abs(spectrums(41:end,grid)),'b') %plot the magitude of the observed spectra
%                     xlim([0,50])                                        %Place a limit on the frequency 
%                     hold off
%                      clear brune %clear the temporary model vector
% clear rmsq nmtc
% clear ref
% end
% 
    %find the corner of the frequency spectrum using another if/then system
    %to detemine the overall change of slope for the spectrum fit 
    %cornerth is the threshold for the corner detection.
   %%
%    cornerth=0.5;
%     for gridt=1:countorig;
%         lza=length(spectrums(:,gridt));
%         [amp(gridt),loc(gridt)]=findpeaks((spectrums(gridt)));
%         figure 
%         hold on
%         plot(ff,spectrums(1:end/2,gridt),'b')
%         xlabel('Frequency (Hertz)')
%         ylabel('Amplitude')
%         title('Sprectums with psuedoenvelope')
%         plot(loc(1:end/2),amp(1:end/2),'g')
%         hold off
%         spectrumssm(:,gridt)=fit(loc,amp,'poly1');
%         for wcorn=1:lza-8
%             y2=spectrumssm(wcorn+8,gridt); y1=spectrumssm(wcorn,gridt); x1=wcorn; x2=wcorn+8;  
%             slopea=(y2-y1)/(x2-x1);
%             if abs(slopea)>=cornerth 
%                y2b=spectrumssm(wcorn+4,gridt); y1b=spectrumssm(wcorn,gridt); x1b=wcorn; x2b=wcorn+4;  
%                slopeb=(y2b-y1b)/(x2b-x2b);
%             if abs(slopeb)>=cornerth 
%                y2c=spectrumssm(wcorn+2,gridt); y1c=spectrumssm(wcorn,gridt); x1c=wcorn; x2c=wcorn+2;  
%                slopec=(y2c-y1c)/(x2c-x1c);
%             if abs(slopec)>=cornerth
%                corner(gridt)=wcorn+1;
%             else
%                corner(gridt)=wcorn+3;
%             end
%             else
%                 y2b2=spectrumssm(wcorn+6,gridt); y1b2=spectrums(wcorn+4,gridt); x1b2=wcorn+4; x2b2=wcorn+6;  
%                 slopeb2=(y2b2-y1b2)/(x2b2-x1b2);
%              if abs(slopeb2)>=cornerth   
%               corner(gridt)=wcorn+5;
%              else
%               corner(gridt)=wcorn+7; 
%              end
%             end 
%             end
%         end
%      end

%Plot  all of the compoenet and observed spectra on one plot with obserevd spectra in blue and modeled spectra in red 
% figure
% hold on 
% plot(ff,abs(spectrums),'b')
% xlabel('Frequencies (Hertz)');
% ylabel('Amplitude')
% title('Total event spectra: observed vs modeled')
% plot(ff,abs(brunea),'r')
% xlim([0,50])
% hold off


%Thanks to Matlab Answers:Dan Ryan and Jill Reese: http://www.mathworks.com/matlabcentral/answers/62382-matrix-multiply-slices-of-3d-matricies?s_tid=srchtitle, 
%Special Thanks to those who contibuted to the Stackoverflow and the
%Mathworks help site for their insight into code functionality and term
%usage. The process of learning about new methods of coding frequently
%referenced the works found in Stackoverflow and Mathworks with the
%knowledge of code functionality then applied in this work.
%--------------------check that these help site are correct. Special thanks
%to other MATLAB help sites such as blogs.Mathworks and answers.MATHWORKS for helping to instruct
%the author in the use of MATLAB and by providing coding examples------. Special thanks to the MathWorks MATLAB
%help library, both online and included with the program as many of its
%examples were learnt form and used in this code. Please note that
%adaptaions of codes displayed in matlab help sites were used in this work.
%Special thanks to Federica Lanza for her help with the coding and by providing the data from which the results of this code were derived.
%Special thanks to Cassandra Javor and Thomas and Cynthia Guettinger for
%their advice and support.
%The method of finding Q such as the equation used in the grid search and
%the were derived from De sienna and Calvert. This code was based off De
%Sienna and Clavert's papers.
%The overall architecture of the code was designed by Dr. Gregory Waite.