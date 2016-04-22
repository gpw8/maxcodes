%final awesome Q code
%see feb15th code for citations if not included here

% _2016_03_03

S=4
codawindow=10
Snumber=S;
nancount=1;
%begin code copied from the feb 15th code
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
%skip a few lines
Aw(:,finddate(2:end)-1)=[];%this removes the zero columns (columns with a single zero value at their begining) that acted as buffers at the end of events-method from the internal program documetation for MATLAB R2015a titled Deleting Data from a Cell Array
Aw(:,end)=[];
%end code copied from the feb 15th code
for mg=1:length(finddate) %this loop should loop from 1 to the length of finddate
finddate=find(cellfun('isempty',regexp(Aw(1,:),'1501[12]'))==0);%this should find the dates assuming that nothing else contains the '15011' or '15012' sequences.-special thanks to 'YYC' answering 'bit-questions' question regarding cells in Stackoverflow that presented regexp as an option to find the cell values with the [12] in the string http://stackoverflow.com/questions/8056131/strfind-for-string-array-in-matlab. Also special thanks to the internal MATLAB help for the R2015a program used at school for helping with its detail on regexp and how to use it to find patterns in strings,  %special thanks to Jonas answering N.C.Rolly's question about finding empty cell arrays on stackoverflow http://stackoverflow.com/questions/3400515/how-do-i-detect-empty-cells-in-a-cell-array   
eventday=Aw(1,finddate(mg)); %use finddate with the mg variable as an index to obtain the event date values in Aw in the order observed in the f24 file----check
timea=Aw(2,finddate(mg)); %use the same method as the previous line but to find the timea variables in Aw(2,)-check 
timeachksec=char(Aw(3,finddate(mg)));            %start time seconds for each event
Aqt=finddate(mg); %set Aqt equal to finddate(mg)
eventdaychk=char(eventday);        %check that the day number for the event is ok-char method can be found in matlab documentaion

%begin code copied from the feb 15th code
yearst=eventdaychk(1:2);      %year of the event 
monthst=eventdaychk(3:4);     %month of the event
dayst=eventdaychk(5:6);       %day of the event
%end code copied from the feb 15th code

timeachk=char(timea); %find the time of the event

%begin code copied from the feb 15th code
strhr=timeachk(1:2);              %start time hour for event
strmin=timeachk(3:4);             %start time minutes for each event
%end code copied from the feb 15th code

%begin code copied from the feb 15th code
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
    if Aqt+Rt<=length(Aw);
    gh(Rt)=strncmp(Aw(1,Aqt+Rt),Aw(1,Aqt),3); %compare the first three string characters of the first cell in each column for Rt colunms for similarity to find the next event time
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
       %end copied section from feb 15th code 
  
  ldnm=strcat('dcnvlv',qs,compnm(1:3),'day',eventdaychk); %use the current station name +the vertical designation (compnm) and the event dat to create a load name with the statement 'dcnvlv' added in
    ldnnm=strcat('/run/media/',ldnm,'.mat'); %create an overall path-loadname for the data-this will differ based on where the data is stored
    
    %begin code copied from the feb 15th code
    load(ldnnm); %load the data
    namd=strcat(qs,compnm); %make a name with just the currentr station name and the vertical designation
    namd=namd(1:7); %cut the "20" out of the namd name
    nmind=strcmp(stacompvec,namd);  %reference the name with the properly ordered (with respect to the Poles Zeros, and Sensitivity arrays) stacompvec to obtain an the array location which corresponds to that of the correct poles, zeros, and sensitivity
    nmind=find(nmind==1);  %find the index value of the previously obtained array location found in the last line
       
    %end copied section from feb 15th code
    data=deconvolvedd; %access the current loaded data structure and rename the data into a new variable, knowledge of how to access data structures from Greg Waite Earthquake seismology course and Mathworks MATLAB help doucmentation
    fs=125; %obtain the sample rate from the currently loaded data structure

    %begin code copied from the feb 15th code
    namdd(:,countt)=namd; 
    %skip some lines in the feb 15th code
    d=data(eventbeginsample-1250:eventbeginsample+((codastart+codawindow)*fs)+1250); %rename the data and add 10 seconds (should be) worth of data before and after the event begins that will later be cut to act as a buffer against edge effects-greg
    %skip some lines in the feb 15th code
    for filt=1:6 %this creates a filter loop which is used to filter the data-deconvolve it and then pick the data according tot he coda window start time (laspe time) and coda window length. It also plots the deconvolved and the picked data as well ast the sqaure of the picked data-calvert?  
     wn1a=[1,2,4,8,16,20]; %lower filter ranges  filters ranges from calvert
     wn2a=[2,4,8,16,32,60]; %upper filter ranges
     centerfreq=[1.5,3,6,12,24,40];
     %end code copied from the feb 15th code
nyquist=fs/2;
     deconvolved=d;
     
    %begin code copied from the feb 15th code
     wn1=wn1a(filt); %select the lower filter according to the loop progression
     wn2=wn2a(filt); % select the upper filter according to the loop progression 
    [FB,FA]=butter(4,[(wn1/nyquist) (wn2/nyquist)],'bandpass');  % the concept and implimentation of the butterworth filters was derived from Dr. Gregory Waite's filt_traces.m function which was based off a function written by Derek Schutt. -filter help by mathworks matalb documetation- greg informed me to set the order to 4 but not really higher to avoid edge effects
    deconvolved(:,1)=filter(FB,FA,deconvolved); % we use a bandpass filter that has variable frequency ranges as dictated by the for loop following the method by Calvert et al. 2013
    %end copied section from feb 15th code
    
    dh=strcat('filt',ldnm,'frequencies',sprintf('%d',wn1),'-',sprintf('%d',wn2));

    %begin code copied from the feb 15th code
    clear FA FB 
    %skip some lines in the feb15th code
    vertdat=deconvolved(codastart*fs+1251:1251+(codastart+codawindow)*fs,1).^2; 
    %skip some lines in the feb15th code
    titlev=dh; %use the name and the event day to make the title for the first plot which should be the vertical
    titleall{1,count}=titlev; %store the title in a cell array of titles that increase with each count   
    %skip some lines in the feb 15th code 
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
       %end copied section from feb 15th code
       clear deconvolved data namd nmind ldnnm pA
       %begin copied code from feb 15th code 
 end
    eventbeginsample=eventbegin*fs;
       clear vertdat cell d
       countt=countt+1;
       if Ki~=(wt(1)-1) %use this if statement to avoid using the last string cell of the column as a station names since it is the travel time of the station at wt-2
       qs=qsta(5:8); %obtain the station name for the next iteration (for the first iteration this replaces the name of the station found in the top cell of the column with the name in the current cell, to be use with the travel time in the next cell lower of the column)
       end 
         end
clear qs %clear qs upon the end of the pick loop
end
count=count-1;

%copy some of the reorganizing variables while removing those that are not
%needed to obtain Q
for m=1:6 %-coding lines dealing with the possible S wave maximum (or maximum value between the event origin time and the start of the coda window) were proposed by Dr. Greg Waite.
pickeventa(:,:,m)=pickevent(:,m:6:end); %-creating multidimensional vectors method found in matlab help. More can be found on Stackoverflow asked by Theodoros Theodoridis http://stackoverflow.com/questions/23376111/multidimensional-arrays-multiplication-in-matlab
vertsavea(:,:,m)=vertsave(:,m:6:end);%reorganize the saved unsmoothed vertdat arrays in a similar fashion to the pickeventa
titlealla(1,:,m)=titleall(1,m:6:end);
end
nstns=length(pickeventa(1,:,1));
kont1=1;

 %begin copied code from feb 15th code 
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
          pnan=find(isnan(pickcheck(:,kount1))==1); %for the particular station and frequency find the NaN values using isnan on pickcheck(:,kount1) and then using find(isnan(pickcheck(:,kount1))==1) to find the indicies of the NaN values and store them in a variable special thanks to Marc for
        %answering Graviton's question on Stackoverflow regarding finding NaN values for it
        %demonstrated how to perform this action with Graviton's method used in
        %this code http://stackoverflow.com/questions/1713724/find-all-nan-elements-inside-an-array
          pickcheck(pnan,kount1)=0; %set all of the NaN values to 0, repeat this method for the logpick and linfittocheck arrays  
          logpick(:,kount1)=Et;
          lonan=find(isnan(logpick(:,kount1))==1);
          logpick(lonan,kount1)=0;
          linfittocheck(:,kount1)=linfit;
          linan=find(isnan(linfittocheck(:,kount1))==1);
          linfittocheck(linan,kount1)=0;
%skip a few lines
          rsend(kount1)=R2;%save the R^2 value in an array which updates with each passing kount1, with the data ultimately being used in the subplots
          titlend(kount1)=titlealla(1,numstations,freqsza); %save the title for the particular station and frequency 
          freqend(kount1)=freqsza; %save the frequency index
          stationend(:,kount1)=namdd(:,numstations); %for a particular station and frequncy save the stationname which will be used later in the subplots, since the namdd variable contains the names of the stations and components in the order in which they wer deconvolved they should be properly indexed by numstations per event 
          %skip a few lines
        kount1=kount1+1; %update the kount1 variable
             else 
                 disp('positive slope for fit, data modified for station iteration')
                 nname=strcat('S: ',sprintf('%d',S),', codawindow: ',sprintf('%d',codawindow),', Freq: ',sprintf('%d',freqsza),' ,stn: ',namdd(:,numstations)');
                 nanproblem{nancount}=nname;
                 nancount=nancount+1;
                 Qtocheck(kount1)=123456789;
                 pickcheck(:,kount1)=vertsavea(:,numstations,freqsza); %use the pickcheck variable to store the vertsavea data for a particular station and frequency
                 pnan=find(isnan(pickcheck(:,kount1))==1); 
                 pickcheck(pnan,kount1)=0;
                 logpick(:,kount1)=repmat(123456789,[length(Et),1]);
                 linfittocheck(:,kount1)=repmat(123456789,[length(K2),1]);
                 rsend(kount1)=123456789;
                 titlend(kount1)=titlealla(1,numstations,freqsza); %save the title for the particular station and frequency 
                 freqend(kount1)=freqsza; %save the frequency index
                 stationend(:,kount1)=namdd(:,numstations); %for a particular station and frequncy save the stationname which will be used later in the subplots, since the namdd variable contains the names of the stations and components in the order in which they wer deconvolved they should be properly indexed by numstations per event 
                %skip a few lines
             kount1=kount1+1;
             end 
        clear R2 ploye Et Qc linfit sr stot K K2 lonan linan pnan nname% clear the variabels that are no loger necessary after this iteration.
%          close Figure 1
     end
%skip a few lines
        Qtochecka(:,kont1)=Qtocheck;
        linfittochecka(:,:,kont1)=linfittocheck;
        pickchecka(:,:,kont1)=pickcheck;
        loggpick(:,:,kont1)=logpick;
        rsenda(:,kont1)=rsend;
        titlenda(:,kont1)=titlend;
        freqenda(:,kont1)=freqend;
        stationenda(:,:,kont1)=stationend;
%skip a few lines
clear kount1 Qtocheck rsend titlend freqend stationend laps windw logpick pickcheck linfittocheck
         kont1=kont1+1;
%         close all %remmeber to disable the close alls when plotting
 end
Q(:,:,mg)=Qtochecka;
linefit{mg}=linfittochecka;
pickdat{mg}=pickchecka;
logdat{mg}=loggpick;
r(:,:,mg)=rsenda;
title{mg}=titlenda;
freqc(:,:,mg)=freqenda;
stnnm{mg}=stationenda;
end