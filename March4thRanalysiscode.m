% _2016_03_28
%March4thRanalysiscode
% %different events and same stations use i=1:6
stns=[11,12,9,8,17,1;11,12,9,8,17,1;8,1,9,11,12,17;17,11,9,1,8,12;8,17,11,12,9,1;17,12,9,11,1,8];
stnindx=[2,4,5,6,7,9;2,3,8,10,11,15;1,3,4,5,6,9;1,4,6,7,8,10;1,3,5,11,12,13;2,3,4,6,7,9];
%different events and different stations use i=1:6
% stns=[15,5,12,13,16,19;13,19,15,5,16,12;16,19,15,5,13,12;5,12,16,13,15,19;12,16,13,15,19,5;12,15,5,13,16,19];
% stnindx=[1,3,4,10,14,17;3,8,9,12,13,15;4,5,6,8,11,13;2,5,12,14,15,18;3,11,13,14,15,18;2,12,13,14,15,17];
%original use i=1:10
%  stns=[17,8,11,9,1,12;9,12,1,17,8,11;11,12,8,17,9,1;8,1,11,12,9,17;8,17,12,11,9,1;11,12,9,1,8,17;11,9,1,8,17,12;11,8,12,1,9,17;1,8,17,11,9,12;1,11,17,9,8,12];
%  stnindx=[1,2,4,5,6,8;1,4,5,8,9,12;1,2,3,4,5,8;1,2,4,5,9,11;2,6,12,13,14,17;3,4,7,8,10,12;1,3,4,7,10,12;1,2,4,6,8,11;1,2,3,4,6,7;2,3,4,6,7,8]
for df=1:6
    for i=1:6 
        for gf=1:6
            stnind=stnindx(i,:); %select the correct list of indices for the stations that was prepared for the feb 15th code
            Rrscan(:,:)=rsscan{i}(stnind(gf),df,:,:); %create a 2 dimensional array or R^2 values by indexing rrscan for a particular event, frequency, and station. 
            if isequal(max(max(Rrscan(:,:))),123456789)==0 % find if the maximum Rrscan value equals 123456789 
            [ro,col]=find(Rrscan==max(max(Rrscan(:,:)))); %if it does not find the row and columns indicies of the maximum value
            Rsn(df,i,gf)=max(max(Rrscan(:,:))); %generate a matrix that stores the maximum Rrscan value
            if length(ro)>1 | length(col)>1 % if the dimensions of the row and column index is greate than 1 display an error 
                error('rows and columns dimensions exceed 1')
            end 
            QfRmax(df,i,gf)=Qscan{i}(stnind(gf),df,ro,col);
            else
                [roww,collie]=find(Rrscan==max(max(Rrscan(:,:)))); %in the case in which the maximim of Rrscan=123456789 we find the row and columns of the values that correspond to this number 
               for ndt=1:length(roww) %for loop created by Greg
                Rrscan(roww(ndt),collie(ndt))=0; %we then set the values in Rrscan whose rows and columns correspond to those of the values which equal 123456789 to zero
               end 
                clear roww collie %we then clear the rows and column indices
                [roww,collie]=find(Rrscan==max(max(Rrscan(:,:)))); %we then find the row and column indicies of the maximimum of the Rrscan matrix now that the values of 123456789 have been set to zero using the same method that we used to find the rows and columns of the previous value 
                if isequal(max(max(Rrscan(:,:))),0)==0 % we neeed to chech wether the maximum of the refined Rrscan matrix is zero
                if length(roww)>1 | length(collie)>1  % if the dimensions of the rows and columns exceed 1 each then we display an error
                error('rows and columns dimensions exceed 1')
                end 
                Rsn(df,i,gf)=max(max(Rrscan(:,:))); %if the maximum of the Rrscan matrix is not zero then we save it and the corresponding Q as well
                QfRmax(df,i,gf)=Qscan{i}(stnind(gf),df,roww,collie);
                else 
                    Rsn(df,i,gf)=0; %if the maximum of the uncorrected Rrscan was 123456789 and then the maximum of the correct Rrscan is 0 then set the Rsn for this iteration equal to zero; 
                    QfRmax(df,i,gf)=0; % if the maximum of Rrscan=0 then we set the corresponding Q value equal to zero, if the rscan is really 0 istead of being chanegd to zero by conversion of the 123456789 values then this may cause some data to be improperly interpreted/disregarded/changed 
                end 
            end 
            dfcount(df,i,gf)=df; %save the df values
            evntcount(df,i,gf)=i; %save the event values
            stncount(df,i,gf)=stns(i,gf); %save the station number 
            clear stnind roww collie ro col rowww collie Rrscan ndt
        end 
    end 
 end 
 clear df i gf
for df=1:6     
    figure
    RRsn(:,:)=Rsn(df,:,:);
    subplot(3,1,1)
    imagesc(RRsn(:,:))
    colorbar('eastoutside')
    ylabel('events')
    xlabel('station index')
    for ab=1:length(RRsn(1,:)) 
        for ae=1:length(RRsn(:,1))
            strR=num2str(RRsn(ae,ab));
            text(ab,ae,strR)
        end 
    end 
    sttitl=strcat('Maximum R for stations and events, Freq',sprintf('%d',df));
    title(sttitl)
    subplot(3,1,2)
    Qfs(:,:)=QfRmax(df,:,:);
    clims=[0 1000];
    imagesc(Qfs(:,:),clims)
    colorbar('eastoutside')
    ylabel('events')
    xlabel('station index')
    for at=1:length(Qfs(1,:)) 
        for ay=1:length(Qfs(:,1))
            strQ=num2str(Qfs(ay,at));
            text(at,ay,strQ)
        end 
    end 
    sttitle=strcat('Q for maximum R for stations and events, Freq:',sprintf('%d',df));
    title(sttitle)
    subplot(3,1,3)
    stnc(:,:)=stncount(df,:,:);
    imagesc(stnc(:,:))
    colorbar('Ticks',[1,8,9,11,12,17])
    ylabel('events')
    xlabel('station index')
        for al=1:length(stnc(1,:)) 
        for ah=1:length(stnc(:,1))
            strnc=num2str(stnc(ah,al));
            text(al,ah,strnc)
        end 
    end 
    title('Station number per event and station index')
end 
 %added march 8th
 for df=1:6
     for i=1:6
         stdR(df,i)=std(Rsn(df,i,:));
         stdQ(df,i)=std(QfRmax(df,i,:));
     end
 end 
 Rtable=table(stdR,'Rownames',{'1';'2';'3';'4';'5';'6'});
 Qtable=table(stdQ,'Rownames',{'1';'2';'3';'4';'5';'6'});
 
 %special thanks to online help resources such as the MATLAB online help documentation from Mathworks (as well as off line help documentation) as
 %well as MATLAB answers and probably to the contributers to Stackoverflow for their
 %help in code functionality. 
 

 