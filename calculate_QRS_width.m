function [medianQRSw,meanQRSw,stdQRSw]=calculate_QRS_width(ecg,posp,fs)
% Computes median, mean and standard deviation (SD) of QRS complex width
%
% INPUT:
% - ecg: preprocessed/denoised ECG
% - posp: instants in seconds of the QRS complexes detected by Hamilton-Tompkins (HT) detector 
% - fs: sampling rate of ECG
% OUTPUT:
% - medianQRSw: median width of the QRS complexes
% - meanQRSw: mean width of the QRS complexes
% - stdQRSw: SD of QRS complex width
%
% Original code by Erik Alonso

aposR=[];aposQ=[];aposS=[];
frac=1.20;
MPH=-0.1;
th=0.05;

decg=diff(ecg);decg=[decg decg(end)];
pos=round(posp*fs);
for cc=1:length(pos)

    [vmin,pmin]=findpeaks(-ecg,'MinPeakHeight',MPH);
    [vmax,pmax]=findpeaks(ecg,'MinPeakHeight',MPH);
    
    % Closest local maximum
    [valMax,posiMax]=min(abs(pmax-pos(cc)));
    cPosMax=pmax(posiMax);
    
    % Closest local minima
    [valMin,posiMin]=min(abs(pmin-pos(cc)));
    cPosMin=pmin(posiMin);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Positive segment: previous local minimum - local maximum closest to QRS instant detected by HT - next local minimum 
    
    % Previous local minimum
    posQ=pmin(pmin<cPosMax);
    if isempty(posQ)
        if cc==1
            [vd,po]=min(ecg(1:cPosMax));
            posQ=po;
        else
            inter=round(frac*(aposR(end)-aposQ(end)));
            [vd,po]=min(ecg(max(1,cPosMax-inter):cPosMax));
            posQ=max(1,cPosMax-inter)+po-1;
        end
    else
        posQ=posQ(end);
    end
    
    % Next local minimum
    posS=pmin(pmin>cPosMax);
    if isempty(posS)
        if cc==length(pos)
            [vd,po]=min(ecg(cPosMax:end));
        else
            inter=round(frac*(aposS(end)-aposR(end)));
            [vd,po]=min(ecg(cPosMax:min(length(ecg),cPosMax+inter)));
        end
        posS=cPosMax+po-1;
    else
        posS=posS(1);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Negative segment: previous local maximum - local minimum closest to QRS instant detected by HT - next local maximum  
    
    % Previous local maximum
    posQ2=pmax(pmax<cPosMin);
    if isempty(posQ2)
        if cc==1
            [vd,po]=max(ecg(1:cPosMin));
            posQ2=po;
        else
            inter=round(frac*(aposR(end)-aposQ(end)));
            [vd,po]=max(ecg(max(1,cPosMin-inter):cPosMin));
            posQ2=max(1,cPosMin-inter)+po-1;
        end
    else
        posQ2=posQ2(end);
    end
    
    % Next local maximum
    posS2=pmax(pmax>cPosMin);
    if isempty(posS2)
        if cc==length(pos)
            [vd,po]=max(ecg(cPosMin:end));
        else
            inter=round(frac*(aposS(end)-aposR(end)));
            [vd,po]=max(ecg(cPosMin:min(length(ecg),cPosMin+inter)));
        end
        posS2=cPosMin+po-1;
    else
        posS2=posS2(1);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % New R wave will be the local minimum or maximum located at the
    % central part of the segment which shows greater slope
    
    if mean(abs(decg(posQ:posS)))>mean(abs(decg(posQ2:posS2)))
        if min(ecg(cPosMax)-ecg(posQ),ecg(cPosMax)-ecg(posS))>th*max(ecg(cPosMax)-ecg(posQ),ecg(cPosMax)-ecg(posS))
            aposR=[aposR cPosMax];
            aposQ=[aposQ posQ];
            aposS=[aposS posS];
        else
            aposR=[aposR cPosMin];
            aposQ=[aposQ posQ2];
            aposS=[aposS posS2];
        end
    else
        if min(ecg(posQ2)-ecg(cPosMin),ecg(posS2)-ecg(cPosMin))>th*max(ecg(posQ2)-ecg(cPosMin),ecg(posS2)-ecg(cPosMin))
            aposR=[aposR cPosMin];
            aposQ=[aposQ posQ2];
            aposS=[aposS posS2];
        else
            aposR=[aposR cPosMax];
            aposQ=[aposQ posQ];
            aposS=[aposS posS];
        end
    end
    
    
end

width=(aposS-aposQ)/fs;
meanQRSw=mean(width);
medianQRSw=median(width);
stdQRSw=std(width);
