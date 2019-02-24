%Prepare to calculate 2nd order structure parameters of temperature, humidity, and
%the temperature humidity cross structure parameter

%Can add to calculate further structure parameters here at a later date.

function [output, raw] = UTES_struct_setup(data, info, output, raw, tableNames, sensorInfo)

if isempty(raw)
    warning('Could not calculate Structure parameters. raw data does not exist')
    return
end

fprintf('\nCalculating Structure Parameters\nC_T^2 C_Tq C_q^2\n');

Rd = 287.058;  % [J/K/kg] Gas constant for air

%Number of Sonics
NumSonic = size(sensorInfo.Tson, 1);
%Number IRGAs
NumIRGA = size(sensorInfo.irgaH2O, 1);

MaxInst = max([NumSonic, NumIRGA]);
MinInst = min([NumSonic, NumIRGA]);


%Iterate through number of instruments (heights & Tables)
cntrIRGA = 1;
for qq=1:MaxInst
    %TableName columns in UTESpac and siteInfo do not match....
    %These lines look for the tablename from siteInfo
    TableNum = sensorInfo.Tson(qq, 1);
    tmp = regexp(info.tableNames, tableNames{TableNum});
    [~, siteInfoTableNum] = max(~cellfun(@isempty, tmp));

    %Calculate number of points per averaging period
    pnts = info.tableScanFrequency(:, siteInfoTableNum)*60*info.avgPer;

    %Create average time vector
    time = raw.t(pnts:pnts:end);

    %Iterate through each aeraging period and calculate strucutre
    %parameters
    cntr1=1;
    cntr2 = pnts;
    for ii=1:length(raw.t)/pnts
        
        SonHeight = sensorInfo.Tson(qq, 3);
        %Sonic Variables for each height
        sigma_u = nanstd(raw.uPF(cntr1:cntr2, qq)); 
        sigma_v = nanstd(raw.vPF(cntr1:cntr2, qq)); 
        sigma_w = nanstd(raw.wPF(cntr1:cntr2, qq)); 
        U = nanmean(sqrt(raw.uPF(cntr1:cntr2, qq).^2+raw.vPF(cntr1:cntr2, qq).^2));
        P = nanmean(data{TableNum}(cntr1:cntr2, sensorInfo.P(siteInfoTableNum, 2)).*1000); %average Pressure [Pa]
        T = raw.fwT(cntr1:cntr2, qq)+273.15; %Temperature [K]
        TSon = raw.sonTs(cntr1:cntr2, qq)+273.15; %Sonic Temperature [k]
        
        %Calculate temperature structure parameter
        [~, Ct2(ii)] = ...
            crossStruct(T, T, 20, 'spatial', sigma_u, sigma_v, sigma_w, U); 
        if Ct2(ii)<0
            Ct2(ii) = nan;
        end
        
        [~, Ctson2(ii)] = ...
            crossStruct(TSon, TSon, 20, 'spatial', sigma_u, sigma_v, sigma_w, U); 
        if Ctson2(ii)<0
            Ctson2(ii) = nan;
        end

        
        if cntrIRGA<=MinInst
            %look for matching IRGA for Sonic height
            IRGAHeight = find(sensorInfo.irgaH2O(:, 3)==SonHeight);

            if ~isempty(IRGAHeight)
                %Abslute humidity [kg m^-3]
                rhov = raw.rhov(cntr1:cntr2, IRGAHeight)./1000;  

                %Calculate mean temperature for period
                T_mean = nanmean(T);
                if isnan(T_mean)
                    T_mean = nanmean(TSon);
                end

                %Specific humidity [kg/kg]
                q = (rhov./(P./(Rd*T_mean)));

                %Calculate temperature humidity cross structure parameter
                [~, Ctq(ii)] = ...
                    crossStruct(T, q, 20, 'spatial', sigma_u, sigma_v, sigma_w, U);
                
                [~, Ctsonq(ii)] = ...
                    crossStruct(TSon, q, 20, 'spatial', sigma_u, sigma_v, sigma_w, U);
                
                %Calculate humidity structure parameter
                [~, Cq2(ii)] = ...
                    crossStruct(q, q, 20, 'spatial', sigma_u, sigma_v, sigma_w, U);
                if Cq2(ii)<0
                    Cq2(ii) = nan;
                end

                %Calculate temperature humidity correlation
                r_tq(ii) = Ctq(ii)./sqrt(Ct2(ii).*Cq2(ii));
                
                %Calculate temperature humidity correlation
                r_tsonq(ii) = Ctsonq(ii)./sqrt(Ctson2(ii).*Cq2(ii));

                IRGAFlag = 1;
            else
                IRGAFlag = 0;
            end
        else
            IRGAFlag = 0;
        end
        %update counters
        cntr1 = cntr1 + pnts;
        cntr2 = cntr2 + pnts; 
    end

    %Save to output structure
    if qq==1
        if IRGAFlag
            %If IRGA Exists
            output.StructParam = [time, Ct2', Ctson2', Ctq', Ctsonq', Cq2', r_tq', r_tsonq'];
            output.StructParamHeader = {'Timestamp', ...
                [num2str(sensorInfo.Tson(qq, end)), 'm Ctfw2: K^2 m^-2/3'],...
                [num2str(sensorInfo.Tson(qq, end)), 'm CtSon2: K^2 m^-2/3'],...
                [num2str(sensorInfo.irgaH2O(cntrIRGA, end)), 'm Ctfwq: K kg kg^-1 m^-2/3'],...
                [num2str(sensorInfo.irgaH2O(cntrIRGA, end)), 'm CtSonq: K kg kg^-1 m^-2/3'],...
                [num2str(sensorInfo.irgaH2O(cntrIRGA, end)), 'm Cq2: kg^2 kg^-2 m^-2/3'],...
                [num2str(sensorInfo.irgaH2O(cntrIRGA, end)), 'm r_tfwq'],...
                [num2str(sensorInfo.irgaH2O(cntrIRGA, end)), 'm r_tSonq'];...
                [],num2str(sensorInfo.Tson(qq, end)),...
                num2str(sensorInfo.Tson(qq, end)),...
                num2str(sensorInfo.irgaH2O(cntrIRGA, end)),...
                num2str(sensorInfo.irgaH2O(cntrIRGA, end)),...
                num2str(sensorInfo.irgaH2O(cntrIRGA, end)),...
                num2str(sensorInfo.irgaH2O(cntrIRGA, end)),...
                num2str(sensorInfo.irgaH2O(cntrIRGA, end))};
            cntrIRGA = cntrIRGA+1;
        else
            %If only Sonics
            output.StructParam = [time, Ct2', CtSon2'];
            output.StructParamHeader = {'Timestamp', ...
                [num2str(sensorInfo.Tson(qq, end)), 'm Ctfw2: K^2 m^-2/3'],...
                [num2str(sensorInfo.Tson(qq, end)), 'm CtSon2: K^2 m^-2/3'];...
                [],num2str(sensorInfo.Tson(qq, end)),...
                num2str(sensorInfo.Tson(qq, end))};
        end
    else
        if IRGAFlag
            %IRGA Exists
            output.StructParam = [output.StructParam, Ct2', Ctson2', Ctq', Ctsonq', Cq2', r_tq',r_tsonq'];
            output.StructParamHeader = [output.StructParamHeader,...
                {[num2str(sensorInfo.Tson(qq, end)), 'm Ctfw2: K^2 m^-2/3'],...
                [num2str(sensorInfo.Tson(qq, end)), 'm CtSon2: K^2 m^-2/3'],...
                [num2str(sensorInfo.irgaH2O(cntrIRGA, end)), 'm Ctfwq: K kg kg^-1 m^-2/3'],...
                [num2str(sensorInfo.irgaH2O(cntrIRGA, end)), 'm CtSonq: K kg kg^-1 m^-2/3'],...
                [num2str(sensorInfo.irgaH2O(cntrIRGA, end)), 'm Cq2: kg^2 kg^-2 m^-2/3'],...
                [num2str(sensorInfo.irgaH2O(cntrIRGA, end)), 'm r_tfwq'],...
                [num2str(sensorInfo.irgaH2O(cntrIRGA, end)), 'm r_tSonq'];...
                num2str(sensorInfo.Tson(qq, end)),...
                num2str(sensorInfo.Tson(qq, end)),...
                num2str(sensorInfo.irgaH2O(cntrIRGA, end)),...
                num2str(sensorInfo.irgaH2O(cntrIRGA, end)),...
                num2str(sensorInfo.irgaH2O(cntrIRGA, end)),...
                num2str(sensorInfo.irgaH2O(cntrIRGA, end)),...
                num2str(sensorInfo.irgaH2O(cntrIRGA, end))}];
            cntrIRGA = cntrIRGA+1;
        else
            %Only Sonics
            output.StructParam = [output.StructParam, Ct2', Ctson2'];
            output.StructParamHeader = [output.StructParamHeader, ...
                {[num2str(sensorInfo.Tson(qq, end)), 'm Ctfw2: K^2 m^-2/3'],...
                [num2str(sensorInfo.Tson(qq, end)), 'm CtSon2: K^2 m^-2/3'];...
                num2str(sensorInfo.Tson(qq, end)),...
                num2str(sensorInfo.Tson(qq, end))}];
        end
    end

    clearvars C* r_*
    
end