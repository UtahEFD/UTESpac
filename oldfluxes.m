function [output, raw] = fluxes(data, rotatedSonicData, info, output, sensorInfo,tableNames)
% fluxes computes all turbulent statistics.  Spike and NaN flags are used from sonic and finewires to nan-out flagged
% statistics.
if isfield(sensorInfo,'u')
    try
        if info.detrendingFormat
            detrendFormat = 'linear';
        else
            detrendFormat = 'constant';
        end
        display(sprintf('Finding turbulent pertubations with %s detrending and computing fluxes.',detrendFormat))
        
        % find pressure and ref pot temp from sonic at pressure level
        if isfield(sensorInfo,'P')
            P0 = 100; % reference pressure in kPa
            tble = sensorInfo.P(1,1);  % Pressure table
            col = sensorInfo.P(1,2);   % Pressure Column
            P = data{1,tble}(:,col);  % Pressure in kPa or mbar
            Ptime = data{1,tble}(:,1); % pressure time stamps
            zRef = sensorInfo.P(1,3);  % Pressure Height
            
            % ensure that P is in kPa and not mbar
            if nanmedian(P) > 200
                P = P*0.1;
            end
            k = 0.286; %R/cp for air.  Arya Pg. 65
            tble = sensorInfo.Tson(sensorInfo.Tson(:,3)==zRef,1);  % Find sonic corresponding to pressure level
            col = sensorInfo.Tson(sensorInfo.Tson(:,3)==zRef,2);
            Tref = data{1,tble}(:,col); % reference sonic temperature in C or K
            if nanmedian(Tref)<200
                Tref = Tref + 273.15;  % ensure that Tref is in K
            end
            Ttime = data{1,tble}(:,1); % time stamps for reference T
            Pinterp = interp1(Ptime(~isnan(P)),P(~isnan(P)),Ttime,'linear','extrap');  % find interpolated pressure values at all sonT
            ThetaRef = Tref.*(P0./Pinterp).^k-273.15;  % in deg C
        else
            Tref = [];
            zRef = [];
            Pinterp = [];
            tble = sensorInfo.u(1,1);
            display('No barometer found. Theta and Thetav = NaN')
        end
        
        % find number of sonics
        numSonics = size(sensorInfo.u,1);
        
        % iterate through all sonics
        for i = 1:numSonics
            if i==1
                
                % find serial date column
                t = data{1,tble}(:,1);
                
                % find total number of averaging periods
                N = round((data{1,tble}(end,1)-data{1,tble}(1,1))/(info.avgPer/(24*60)));
                
                % find breakpoints for detrending.  size(bp,1) = size(N,1) + 1
                bp = round(linspace(0,size(data{1,tble},1),N+1));
                
                % initialize flux matrices
                H = nan(N,1);  % kinematic sensible heat flux
                Hlat = nan(N,1); % kinematic sensible, lateral heat flux
                tau = nan(N,1); % momentum flux
                tke = nan(N,1); % turbulent kinetic energy
                LHflux = nan(N,1); % latent heat flux
                CO2flux = nan(N,1); % CO2 flux
                derivedT = nan(N,1);  % matrix for derived temperatures
                specHum = nan(N,1);  % specific humidity in kg water vapor/ kg dry air
                L = nan(N,1); % Obukhov Length
                sigma = nan(N,1); % standard deviations (u, v, w)
                
                % raw flux matrices
                if info.saveRawConditionedData
                    raw.u = nan(bp(end),numSonics); % u
                    raw.v = nan(bp(end),numSonics); % v
                    raw.w = nan(bp(end),numSonics); % w
                    raw.uPF = nan(bp(end),numSonics);  % u planar fit and yaw corrected
                    raw.vPF = nan(bp(end),numSonics);  % v planar fit and yaw corrected
                    raw.wPF = nan(bp(end),numSonics);  % w planar fit and yaw corrected
                    raw.sonTs = nan(bp(end),numSonics); % Temperature from Sonic
                    raw.uPF_Prime = nan(bp(end),numSonics);  % u' planar fit and yaw corrected
                    raw.vPF_Prime = nan(bp(end),numSonics); % v' planar fit and yaw corrected
                    raw.wPF_Prime = nan(bp(end),numSonics); % w' planar fit and yaw corrected
                    raw.sonTsPrime = nan(bp(end),numSonics); % Ts' from sonic
                    raw.t = t; % serial time stamp
                    raw.z = nan(1,numSonics); % sonic heights
                    
                    % finewires
                    if isfield(sensorInfo,'fw')
                        raw.fwTh = nan(bp(end),numSonics); % theta from FW
                        raw.fwT = nan(bp(end),numSonics); % T from FW
                        faw.fwTPrime = nan(bp(end),numSonics); % T' from FW
                        raw.fwThPrime = nan(bp(end),numSonics); % theta' from finewire
                    end
                    
                    
                    % H2O
                    if isfield(sensorInfo,'irgaH2O')
                        numH2OSensors = size(sensorInfo.irgaH2O,1);
                        raw.H2O = nan(bp(end),numH2OSensors); % H2O
                        raw.H2OPrime = nan(bp(end),numH2OSensors); % H2O'
                    elseif isfield(sensorInfo,'KH2O')
                        numH2OSensors = size(sensorInfo.KH2O,1);
                        raw.H2O = nan(bp(end),numH2OSensors); % H2O
                        raw.H2OPrime = nan(bp(end),numH2OSensors); % H2O'
                    end
                    
                    % CO2
                    if isfield(sensorInfo,'irgaCO2')
                        numCO2Sensors = size(sensorInfo.irgaCO2,1);
                        raw.CO2 = nan(bp(end),numCO2Sensors); % CO2
                        raw.CO2Prime = nan(bp(end),numCO2Sensors); % CO2'
                    end
                    
                else
                    raw = [];
                end
                
                % initialize headers
                Hheader = cell(1);
                HlatHeader = cell(1);
                tauHeader = cell(1);
                tkeHeader = cell(1);
                LHfluxHeader = cell(1);
                CO2fluxHeader = cell(1);
                derivedTheader = cell(1);
                specHumHeader = cell(1);
                Lheader = cell(1);
                sigmaHeader = cell(1);
            end
            try
                % find sonic information
                tble = sensorInfo.u(i,1);
                sonHeight = sensorInfo.u(i,3);
                uCol = sensorInfo.u(sensorInfo.u(:,3)==sonHeight,2);
                vCol = sensorInfo.v(sensorInfo.v(:,3)==sonHeight,2);
                wCol = sensorInfo.w(sensorInfo.v(:,3)==sonHeight,2);
                TsCol = sensorInfo.Tson(sensorInfo.Tson(:,3)==sonHeight,2);
                
                % find sonic flag information (averaged values!)
                nanFlagTableName = [tableNames{tble},'NanFlag'];
                spikeFlagTableName = [tableNames{tble},'SpikeFlag'];
                uNanFlag = output.(nanFlagTableName)(:,uCol);
                uSpikeFlag = output.(spikeFlagTableName)(:,uCol);
                vNanFlag = output.(nanFlagTableName)(:,vCol);
                vSpikeFlag = output.(spikeFlagTableName)(:,vCol);
                wNanFlag = output.(nanFlagTableName)(:,wCol);
                wSpikeFlag = output.(spikeFlagTableName)(:,wCol);
                TsNanFlag = output.(nanFlagTableName)(:,TsCol);
                TsSpikeFlag = output.(nanFlagTableName)(:,TsCol);
                % check for diagnostic
                if isfield(sensorInfo,'sonDiagnostic')
                    sonicDiagnosticCol = sensorInfo.sonDiagnostic(sensorInfo.sonDiagnostic(:,3)==sonHeight,2);
                    if isempty(sonicDiagnosticCol)
                        sonicDiagnosticFlag = zeros(size(uNanFlag));
                    else
                        sonicDiagnosticFlag = output.(tableNames{tble})(:,sonicDiagnosticCol);
                        sonicDiagnosticFlag(sonicDiagnosticFlag < info.diagnosticTest.meanSonicDiagnosticLimit) = 0;
                    end
                end
                sonicDiagnosticFlag(isnan(sonicDiagnosticFlag)) = 0;
                
                % find unrotated sonic values
                u = data{1,tble}(:,uCol);
                v = data{1,tble}(:,vCol);
                w = data{1,tble}(:,wCol);
                unrotatedSonFlag = logical(wNanFlag+wSpikeFlag+sonicDiagnosticFlag); % total sonic flag for unrotated calculations
                
                % find rotated sonic columns
                uCol = strcmp(output.rotatedSonicHeader,strcat(num2str(sonHeight),'m:u'));
                vCol = strcmp(output.rotatedSonicHeader,strcat(num2str(sonHeight),'m:v'));
                wCol = strcmp(output.rotatedSonicHeader,strcat(num2str(sonHeight),'m:w'));
                
                % load rotated sonic values
                uPF = rotatedSonicData(:,uCol);
                vPF = rotatedSonicData(:,vCol);
                wPF = rotatedSonicData(:,wCol);
                rotatedSonFlag = logical(uNanFlag+uSpikeFlag+vNanFlag+vSpikeFlag+wNanFlag+wSpikeFlag+sonicDiagnosticFlag); % total sonic flag for rotated calculations
                
                % find sonic temperature
                Tson = data{1,tble}(:,TsCol);
                if nanmedian(Tson) > 250
                    Tson = Tson - 273.15;
                end
                TsonFlag = logical(TsNanFlag+TsSpikeFlag);
                
                % find slow RH and T measurment to get specific humidity
                if isfield(sensorInfo,'RH') && isfield(sensorInfo,'P')
                    RHtble = sensorInfo.RH(sensorInfo.RH(:,3)==sonHeight,1);
                    RHcol = sensorInfo.RH(sensorInfo.RH(:,3)==sonHeight,2);
                    Ttble = sensorInfo.T(sensorInfo.T(:,3)==sonHeight,1);
                    Tcol = sensorInfo.T(sensorInfo.T(:,3)==sonHeight,2);
                    
                    % if RH doesn't exist at sonic height, use first RH sensor
                    if isempty(RHcol)
                        RHtble = sensorInfo.RH(1,1);
                        RHcol = sensorInfo.RH(1,2);
                        RH = abs(data{1,RHtble}(:,RHcol));
                        Ttble = sensorInfo.T(1,1);
                        Tcol = sensorInfo.T(1,2);
                        T = abs(data{1,Ttble}(:,Tcol));
                    else
                        RH = abs(data{1,RHtble}(:,RHcol));
                        T = abs(data{1,Ttble}(:,Tcol));
                    end
                    
                    % ensure that RH is in percent and not a fraction
                    if nanmedian(RH) < 1
                        RH = RH*100;
                    end
                    
                    % ensure that T is in deg C and not Kelvin
                    if nanmedian(T) > 200
                        T = T - 273.15;
                    end
                    
                    % interp RH to the same time stamps as the sonic
                    RHinterp = interp1(data{1,RHtble}(~isnan(RH),1),RH(~isnan(RH)),data{1,tble}(:,1),'linear','extrap');
                    Tinterp = interp1(data{1,Ttble}(~isnan(T),1),T(~isnan(T)),data{1,tble}(:,1),'linear','extrap');
                    
                    % convert RH to q See: http://www.mech.utah.edu/~pardyjak/efd/pottempcalc.pdf
                    l_v = 2.5e6;    %Latent Heat of Vaporization (J/kg)
                    R_v = 461.5;    %Gas Constant for water (J-K/kg)
                    T_0 = 273.15;   %Reference Temperature
                    es0 = 6.11;     %Reference Vapor Pressure at 273.15 (hPa)
                    
                    %Find Saturated Vapor Pressure
                    es = es0.*exp((l_v/R_v).*(1/T_0-1./(Tinterp+273.15)));
                    e = RHinterp./100.*es; % vapor pressure
                    q = 0.622.*(e./(Pinterp.*10)); %(kg/kg) see Thermo by Cengel Ch. 14
                    % check output here: http://www.rotronic.com/humidity_measurement-feuchtemessung-mesure_de_l_humidite/humidity-calculator-feuchterechner-mr
                    
                    % store qRef
                    qRef = q;
                    
                    % average q down to averaging periods and store in output
                    qBar = simpleAvg([q t],info.avgPer);
                    
                    if i == 1
                        specHum(:,1:2) = fliplr(qBar);
                        specHumHeader{1} = 'time';
                        specHumHeader{2} = strcat(num2str(sonHeight),' m: q(g/g)');
                    else
                        specHum(:,end+1) = qBar(:,1);
                        specHumHeader{end+1} = strcat(num2str(sonHeight),' m: q(g/g)');
                    end
                elseif ~isfield(sensorInfo,'RH') && isfield(sensorInfo,'P') && isfield(sensorInfo,'irgaH2O')  % EC150
                    rhoH2O_tble = sensorInfo.irgaH2O(1,1);
                    rhoH2O_col = sensorInfo.irgaH2O(1,2);
                    rhoH2O = abs(data{1,rhoH2O_tble}(:,rhoH2O_col));
                    R = 287.058;  % Gas constant for air
                    q = (rhoH2O*1/1000)./(Pinterp*1000./(R*Tref)); % kg/kg  rhoH2O/rhoAir
                    qRef = q;
                elseif ~isfield(sensorInfo,'RH') && isfield(sensorInfo,'P') && isfield(sensorInfo,'KH2O') % KH2O
                    rhoH2O_tble = sensorInfo.KH2O(1,1);
                    rhoH2O_col = sensorInfo.KH2O(1,2);
                    rhoH2O = abs(data{1,rhoH2O_tble}(:,rhoH2O_col));
                    R = 287.058;  % Gas constant for air
                    q = (rhoH2O*1/1000)./(Pinterp*1000./(R*Tref)); % kg/kg  rhoH2O/rhoAir
                    qRef = q;
                else
                    q = [];
                    qRef = [];
                end
                
                % find fw
                if isfield(sensorInfo,'fw')
                    fwCol = sensorInfo.fw(sensorInfo.fw(:,3)==sonHeight,2);
                    fw = data{1,tble}(:,fwCol);
                    nanFlagTableName = [tableNames{tble},'NanFlag'];
                    spikeFlagTableName = [tableNames{tble},'SpikeFlag'];
                    fwNanFlag = output.(nanFlagTableName)(:,fwCol);
                    fwSpikeFlag = output.(spikeFlagTableName)(:,fwCol);
                    fwFlag = logical(fwNanFlag+fwSpikeFlag);
                else
                    fw = [];
                    fwFlag = [];
                end
                
                % find fw pot temp
                Gamma = 0.0098; % Lapse Rate. K/m
                if ~isempty(fw) && ~isempty(Tref)
                    thetaFw = fw - (Tref-273.15) + Gamma*(sonHeight - zRef) + ThetaRef;
                    temp = simpleAvg([thetaFw t],info.avgPer);
                    derivedT(:,end+1) = temp(:,1);
                    derivedTheader{end+1} = strcat(num2str(sonHeight),' m: theta_fw');
                else
                    thetaFw = [];
                end
                
                % find sonic pot temp
                if ~isempty(Tref)
                    thetaSon = Tson - (Tref-273.15) + Gamma*(sonHeight - zRef) + ThetaRef;
                    temp = simpleAvg([thetaSon t],info.avgPer);
                    derivedT(:,end+1) = temp(:,1);
                    derivedTheader{end+1} = strcat(num2str(sonHeight),' m: theta_s_son');
                else
                    thetaSon = [];
                end
                
                % find fw virt, pot temp
                if ~isempty(thetaFw) && ~isempty(q)
                    VthetaFw = thetaFw.*(1+0.61*q); % stull pg 7
                    temp = simpleAvg([VthetaFw t],info.avgPer);
                    derivedT(:,end+1) = temp(:,1);
                    derivedTheader{end+1} = strcat(num2str(sonHeight),' m: theta_v_fw');
                else
                    VthetaFw = [];
                end
                
                % find pot temp from sonic
                if ~isempty(thetaSon) && ~isempty(q)
                    thetaSonAir = thetaSon./(1+0.51*q); % stull pg 7
                    temp = simpleAvg([thetaSonAir t],info.avgPer);
                    derivedT(:,end+1) = temp(:,1);
                    derivedTheader{end+1} = strcat(num2str(sonHeight),' m: theta_son');
                else
                    thetaSonAir = [];
                end
                
                % find H2O and CO2 columns if they exist
                if isfield(sensorInfo,'irgaH2O')
                    irgaH2Ocol = sensorInfo.irgaH2O(sensorInfo.irgaH2O(:,3)==sonHeight,2);
                    irgaGasDiagCol = sensorInfo.irgaGasDiag(sensorInfo.irgaGasDiag(:,3)==sonHeight,2);
                    irgaH2OSigCol = sensorInfo.irgaH2OsigStrength(sensorInfo.irgaH2OsigStrength(:,3)==sonHeight,2);
                    irgaH2O = data{1,tble}(:,irgaH2Ocol);
                    nanFlagTableName = [tableNames{tble},'NanFlag'];
                    spikeFlagTableName = [tableNames{tble},'SpikeFlag'];
                    H2ONanFlag = output.(nanFlagTableName)(:,irgaH2Ocol);
                    H2OSpikeFlag = output.(spikeFlagTableName)(:,irgaH2Ocol);
                    H2OsigFlag = output.(tableNames{tble})(:,irgaH2OSigCol);
                    H2OsigFlag(H2OsigFlag > info.diagnosticTest.H2OminSignal) = 0;
                    H2OsigFlag(isnan(H2OsigFlag)) = 0;
                    gasDiagFlag = output.(tableNames{tble})(:,irgaGasDiagCol);
                    gasDiagFlag(gasDiagFlag < info.diagnosticTest.meanGasDiagnosticLimit) = 0;
                    gasDiagFlag(isnan(gasDiagFlag)) = 0;
                    H2OFlag = logical(H2ONanFlag+H2OSpikeFlag+H2OsigFlag+gasDiagFlag);
                else
                    irgaH2O = [];
                    H2OFlag = [];
                end
                if isfield(sensorInfo,'irgaCO2')
                    irgaCO2col = sensorInfo.irgaCO2(sensorInfo.irgaCO2(:,3)==sonHeight,2);
                    irgaCO2SigCol = sensorInfo.irgaCO2sigStrength(sensorInfo.irgaCO2sigStrength(:,3)==sonHeight,2);
                    irgaCO2 = data{1,tble}(:,irgaCO2col);
                    CO2sigFlag = output.(tableNames{tble})(:,irgaCO2SigCol);
                    CO2sigFlag(CO2sigFlag > info.diagnosticTest.CO2minSignal) = 0;
                    CO2sigFlag(isnan(CO2sigFlag)) = 0;
                    nanFlagTableName = [tableNames{tble},'NanFlag'];
                    spikeFlagTableName = [tableNames{tble},'SpikeFlag'];
                    CO2NanFlag = output.(nanFlagTableName)(:,irgaCO2col);
                    CO2SpikeFlag = output.(spikeFlagTableName)(:,irgaCO2col);
                    CO2Flag = logical(CO2NanFlag+CO2SpikeFlag+gasDiagFlag+CO2sigFlag);
                else
                    irgaCO2 = [];
                    CO2Flag = [];
                end
                if isfield(sensorInfo,'KH2O')
                    KH2Ocol = sensorInfo.KH2O(sensorInfo.KH2O(:,3)==sonHeight,2);
                    KH2O = data{1,tble}(:,KH2Ocol);
                    nanFlagTableName = [tableNames{tble},'NanFlag'];
                    spikeFlagTableName = [tableNames{tble},'SpikeFlag'];
                    H2ONanFlag = output.(nanFlagTableName)(:,KH2Ocol);
                    H2OSpikeFlag = output.(spikeFlagTableName)(:,KH2Ocol);
                    H2OFlag = logical(H2ONanFlag+H2OSpikeFlag);
                elseif isfield(sensorInfo,'irgaH2O') % don't delete H2O flag if irga exists!
                    KH2O = [];
                else
                    KH2O = [];
                    H2OFlag = [];
                end
                
                % iterate though all time stamps and find fluxes
                for j = 1:N
                    if i == 1
                        % place time stamp in column 1
                        H(j,1) = t(bp(j+1));
                        Hlat(j,1) = t(bp(j+1));
                        tau(j,1) = t(bp(j+1));
                        tke(j,1) = t(bp(j+1));
                        LHflux(j,1) = t(bp(j+1));
                        CO2flux(j,1) = t(bp(j+1));
                        derivedT(j,1) = t(bp(j+1));
                        L(j,1) = t(bp(j+1));
                        sigma(j,1) = t(bp(j+1));
                        if j == 1
                            Hheader{1} = 'time';
                            HlatHeader{1} = 'time';
                            tauHeader{1} = 'time';
                            tkeHeader{1} = 'time';
                            derivedTheader{1} = 'time';
                            LHfluxHeader{1} = 'time';
                            CO2fluxHeader{1} = 'time';
                            Lheader{1} = 'time';
                            sigmaHeader{1} = 'time';
                        end
                    end
                    
                    % find air density
                    Rd = 287.058;  % Gas constant for air
                    Rv = 461.495;  % Gas constant for water vapor
                    
                    % find moist air density with qRef and Barometer
                    if ~isempty(Pinterp) && ~isempty(qRef) 
                       
                        T = nanmedian(Tref(bp(j)+1:bp(j+1)));  % median temperature (K)
                        q = nanmedian(qRef(bp(j)+1:bp(j+1)));  % median specific humidity (kg_v/kg_a)
                        Pv = q*P/0.622;  % partial pressure of water vapor (Pa)
                        Pd = P-Pv;  % partial pressure of dry air
                        rho = Pd/(Rd*T)+Pv/(Rv*T); % density of moist air (kg/m^3)
                        rhod = Pd/(Rd*T); % density of dry air (kg/m^3)
                        rhov = Pv/(Rv*T); % density of water vapor (kg/m^3)
                        
                    % find moist air density with qRef only
                    elseif ~isempty(qRef)  
                        P = info.Pref*1000; % reference pressure defined in INFORMATION Block (Pa)
                        T = nanmedian(Tref(bp(j)+1:bp(j+1)));  % median temperature (K)
                        q = nanmedian(qRef(bp(j)+1:bp(j+1)));  % median specific humidity (kg_v/kg_a)
                        Pv = q*P/0.622;  % partial pressure of water vapor (Pa)
                        Pd = P-Pv;  % partial pressure of dry air
                        rho = Pd/(Rd*T)+Pv/(Rv*T);
                        rhod = Pd/(Rd*T); % density of dry air (kg/m^3)
                        rhov = Pv/(Rv*T); % density of water vapor (kg/m^3)
                        
                    % find air density, ignore moisture    
                    elseif ~isempty(Pinterp) 
                        P = nanmedian(Pinterp(bp(j)+1:bp(j+1)))*1000; % median air pressure (Pa)
                        T = nanmedian(Tref(bp(j)+1:bp(j+1)));  % median temperature (K)
                        rho = P/(Rd*T);
                        
                    % if neither barometer, nor qRef exist  
                    else 
                        P = info.Pref*1000; % median air pressure (Pa)
                        T = nanmedian(Tref(bp(j)+1:bp(j+1)));  % median temperature (K)
                        rho = P/(Rd*T);
                    end
                    
                    % place rho and cp in columns 2 and 3 of H
                    if sonHeight == zRef  % zRef is height of barometer!
                        Hheader{2} = 'rho';
                        H(j,2) = rho;
                        
                        % place cp in column 3 of H
                        Hheader{3} = 'cp';
                        if ~isempty(qRef)
                            H(j,3) = 1004.67*(1+0.84*nanmedian(qRef(bp(j)+1:bp(j+1)))); % pg. 640 of Stull
                        else
                            H(j,3) = 1005;
                        end
                    end
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    % find unrotated sonic pertubations
                    uP = nandetrend(u(bp(j)+1:bp(j+1)),detrendFormat);
                    vP = nandetrend(v(bp(j)+1:bp(j+1)),detrendFormat);
                    wP = nandetrend(w(bp(j)+1:bp(j+1)),detrendFormat);
                    TsonP = nandetrend(Tson(bp(j)+1:bp(j+1)),detrendFormat);
                    
                    % find rotated sonic pertubations
                    uPF_P = nandetrend(uPF(bp(j)+1:bp(j+1)),detrendFormat);
                    vPF_P = nandetrend(vPF(bp(j)+1:bp(j+1)),detrendFormat);
                    wPF_P = nandetrend(wPF(bp(j)+1:bp(j+1)),detrendFormat);
                    
                    % store rotated data
                    if info.saveRawConditionedData
                        raw.u(bp(j)+1:bp(j+1),i) = u(bp(j)+1:bp(j+1));
                        raw.v(bp(j)+1:bp(j+1),i) = v(bp(j)+1:bp(j+1));
                        raw.w(bp(j)+1:bp(j+1),i) = w(bp(j)+1:bp(j+1));
                        raw.uPF(bp(j)+1:bp(j+1),i) = uPF(bp(j)+1:bp(j+1));
                        raw.vPF(bp(j)+1:bp(j+1),i) = vPF(bp(j)+1:bp(j+1));
                        raw.wPF(bp(j)+1:bp(j+1),i) = wPF(bp(j)+1:bp(j+1));
                        raw.uPF_Prime(bp(j)+1:bp(j+1),i) = uPF_P;
                        raw.vPF_Prime(bp(j)+1:bp(j+1),i) = vPF_P;
                        raw.wPF_Prime(bp(j)+1:bp(j+1),i) = wPF_P;
                        raw.z(i) = sonHeight;
                    end
                    
                    % find standard deviations of wind vector, sonic and finwire temperature
                    if ~isempty(fw) % check for fw first
                        numSigmaVariables = 9;
                        
                        sigma(j,10+(i-1)*numSigmaVariables) = nanstd(fw(bp(j)+1:bp(j+1))); sigmaHeader{10+(i-1)*numSigmaVariables} = strcat(num2str(sonHeight),'m :sigma_TFW');
                        if fwFlag(j); sigma(j,10+(i-1)*numSigmaVariables) = nan; end;
                        
                    else
                        numSigmaVariables = 8;
                    end
                    
                    sigma(j,2+(i-1)*numSigmaVariables) = nanstd(u(bp(j)+1:bp(j+1))); sigmaHeader{2+(i-1)*numSigmaVariables} = strcat(num2str(sonHeight),'m :sigma_u');
                    if rotatedSonFlag(j); sigma(j,2+(i-1)*numSigmaVariables) = nan; end;
                    
                    sigma(j,3+(i-1)*numSigmaVariables) = nanstd(v(bp(j)+1:bp(j+1))); sigmaHeader{3+(i-1)*numSigmaVariables} = strcat(num2str(sonHeight),'m :sigma_v');
                    if rotatedSonFlag(j); sigma(j,3+(i-1)*numSigmaVariables) = nan; end;
                    
                    sigma(j,4+(i-1)*numSigmaVariables) = nanstd(w(bp(j)+1:bp(j+1))); sigmaHeader{4+(i-1)*numSigmaVariables} = strcat(num2str(sonHeight),'m :sigma_w');
                    if unrotatedSonFlag(j); sigma(j,4+(i-1)*numSigmaVariables) = nan; end;
                    
                    sigma(j,5+(i-1)*numSigmaVariables) = nanstd(uPF(bp(j)+1:bp(j+1))); sigmaHeader{5+(i-1)*numSigmaVariables} = strcat(num2str(sonHeight),'m :sigma_uPF');
                    if rotatedSonFlag(j); sigma(j,5+(i-1)*numSigmaVariables) = nan; end;
                    
                    sigma(j,6+(i-1)*numSigmaVariables) = nanstd(vPF(bp(j)+1:bp(j+1))); sigmaHeader{6+(i-1)*numSigmaVariables} = strcat(num2str(sonHeight),'m :sigma_vPF');
                    if rotatedSonFlag(j); sigma(j,6+(i-1)*numSigmaVariables) = nan; end;
                    
                    sigma(j,7+(i-1)*numSigmaVariables) = nanstd(wPF(bp(j)+1:bp(j+1))); sigmaHeader{7+(i-1)*numSigmaVariables} = strcat(num2str(sonHeight),'m :sigma_wPF');
                    if rotatedSonFlag(j); sigma(j,7+(i-1)*numSigmaVariables) = nan; end;
                    
                    sigma(j,8+(i-1)*numSigmaVariables) = nanstd(Tson(bp(j)+1:bp(j+1))); sigmaHeader{8+(i-1)*numSigmaVariables} = strcat(num2str(sonHeight),'m :sigma_Tson');
                    if TsonFlag(j); sigma(j,8+(i-1)*numSigmaVariables) = nan; end;
                    
                    sigma(j,9+(i-1)*numSigmaVariables) = nanmean(wPF_P.*TsonP.^2); sigmaHeader{9+(i-1)*numSigmaVariables} = strcat(num2str(sonHeight),'m :wPFP_TsonP_TsonP');
                    if TsonFlag(j) || rotatedSonFlag(j); sigma(j,9+(i-1)*numSigmaVariables) = nan; end;
                    
                    % find rotated and unrated momentum flux
                    tau(j,2+(i-1)*3) = sqrt(nanmean(uP.*wP)^2+nanmean(vP.*wP)^2); tauHeader{2+(i-1)*3} = strcat(num2str(sonHeight),'m :sqrt(u''w''^2+v''w''^2)');
                    if rotatedSonFlag(j); tau(j,2+(i-1)*3) = nan; end;
                    
                    tau(j,3+(i-1)*3) = sqrt(nanmean(uPF_P.*wPF_P)^2+nanmean(vPF_P.*wPF_P)^2); tauHeader{3+(i-1)*3} = strcat(num2str(sonHeight),'m :sqrt(uPF''wPF''^2+vPF''wPF''^2)');
                    if rotatedSonFlag(j); tau(j,3+(i-1)*3) = nan; end;
                    
                    tau(j,4+(i-1)*3) = nanmean(uPF_P.*wPF_P); tauHeader{4+(i-1)*3} = strcat(num2str(sonHeight),'m :uPF''wPF''');
                    if rotatedSonFlag(j); tau(j,4+(i-1)*3) = nan; end;
                    
                    % find tke
                    tke(j,1+i) = 1/2*(nanmean(uP.^2)+nanmean(vP.^2)+nanmean(wP.^2)); tkeHeader{1+i} = strcat(num2str(sonHeight),'m :0.5(u''^2+v''^2+w''^2)');
                    if rotatedSonFlag(j); tke(j,1+i) = nan; end;
                    
                    % find Ts'w' and Ts'wPF' from sonic
                    
                    H(j,4+(i-1)*12) = nanmean(wP.*TsonP); Hheader{4+(i-1)*12} = strcat(num2str(sonHeight),'m son:Ts''w''');
                    if unrotatedSonFlag(j)||TsonFlag(j); H(j,4+(i-1)*12) = nan; end;
                    
                    H(j,5+(i-1)*12) = nanmean(wPF_P.*TsonP);Hheader{5+(i-1)*12} = strcat(num2str(sonHeight),'m son:Ts''wPF''');
                    if rotatedSonFlag(j)||TsonFlag(j); H(j,5+(i-1)*12) = nan; end;
                    
                    Hlat(j,2+(i-1)*12) = nanmean(uP.*TsonP);HlatHeader{2+(i-1)*12} = strcat(num2str(sonHeight),'m son:Ts''u''');
                    if rotatedSonFlag(j)||TsonFlag(j); Hlat(j,2+(i-1)*12) = nan; end;
                    
                    Hlat(j,3+(i-1)*12) = nanmean(uPF_P.*TsonP);HlatHeader{3+(i-1)*12} = strcat(num2str(sonHeight),'m son:Ts''uPF''');
                    if rotatedSonFlag(j)||TsonFlag(j); Hlat(j,3+(i-1)*12) = nan; end;
                    
                    Hlat(j,4+(i-1)*12) = nanmean(vP.*TsonP);HlatHeader{4+(i-1)*12} = strcat(num2str(sonHeight),'m son:Ts''v''');
                    if rotatedSonFlag(j)||TsonFlag(j); Hlat(j,4+(i-1)*12) = nan; end;
                    
                    Hlat(j,5+(i-1)*12) = nanmean(vPF_P.*TsonP);HlatHeader{5+(i-1)*12} = strcat(num2str(sonHeight),'m son:Ts''vPF''');
                    if rotatedSonFlag(j)||TsonFlag(j); Hlat(j,5+(i-1)*12) = nan; end;
                    
                    % find obukhov length
                    kappa = 0.4;
                    g = 9.81;
                    T0_L = nanmedian(Tson(bp(j)+1:bp(j+1))) + 273.15; % Tson in K
                    uStarCubed = tau(j,3+(i-1)*3).^(3/2); % sqrt(uPF'*wPF')
                    H0_L = H(j,5+(i-1)*12); % wPF'.*Tson'
                    
                    L(j,2+(i-1)) = -uStarCubed/(kappa*g/T0_L*H0_L); Lheader{2+(i-1)} = strcat(num2str(sonHeight),'m L:sqrt(uPF''wPF'')^3/2*T_S/(k*g*wPF''Ts'')');
                    if rotatedSonFlag(j)||TsonFlag(j); L(j,2+(i-1)) = nan; end;
                    
                    % store thetaSonP in raw structure
                    if info.saveRawConditionedData
                        raw.sonTs(bp(j)+1:bp(j+1),i) = Tson(bp(j)+1:bp(j+1));
                        raw.sonTsPrime(bp(j)+1:bp(j+1),i) = TsonP;
                    end
                    
                    % find Th'w' and Th'wPF' from sonic
                    if ~isempty(thetaSonAir)
                        thetaSonAirP = nandetrend(thetaSonAir(bp(j)+1:bp(j+1)),detrendFormat);
                        
                        H(j,6+(i-1)*12) = nanmean(wP.*thetaSonAirP); Hheader{6+(i-1)*12} = strcat(num2str(sonHeight),'m son:Th''w''');
                        if unrotatedSonFlag(j)||TsonFlag(j); H(j,6+(i-1)*12) = nan; end;
                        
                        H(j,7+(i-1)*12) = nanmean(wPF_P.*thetaSonAirP);Hheader{7+(i-1)*12} = strcat(num2str(sonHeight),'m son:Th''wPF''');
                        if rotatedSonFlag(j)||TsonFlag(j); H(j,7+(i-1)*12) = nan; end;
                    end
                    
                    % find T'w' and T'wPF' from fw
                    if ~isempty(fw)
                        fwP = nandetrend(fw(bp(j)+1:bp(j+1)),detrendFormat);
                        
                        H(j,8+(i-1)*12) = nanmean(wP.*fwP); Hheader{8+(i-1)*12} = strcat(num2str(sonHeight),'m fw:T''w''');
                        if unrotatedSonFlag(j)||fwFlag(j); H(j,8+(i-1)*12) = nan; end;
                        
                        H(j,9+(i-1)*12) = nanmean(wPF_P.*fwP);Hheader{9+(i-1)*12} = strcat(num2str(sonHeight),'m fw:T''wPF''');
                        if rotatedSonFlag(j)||fwFlag(j); H(j,9+(i-1)*12) = nan; end;
                    end
                    
                    % find Th_s'w' and Th_s'wPF' from sonics
                    if ~isempty(thetaSon)
                        thetaSonP = nandetrend(thetaSon(bp(j)+1:bp(j+1)),detrendFormat);
                        
                        H(j,10+(i-1)*12) = nanmean(wP.*thetaSonP); Hheader{10+(i-1)*12} = strcat(num2str(sonHeight),'m son:Th_s''w''');
                        if unrotatedSonFlag(j)||TsonFlag(j); H(j,10+(i-1)*12) = nan; end;
                        
                        H(j,11+(i-1)*12) = nanmean(wPF_P.*thetaSonP);Hheader{11+(i-1)*12} = strcat(num2str(sonHeight),'m son:Th_s''wPF''');
                        if rotatedSonFlag(j)||TsonFlag(j); H(j,11+(i-1)*12) = nan; end;
                        
                        Hlat(j,6+(i-1)*12) = nanmean(uP.*thetaSonP); HlatHeader{6+(i-1)*12} = strcat(num2str(sonHeight),'m son:Th_s''u''');
                        if rotatedSonFlag(j)||TsonFlag(j); Hlat(j,6+(i-1)*12) = nan; end;
                        
                        Hlat(j,7+(i-1)*12) = nanmean(uPF_P.*thetaSonP);HlatHeader{7+(i-1)*12} = strcat(num2str(sonHeight),'m son:Th_s''uPF''');
                        if rotatedSonFlag(j)||TsonFlag(j); Hlat(j,7+(i-1)*12) = nan; end;
                        
                        Hlat(j,8+(i-1)*12) = nanmean(vP.*thetaSonP); HlatHeader{8+(i-1)*12} = strcat(num2str(sonHeight),'m son:Th_s''v''');
                        if rotatedSonFlag(j)||TsonFlag(j); Hlat(j,8+(i-1)*12) = nan; end;
                        
                        Hlat(j,9+(i-1)*12) = nanmean(vPF_P.*thetaSonP);HlatHeader{9+(i-1)*12} = strcat(num2str(sonHeight),'m son:Th_s''vPF''');
                        if rotatedSonFlag(j)||TsonFlag(j); Hlat(j,9+(i-1)*12) = nan; end;
                    end
                    
                    % find Th'w' and Th'wPF' from fw
                    if ~isempty(thetaFw)
                        thetaFwP = nandetrend(thetaFw(bp(j)+1:bp(j+1)),detrendFormat);
                        
                        H(j,12+(i-1)*12) = nanmean(wP.*thetaFwP); Hheader{12+(i-1)*12} = strcat(num2str(sonHeight),'m fw:Th''w''');
                        if unrotatedSonFlag(j)||fwFlag(j); H(j,12+(i-1)*12) = nan; end;
                        
                        H(j,13+(i-1)*12) = nanmean(wPF_P.*thetaFwP);Hheader{13+(i-1)*12} = strcat(num2str(sonHeight),'m fw:Th''wPF''');
                        if rotatedSonFlag(j)||fwFlag(j); H(j,13+(i-1)*12) = nan; end;
                        
                        Hlat(j,10+(i-1)*12) = nanmean(uP.*thetaFwP); HlatHeader{10+(i-1)*12} = strcat(num2str(sonHeight),'m fw:Th''u''');
                        if rotatedSonFlag(j)||fwFlag(j); Hlat(j,10+(i-1)*12) = nan; end;
                        
                        Hlat(j,11+(i-1)*12) = nanmean(uPF_P.*thetaFwP);HlatHeader{11+(i-1)*12} = strcat(num2str(sonHeight),'m fw:Th''uPF''');
                        if rotatedSonFlag(j)||fwFlag(j); Hlat(j,11+(i-1)*12) = nan; end;
                        
                        Hlat(j,12+(i-1)*12) = nanmean(vP.*thetaFwP); HlatHeader{12+(i-1)*12} = strcat(num2str(sonHeight),'m fw:Th''v''');
                        if rotatedSonFlag(j)||fwFlag(j); Hlat(j,12+(i-1)*12) = nan; end;
                        
                        Hlat(j,13+(i-1)*12) = nanmean(vPF_P.*thetaFwP);HlatHeader{13+(i-1)*12} = strcat(num2str(sonHeight),'m fw:Th''vPF''');
                        if rotatedSonFlag(j)||fwFlag(j); Hlat(j,13+(i-1)*12) = nan; end;
                        
                        % store fwthetaP in raw structure
                        if info.saveRawConditionedData
                            raw.fwThPrime(bp(j)+1:bp(j+1),i) = thetaFwP;
                            raw.fwTh(bp(j)+1:bp(j+1),i) = thetaFw(bp(j)+1:bp(j+1));
                            raw.fwT(bp(j)+1:bp(j+1),i) = fw(bp(j)+1:bp(j+1));
                        end
                    end
                    
                    % find VTh'w' and VTh'wPF' from fw
                    if ~isempty(VthetaFw)
                        VthetaFwP = nandetrend(VthetaFw(bp(j)+1:bp(j+1)),detrendFormat);
                        
                        H(j,14+(i-1)*12) = nanmean(wP.*VthetaFwP); Hheader{14+(i-1)*12} = strcat(num2str(sonHeight),'m fw:VTh''w''');
                        if unrotatedSonFlag(j)||fwFlag(j); H(j,14+(i-1)*12) = nan; end;
                        
                        H(j,15+(i-1)*12) = nanmean(wPF_P.*VthetaFwP);Hheader{15+(i-1)*12} = strcat(num2str(sonHeight),'m fw:VTh''wPF''');
                        if rotatedSonFlag(j)||fwFlag(j); H(j,15+(i-1)*12) = nan; end;
                    end
                    
                    
                    
                    
                    
                    
                    
                    % find latent heat flux if it exists at height
                    if ~isempty(irgaH2O) || ~isempty(KH2O)
                        
                        % find H2O measurement
                        if ~isempty(irgaH2O)
                            rho_v = irgaH2O./1000; % kg/m^3
                            H2Otype = 0;
                            rho_CO2 = irgaCO2./1000000; % kg/m^3
                            H2OsensorNumber = find(abs(sensorInfo.irgaH2O(:,3)-sonHeight)<0.2);  % find IRGA number to store raw pertubations
                            CO2sensorNumber = H2OsensorNumber;
                        else
                            rho_v = KH2O./MH2O; % kg/m^3
                            H2Otype = 1;
                            H2OsensorNumber = find(abs(sensorInfo.KH2O(:,3)-sonHeight)<0.2);  % find KH2O number to store raw pertubations
                        end
                        
                        % declare constants
                        MH2O = 18.0153;       %Molar Mass of H2O (g/mol)
                        MdryAir = 28.97;     %Molar Mass of Dry Air (g/mol)
                        MCO2 = 44.01;     %Molar Mass of Dry Air (g/mol)
                        
                        % find averaged variables
                        rho_v_avg = nanmean(rho_v(bp(j)+1:bp(j+1))); %(kg/m^3)
                        TrefAvg = nanmean(Tref(bp(j)+1:bp(j+1)));
                        kinSenFlux = H(j,5+(i-1)*12); % wPF_P*TsonP
                        if H(j,2) > 0.95 && H(j,2) < 1.3  % make sure rho is stored in column 2 of H.  If not, use rho = 1.05
                            rho_a_avg = H(j,2);
                        else
                            rho_a_avg = 1.05;
                        end
                        rho_d_avg = rho_a_avg - rho
                        
                        Lv = (2.501-0.00237*(TrefAvg-273.15))*10^3*MH2O; % Latent Heat of Vaporization in Moles (J/mol) Stoll P. 641
                        % Xv = (H2Oavg/MH2O)/((rho*1000-H2Oavg*MH2O)/MdryAir); %Molar Mixing Ratio (mol_H2O/mol_dry_air)
                        Xv = (rho_v_avg)/((rho_a_avg*1000-rho_v_avg*MH2O)/MdryAir); %Molar Mixing Ratio (mol_H2O/mol_dry_air)
                        
                        % find H2O pertubations in mol/m^3
                        H2Op = nandetrend(rho_v(bp(j)+1:bp(j+1)),detrendFormat); %(mol/m^3)
                        
                        % store H2OP in raw structure
                        if info.saveRawConditionedData
                            raw.H2OPrime(bp(j)+1:bp(j+1),H2OsensorNumber) = H2Op; % (mol/m^3)
                            raw.H2O(bp(j)+1:bp(j+1),H2OsensorNumber) = rho_v(bp(j)+1:bp(j+1)); % (mol/m^3)
                        end
                        
                        
                        E = nanmean(wP.*H2Op);     % find evaporation flux from unrotated data
                        ER = nanmean(wPF_P.*H2Op);  % find evaporation flux from rotated data
                        
                        LHflux(j,2+(i-1)*7) = Lv; LHfluxHeader{2+(i-1)*7} = strcat(num2str(sonHeight),'m Lv(J/mol)');
                        
                        LHflux(j,3+(i-1)*7) = E; LHfluxHeader{3+(i-1)*7} = strcat(num2str(sonHeight),'m w'':E(mol/m^2s)');
                        if unrotatedSonFlag(j)||H2OFlag(j); LHflux(j,3+(i-1)*7) = nan; end;
                        
                        LHflux(j,4+(i-1)*7) = ER; LHfluxHeader{4+(i-1)*7} = strcat(num2str(sonHeight),'m wPF'':E(mol/m^2s)');
                        if rotatedSonFlag(j)||H2OFlag(j); LHflux(j,4+(i-1)*7) = nan; end;
                        
                        LHflux(j,5+(i-1)*7) = Lv*(1+Xv)*(E+(rho_v_avg/TrefAvg)*kinSenFlux); LHfluxHeader{5+(i-1)*7} = strcat(num2str(sonHeight),'m WPL, w'' (W/m^2)');
                        if unrotatedSonFlag(j)||H2OFlag(j); LHflux(j,5+(i-1)*7) = nan; end;
                        
                        LHflux(j,6+(i-1)*7) = Lv*(1+Xv)*(ER+(rho_v_avg/TrefAvg)*kinSenFlux); LHfluxHeader{6+(i-1)*7} = strcat(num2str(sonHeight),'m WPL, wPF'' (W/m^2)');
                        if rotatedSonFlag(j)||H2OFlag(j); LHflux(j,6+(i-1)*7) = nan; end;
                        
                        if H2Otype %Peform O2 corrections for KH2O (E.C. by Marc Aubinet Pg. 104)
                            rhoDryAir = rho_a_avg*1000-rho_v_avg*(MH2O);
                            mixRatConvert = rhoDryAir/(MH2O);
                            k0 = -0.0045; %Extinction Coefficient for O2
                            kw = -0.159;  %Extinction Coefficient for H2O
                            Ck0 = 0.23*k0/kw;
                            correction = Ck0*(rho_v_avg*MH2O)/(TrefAvg)*kinSenFlux*mixRatConvert;
                            LHflux(j,7+(i-1)*7) = LHflux(j,5+(i-1)*7)+Lv*correction; LHfluxHeader{7+(i-1)*7} = strcat(num2str(sonHeight),'m WPL,O2,w'' (W/m^2)');
                            LHflux(j,8+(i-1)*7) = LHflux(j,6+(i-1)*7)+Lv*correction; LHfluxHeader{8+(i-1)*7} = strcat(num2str(sonHeight),'m WPL,O2,wPF'' (W/m^2)');
                        else
                            CO2p = nandetrend(rho_CO2(bp(j)+1:bp(j+1)),detrendFormat);
                            CO2avg = nanmean(rho_CO2(bp(j)+1:bp(j+1)));
                            latentFlux = LHflux(j,5+(i-1)*7);
                            Cbar = 1000*rho_a_avg./(MdryAir); % air concentration
                            
                            CO2flux(j,2+(i-1)*4) = nanmean(wP.*CO2p); CO2fluxHeader{2+(i-1)*4} = strcat(num2str(sonHeight),'m w'':CO2(mol/m^2s)');
                            if unrotatedSonFlag(j)||CO2Flag(j); CO2flux(j,2+(i-1)*4) = nan; end;
                            
                            CO2flux(j,3+(i-1)*4) = nanmean(wPF_P.*CO2p); CO2fluxHeader{3+(i-1)*4} = strcat(num2str(sonHeight),'m wPF'':CO2(mol/m^2s)');
                            if rotatedSonFlag(j)||CO2Flag(j); CO2flux(j,2+(i-1)*4) = nan; end;
                            
                            CO2flux(j,4+(i-1)*4) = CO2flux(j,2+(i-1)*4) + CO2avg.*(latentFlux/(Lv.*Cbar)+kinSenFlux/TrefAvg); CO2fluxHeader{4+(i-1)*4} = strcat(num2str(sonHeight),'m WPL,w'':CO2(mol/m^2s)');
                            
                            CO2flux(j,5+(i-1)*4) = CO2flux(j,3+(i-1)*4) + CO2avg.*(latentFlux/(Lv.*Cbar)+kinSenFlux/TrefAvg); CO2fluxHeader{5+(i-1)*4} = strcat(num2str(sonHeight),'m WPL,wPF'':CO2(mol/m^2s)');
                            
                            % store CO2P in raw structure
                            if info.saveRawConditionedData
                                raw.CO2Prime(bp(j)+1:bp(j+1),CO2sensorNumber) = CO2p;
                                raw.CO2(bp(j)+1:bp(j+1),CO2sensorNumber) = rho_CO2(bp(j)+1:bp(j+1));
                            end
                        end
                    end
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    %                     % find latent heat flux if it exists at height
                    %                     if ~isempty(irgaH2O) || ~isempty(KH2O)
                    %                         MH2O = 18.0153;       %Molar Mass of H2O (g/mol)
                    %                         MdryAir = 28.97;     %Molar Mass of Dry Air (kg/mol)
                    %                         MCO2 = 44.01;     %Molar Mass of Dry Air (g/mol)
                    %                         if ~isempty(irgaH2O)
                    %                             H2O = irgaH2O./MH2O; % mol/m^3
                    %                             H2Otype = 0;
                    %                             CO2 = 1000*irgaCO2./MCO2; % mol/m^3
                    %                             H2OsensorNumber = find(abs(sensorInfo.irgaH2O(:,3)-sonHeight)<0.2);  % find IRGA number to store raw pertubations
                    %                             CO2sensorNumber = H2OsensorNumber;
                    %                         else
                    %                             H2O = KH2O./MH2O; % mol/m^3
                    %                             H2Otype = 1;
                    %                             H2OsensorNumber = find(abs(sensorInfo.KH2O(:,3)-sonHeight)<0.2);  % find KH2O number to store raw pertubations
                    %                         end
                    %
                    %                         % find H2O pertubations in mol/m^3
                    %                         H2Op = nandetrend(H2O(bp(j)+1:bp(j+1)),detrendFormat); %(mol/m^3)
                    %
                    %                         % store H2OP in raw structure
                    %                         if info.saveRawConditionedData
                    %                             raw.H2OPrime(bp(j)+1:bp(j+1),H2OsensorNumber) = H2Op; % (mol/m^3)
                    %                             raw.H2O(bp(j)+1:bp(j+1),H2OsensorNumber) = H2O(bp(j)+1:bp(j+1)); % (mol/m^3)
                    %                         end
                    %
                    %                         % find avgH2O, Tref and rho
                    %                         H2Oavg = nanmean(H2O(bp(j)+1:bp(j+1))); %(mol/m^3)
                    %                         TrefAvg = nanmean(Tref(bp(j)+1:bp(j+1)));
                    %                         kinSenFlux = H(j,5+(i-1)*12); % wPF_P*TsonP
                    %                         if H(j,2) > 0.95 && H(j,2) < 1.3  % make sure rho is stored in column 2 of H.  If not, use rho = 1.05
                    %                             rho = H(j,2);
                    %                         else
                    %                             rho = 1.05;
                    %                         end
                    %
                    %                         Lv = (2.501-0.00237*(TrefAvg-273.15))*10^3*MH2O; % Latent Heat of Vaporization in Moles (J/mol) Stoll P. 641
                    %                         % Xv = (H2Oavg/MH2O)/((rho*1000-H2Oavg*MH2O)/MdryAir); %Molar Mixing Ratio (mol_H2O/mol_dry_air)
                    %                         Xv = (H2Oavg)/((rho*1000-H2Oavg*MH2O)/MdryAir); %Molar Mixing Ratio (mol_H2O/mol_dry_air)
                    %                         E = nanmean(wP.*H2Op);     % find evaporation flux from unrotated data
                    %                         ER = nanmean(wPF_P.*H2Op);  % find evaporation flux from rotated data
                    %
                    %                         LHflux(j,2+(i-1)*7) = Lv; LHfluxHeader{2+(i-1)*7} = strcat(num2str(sonHeight),'m Lv(J/mol)');
                    %
                    %                         LHflux(j,3+(i-1)*7) = E; LHfluxHeader{3+(i-1)*7} = strcat(num2str(sonHeight),'m w'':E(mol/m^2s)');
                    %                         if unrotatedSonFlag(j)||H2OFlag(j); LHflux(j,3+(i-1)*7) = nan; end;
                    %
                    %                         LHflux(j,4+(i-1)*7) = ER; LHfluxHeader{4+(i-1)*7} = strcat(num2str(sonHeight),'m wPF'':E(mol/m^2s)');
                    %                         if rotatedSonFlag(j)||H2OFlag(j); LHflux(j,4+(i-1)*7) = nan; end;
                    %
                    %                         LHflux(j,5+(i-1)*7) = Lv*(1+Xv)*(E+(H2Oavg/TrefAvg)*kinSenFlux); LHfluxHeader{5+(i-1)*7} = strcat(num2str(sonHeight),'m WPL, w'' (W/m^2)');
                    %                         if unrotatedSonFlag(j)||H2OFlag(j); LHflux(j,5+(i-1)*7) = nan; end;
                    %
                    %                         LHflux(j,6+(i-1)*7) = Lv*(1+Xv)*(ER+(H2Oavg/TrefAvg)*kinSenFlux); LHfluxHeader{6+(i-1)*7} = strcat(num2str(sonHeight),'m WPL, wPF'' (W/m^2)');
                    %                         if rotatedSonFlag(j)||H2OFlag(j); LHflux(j,6+(i-1)*7) = nan; end;
                    %
                    %                         if H2Otype %Peform O2 corrections for KH2O (E.C. by Marc Aubinet Pg. 104)
                    %                             rhoDryAir = rho*1000-H2Oavg*(MH2O);
                    %                             mixRatConvert = rhoDryAir/(MH2O);
                    %                             k0 = -0.0045; %Extinction Coefficient for O2
                    %                             kw = -0.159;  %Extinction Coefficient for H2O
                    %                             Ck0 = 0.23*k0/kw;
                    %                             correction = Ck0*(H2Oavg*MH2O)/(TrefAvg)*kinSenFlux*mixRatConvert;
                    %                             LHflux(j,7+(i-1)*7) = LHflux(j,5+(i-1)*7)+Lv*correction; LHfluxHeader{7+(i-1)*7} = strcat(num2str(sonHeight),'m WPL,O2,w'' (W/m^2)');
                    %                             LHflux(j,8+(i-1)*7) = LHflux(j,6+(i-1)*7)+Lv*correction; LHfluxHeader{8+(i-1)*7} = strcat(num2str(sonHeight),'m WPL,O2,wPF'' (W/m^2)');
                    %                         else
                    %                             CO2p = nandetrend(CO2(bp(j)+1:bp(j+1)),detrendFormat);
                    %                             CO2avg = nanmean(CO2(bp(j)+1:bp(j+1)));
                    %                             latentFlux = LHflux(j,5+(i-1)*7);
                    %                             Cbar = 1000*rho./(MdryAir); % air concentration
                    %
                    %                             CO2flux(j,2+(i-1)*4) = nanmean(wP.*CO2p); CO2fluxHeader{2+(i-1)*4} = strcat(num2str(sonHeight),'m w'':CO2(mol/m^2s)');
                    %                             if unrotatedSonFlag(j)||CO2Flag(j); CO2flux(j,2+(i-1)*4) = nan; end;
                    %
                    %                             CO2flux(j,3+(i-1)*4) = nanmean(wPF_P.*CO2p); CO2fluxHeader{3+(i-1)*4} = strcat(num2str(sonHeight),'m wPF'':CO2(mol/m^2s)');
                    %                             if rotatedSonFlag(j)||CO2Flag(j); CO2flux(j,2+(i-1)*4) = nan; end;
                    %
                    %                             CO2flux(j,4+(i-1)*4) = CO2flux(j,2+(i-1)*4) + CO2avg.*(latentFlux/(Lv.*Cbar)+kinSenFlux/TrefAvg); CO2fluxHeader{4+(i-1)*4} = strcat(num2str(sonHeight),'m WPL,w'':CO2(mol/m^2s)');
                    %
                    %                             CO2flux(j,5+(i-1)*4) = CO2flux(j,3+(i-1)*4) + CO2avg.*(latentFlux/(Lv.*Cbar)+kinSenFlux/TrefAvg); CO2fluxHeader{5+(i-1)*4} = strcat(num2str(sonHeight),'m WPL,wPF'':CO2(mol/m^2s)');
                    %
                    %                             % store CO2P in raw structure
                    %                             if info.saveRawConditionedData
                    %                                 raw.CO2Prime(bp(j)+1:bp(j+1),CO2sensorNumber) = CO2p;
                    %                                 raw.CO2(bp(j)+1:bp(j+1),CO2sensorNumber) = CO2(bp(j)+1:bp(j+1));
                    %                             end
                    %                         end
                    %                     end
                    
                    
                    
                    
                    
                    
                    
                    
                    
                end
            catch err
                message = strcat(err.message,'@ line',num2str(err.stack.line),' Problem with sonic at ',num2str(sonHeight),' m will be skipped');
                warning(message)
                if isempty(output.warnings{1})
                    output.warnings{1,1} = message;
                else output.warnings{end+1,1} = message;
                end
            end
        end
        % store outputs.  Use flag to delete empty columns!
        flag = logical(any(H,1)+isnan(H(1,:)));
        output.H = H(:,flag);
        output.Hheader = Hheader(flag);
        
        flag = logical(any(Hlat,1)+isnan(Hlat(1,:)));
        output.Hlat = Hlat(:,flag);
        output.HlatHeader = HlatHeader(flag);
        
        flag = logical(any(tau,1)+isnan(tau(1,:)));
        output.tau = tau(:,flag);
        output.tauHeader = tauHeader(flag);
        
        flag = logical(any(tke,1)+isnan(tke(1,:)));
        output.tke = tke(:,flag);
        output.tkeHeader = tkeHeader(flag);
        
        flag = logical(any(sigma,1)+isnan(sigma(1,:)));
        output.sigma = sigma(:,flag);
        output.sigmaHeader = sigmaHeader(flag);
        
        flag = logical(any(L,1)+isnan(L(1,:)));
        output.L = L(:,flag);
        output.Lheader = Lheader(flag);
        
        if size(specHum,2) > 1;
            flag = logical(any(specHum,1)+isnan(specHum(1,:)));
            output.specificHum = specHum(:,flag);
            output.specificHumHeader = specHumHeader(flag);
        end
        
        if size(derivedT,2) > 1;
            flag = logical(any(derivedT,1)+isnan(derivedT(1,:)));
            output.derivedT = derivedT(:,flag);
            output.derivedTheader = derivedTheader(flag);
        end
        
        if size(LHflux,2) > 1;
            flag = logical(any(LHflux,1)+isnan(LHflux(1,:)));
            output.LHflux = LHflux(:,flag);
            output.LHfluxHeader = LHfluxHeader(flag);
        end
        
        if size(CO2flux,2) > 1;
            flag = logical(any(CO2flux,1)+isnan(CO2flux(1,:)));
            output.CO2flux = CO2flux(:,flag);
            output.CO2fluxHeader = CO2fluxHeader(flag);
        end
    catch err
        message = strcat(err.message,'@ line ',num2str(err.stack.line),' UNABLE TO FIND FLUXES AT All HEIGHTS');
        warning(message)
        raw = [];
        if isempty(output.warnings{1})
            output.warnings{1,1} = message;
        else output.warnings{end+1,1} = message;
        end
    end
end
pause(4)
clc
end