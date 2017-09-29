function build_data(sectionType,memberType,startSpecimen)

if nargin == 0
    build_data('CCFT','C+PBC')
    build_data('RCFT','C+PBC')
    build_data('SRC' ,'C+PBC')
    build_data('CCFT','Beams')
    build_data('RCFT','Beams')
    build_data('CCFT','Other')
    build_data('RCFT','Other')
    return
end

if nargin < 3
    startSpecimen = 1;
end

% Units
dbUnits = 'US';
dbUnitSystem = unitSystem(dbUnits);

% Which strengths to compute
compute_AISC2016     = false;
compute_PSD          = false;
compute_Analysis_PfD = false;

%% Read and Adjust Data
filename = sprintf('%s_%s.csv',sectionType,memberType);
[~,S] = csvread2(filename);
numData = length(S.Author);

data(numData) = struct;

% Trim whitespace from text fields
S.Author = strtrim(S.Author);
S.Specimen = strtrim(S.Specimen);

% Make sure things that are supposed to be cell arrays are cell arrays 
% even if they are blank
if isnumeric(S.Year)
    S.Year = cellfun(@int2str,num2cell(S.Year),'UniformOutput',false);
end
if isnumeric(S.Tags)
    S.Tags = cell(size(S.Tags));
end
if isnumeric(S.Notes)
    S.Notes = cell(size(S.Notes));
end
if strcmp(memberType,'C+PBC')
    if isnumeric(S.eb_units)
        S.eb_units = cell(size(S.eb_units));
    end
end
if strcmp(memberType,'Other')
    if isnumeric(S.M1exp_units)
        S.M1exp_units = cell(size(S.M1exp_units));
    end
end


%% Interpret Data
for i = [startSpecimen:numData 1:(startSpecimen-1)]
    
    try
        % General Information   
        data(i).Author      = S.Author{i};
        data(i).Year        = S.Year{i};
        data(i).Reference   = sprintf('%s %s',S.Author{i},S.Year{i});
        data(i).Specimen    = S.Specimen{i};
        data(i).Tags        = S.Tags{i};
        data(i).Notes       = S.Notes{i};


        % Steel Strength
        data(i).Fy     = unitConvert('stress',S.Fy(i),S.Fy_units{i},dbUnitSystem);
        data(i).Fu     = unitConvert('stress',S.Fu(i),S.Fu_units{i},dbUnitSystem);

        % Concrete Strength
        switch lower(S.fc_type{i})
            case {'cube','cube/100mm','cube/150mm','cube/200mm'}
                data(i).fc = 0.71*unitConvert('stress',S.fc(i),S.fc_units{i},dbUnitSystem);
            case {'','cylinder','cylinder/100mm','cylinder/150mm'}
                data(i).fc = unitConvert('stress',S.fc(i),S.fc_units{i},dbUnitSystem);
            otherwise
                error('Unknown fc type: %s',S.fc_type{i});
        end

        % Data Specific to Section Type
        switch sectionType
            case 'CCFT'
                data(i).D = unitConvert('length',S.D(i),S.D_units{i},dbUnitSystem);
                data(i).t = unitConvert('length',S.t(i),S.t_units{i},dbUnitSystem);

            case 'RCFT'
                data(i).B = unitConvert('length',S.B(i),S.B_units{i},dbUnitSystem);
                data(i).H = unitConvert('length',S.H(i),S.H_units{i},dbUnitSystem);
                data(i).t = unitConvert('length',S.t(i),S.t_units{i},dbUnitSystem);
                data(i).SteelTubeType = S.SteelTubeType{i};

            case 'SRC'
                % Gross Section Dimensions
                data(i).H = unitConvert('length',S.H(i),S.H_units{i},dbUnitSystem);
                data(i).B = unitConvert('length',S.B(i),S.B_units{i},dbUnitSystem);

                % Steel Shape Data
                data(i).ShapeName = S.ShapeName{i};
                data(i).d  = unitConvert('length',S.d(i) ,S.d_units{i} ,dbUnitSystem);
                data(i).tw = unitConvert('length',S.tw(i),S.tw_units{i},dbUnitSystem);
                data(i).bf = unitConvert('length',S.bf(i),S.bf_units{i},dbUnitSystem);
                data(i).tf = unitConvert('length',S.tf(i),S.tf_units{i},dbUnitSystem);

                % Longitudinal Reinforcement
                switch lower(S.config_longitudinal{i})
                    case 'none'
                        data(i).config = 'none';
                        data(i).db     = 0;
                        data(i).Fylr   = 0;
                    case '2x-2y'
                        data(i).config = S.config_longitudinal{i};
                        db = str2double(S.db{i});
                        if ~isempty(db) && ~isnan(db)
                            data(i).db = unitConvert('length',db,S.db_units{i},dbUnitSystem);
                        else
                            switch S.db{i}
                                case '#5'
                                    data(i).db = unitConvert('length',5/8,'in',dbUnitSystem);
                                case '#7'
                                    data(i).db = unitConvert('length',7/8,'in',dbUnitSystem);
                                otherwise
                                    error('Unknown db: %s',S.db{i});
                            end
                        end
                        data(i).Fylr    = unitConvert('stress',S.Fylr(i),S.Fylr_units{i},dbUnitSystem);
                    otherwise

                end

                % Lateral Reinforcement
                switch lower(S.config_lateral{i})
                    case 'none'
                        data(i).dbTies = 0;
                        data(i).s      = 0;
                        data(i).Fytr   = 0;
                        data(i).cover  = 0;

                    case 'ties'
                        dbTies = str2double(S.dbTies{i});
                        if ~isempty(dbTies) && ~isnan(dbTies)
                            data(i).dbTies = unitConvert('length',dbTies,S.dbTies_units{i},dbUnitSystem);
                        else
                            switch S.dbTies{i}
                                case '#3'
                                    data(i).dbTies = unitConvert('length',3/8,'in',dbUnitSystem);
                                otherwise
                                    error('Unknown dbTies: %s',S.dbTies{i});
                            end
                        end
                        data(i).s       = unitConvert('length',S.s(i),S.s_units{i},dbUnitSystem);
                        data(i).Fytr    = unitConvert('stress',S.Fytr(i),S.Fytr_units{i},dbUnitSystem);
                        data(i).cover   = unitConvert('length',S.cover(i),S.cover_units{i},dbUnitSystem) - data(i).dbTies - data(i).db/2;
                    otherwise

                end

            otherwise
                error('Unknown section type %s',sectionType)
        end

        % Bending Axis
        if strcmp(sectionType,'CCFT')
            data(i).axis = 'weak';
        else
            data(i).axis = S.BendingAxis{i};
        end

        % Load Deformation
        data(i).LoadDeformation = S.LoadDeformation{i};

        % Data Specific to Member Type
        switch memberType
            case 'C+PBC'
                switch S.L_units{i}
                    case 'ratio_D'
                        data(i).L = S.L(i)*data(i).D;
                    case 'ratio_H'
                        data(i).L = S.L(i)*data(i).H;
                    otherwise
                        data(i).L = unitConvert('length',S.L(i),S.L_units{i},dbUnitSystem);
                end

                data(i).et = unitConvert('length',S.et(i),S.et_units{i},dbUnitSystem);

                if isnan(S.eb(i)) || isempty(S.eb(i))
                    data(i).eb = data(i).et;
                else
                    data(i).eb = unitConvert('length',S.eb(i),S.eb_units{i},dbUnitSystem);
                end

                data(i).isColumn = (data(i).et == 0 && data(i).eb == 0);

                data(i).Pexp = unitConvert('force',S.Pexp(i),S.Pexp_units{i},dbUnitSystem);

            case 'Beams'
                data(i).Mexp = unitConvert('moment',S.Mexp(i),S.Mexp_units{i},dbUnitSystem);

            case 'Other'
                data(i).KL      = unitConvert('length',S.KL(i),S.KL_units{i},dbUnitSystem);
                data(i).Pexp    = unitConvert( 'force',S.Pexp(i),S.Pexp_units{i},dbUnitSystem);
                data(i).M1exp   = unitConvert('moment',S.M1exp(i),S.M1exp_units{i},dbUnitSystem);
                data(i).M2exp   = unitConvert('moment',S.M2exp(i),S.M2exp_units{i},dbUnitSystem);

            otherwise
                error('Unknown member type %s',memberType);
        end

        % Computed Section Data
        switch sectionType
            case 'CCFT'
                section = CCFT(data(i).D,data(i).t,data(i).Fy,data(i).fc,dbUnits);

                data(i).compactness = section.slendernessInCompression;
                data(i).rho_s       = section.As/section.Ag;
                data(i).rho_sr      = 0;

            case 'RCFT'
                section = RCFT(data(i).H,data(i).B,data(i).t,data(i).Fy,data(i).fc,dbUnits);

                data(i).compactness = section.slendernessInCompression;
                data(i).rho_s       = section.As/section.Ag;
                data(i).rho_sr      = 0;

            case 'SRC'
                section = SRC(data(i).d,data(i).tw,data(i).bf,data(i).tf,data(i).Fy,...
                    data(i).H,data(i).B,data(i).fc,...
                    data(i).db,data(i).config,data(i).Fylr,...
                    data(i).dbTies,data(i).s,data(i).Fytr,data(i).cover,dbUnits);     

                data(i).compactness = 'Not Applicable';
                data(i).rho_s       = section.As/section.Ag;
                data(i).rho_sr      = section.Asr/section.Ag;

            otherwise
                error('Unknown section type: %s',sectionType)
        end

        % Member Data
        switch memberType
            case 'C+PBC'

                % Assign Length
                section.Lx = data(i).L;
                section.Ly = data(i).L;
                section.Kx = 1;
                section.Ky = 1;
                
                % Store Basic Data
                data(i).depth = section.depth(data(i).axis);            

                % Get Relevent Data
                e = max(abs([data(i).et data(i).eb]));
                if data(i).isColumn
                    beta = 1;
                else
                    if abs(data(i).et) > abs(data(i).eb)
                        beta = data(i).eb/data(i).et;
                    else
                        beta = data(i).et/data(i).eb;
                    end
                end

                % AISC 2016
                if compute_AISC2016
                    section.option_EI = 'AISC2016';
                    Pno       = section.Pnco;
                    EIelastic = 0.64*section.EIeff(data(i).axis);
                    tauType   = 'Composite';
                    [idP,idM] = section.beamColumnInteraction2d(data(i).axis,'AISC','CompPos');
                    BA_Elastic = BenchmarkAnalysis2d_Elastic_Sidesway_Inhibited(...
                        EIelastic,data(i).L,beta,0);
                    [P,M2] = BA_Elastic.determinePeakLoadWithEccentricity(idM,idP,e,Pno,tauType);

                    data(i).AISC2016_EIelastic           = EIelastic;
                    data(i).AISC2016_Pno                 = Pno;
                    data(i).AISC2016_P                   = -P;
                    data(i).AISC2016_M1                  = -P*e;
                    data(i).AISC2016_M2                  = M2;
                    data(i).AISC2016_test_to_predicted   = -data(i).Pexp/P;
                end
                
                % Plastic Stress Distribution
                if compute_PSD
                    psd = section.plasticStressDistributionObject();
                    num_points = 50;
                    switch lower(data(i).axis)
                        case 'strong'
                            [P,M,~] = psd.interactionSweep(0,num_points);
                        case 'weak'
                            [P,~,M] = psd.interactionSweep(pi/2,num_points);
                        otherwise
                            error('Bad axis: %s',data(i).axis);
                    end
                    
                    id = interactionDiagram2d(M,-P);
                    Pmax = 1.1*max(-P);
                    [~,Ppsd] = id.findIntersection(linspace(0,e*Pmax,100),linspace(0,Pmax,100));
                    
                    data(i).PSD_P = P;
                    data(i).PSD_M = M;
                    data(i).PSD_test_to_predicted = data(i).Pexp/Ppsd;
                end                
                
                % Nonlinear Analysis - PfD
                if compute_Analysis_PfD
                    fsDefOpts = struct;
                    fsDefOpts.nf1 = 30;
                    fsDefOpts.SteelMaterialType = 'AbdelRahman';
                    fsDefOpts.ConcreteMaterialType = 'ProposedForDesign';
                    fsDefOpts.includePackageDefinition = false;
                    
                    sectionDef = FiberSectionDefinition(section,data(i).axis,1,1,fsDefOpts);
                    
                    analysisOpts = struct;
                    analysisOpts.absoluteStrainLimit = 0.05;
                    analysisOpts.deleteFilesAfterAnalysis = false;

                    data2 = data(i);
                    data2.frame_type = 'Sidesway_Inhibited';
                    data2.section    = section;
                    data2.beta       = beta;
                    data2.delta0     = data(i).L/1000;
                    BA_OpenSees = BenchmarkAnalysis2d_OpenSees(data2,sectionDef,analysisOpts);

                    if data(i).L/data(i).depth <= 3.5
                        iResults = BA_OpenSees.runSectionAnalysis('LimitPoint_Proportional',max([e data(i).L/1000]),[],1);

                        if ~iResults.limitPoint.good
                            iResults = BA_OpenSees.runSectionAnalysis('LimitPoint_Proportional2',max([e data(i).L/1000]),[],1);
                        end
                        
                        if ~iResults.limitPoint.good
                            error('PfD Limit Point Not Determined for Specimen %i - %s - %s\n',i,sectionType,memberType);
                        end
                    else
                        iResults = BA_OpenSees.runAnalysis('LimitPoint_Proportional',e,[],1);

                        if ~iResults.limitPoint.good
                            iResults = BA_OpenSees.runAnalysis('LimitPoint_Proportional',e,[],2);
                        end

                        if ~iResults.limitPoint.good
                            iResults = BA_OpenSees.runAnalysis('LimitPoint_Proportional',e,[],3);
                        end

                        if ~iResults.limitPoint.good
                            error('PfD Limit Point Not Determined for Specimen %i - %s - %s\n',i,sectionType,memberType);
                        end
                    end                

                    data(i).Analysis_PfD_P  = -iResults.limitPoint.P1;
                    data(i).Analysis_PfD_M2 =  iResults.limitPoint.M2;
                    data(i).Analysis_PfD_test_to_predicted = -data(i).Pexp/iResults.limitPoint.P1;
                end
                
            case 'Beams'

            case 'Other'

            otherwise
                error('Unknown member type: %s',memberType);
        end
        
    catch exception 
        fprintf('Error in building data for specicmen %i (%s, %s)\n',i,sectionType,memberType);
        rethrow(exception);
    end
        
end

%% Save Data
save(sprintf('%s_%s.mat',sectionType,memberType),'data')
