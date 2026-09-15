clear; clc;

%Author: Martina Drecogna, U of Pavia and Padova
%Date July 31, 2026
%martina.drecogna@unipd.it

%% MATLAB computation

inputFileName = 'OGTT_input_sample10.csv';
opts = detectImportOptions(inputFileName);
opts.VariableNamingRule = 'preserve';
MATData = readtable(inputFileName, opts);

% Glucose (mg/dL)
G0   = MATData.glu0;
G30  = MATData.glu30;
G60  = MATData.glu60;
G90  = MATData.glu90;
G120 = MATData.glu120;

% Insulin (uU/mL)
I0   = MATData.ins0;
I30  = MATData.ins30;
I60  = MATData.ins60;
I90  = MATData.ins90;
I120 = MATData.ins120;

% C-peptide - two representations are kept, for different purposes:
%
%  1) *_ng: "raw" ng/mL, exactly as reported in the file. Used for the
%     ORIGINAL published Index60/DPTRS formula (Sosenko 2015; confirmed
%     by Redondo et al., Diabetologia 2021), which is natively defined
%     in ng/mL - no conversion needed here.
%  2) CP0/CP30/CP60/CP120: pmol/L, converted by multiplying by 331.
%     Used to compare against the AOMM columns in pmol/L
%     (C_peptide_Index, AUC_CP_incremental).

CP0_ng   = MATData.cpep0;
CP30_ng  = MATData.cpep30;
CP60_ng  = MATData.cpep60;
CP120_ng = MATData.cpep120;

CP0   = CP0_ng   .* 331;
CP30  = CP30_ng  .* 331;
CP60  = CP60_ng  .* 331;
CP120 = CP120_ng .* 331;

% Anthropometrics and demographics
age = MATData.age_at_visit;
BMI = MATData.bmi;

%% Static indices
MATData.HOMA_IR = (G0 .* I0) / 405;
MATData.HOMA_B  = (360 .* I0) ./ (G0 - 63);
MATData.QUICKI  = 1 ./ (log10(I0) + log10(G0));

%% Empirical dynamic indices
% Mean glucose and insulin during the OGTT
G_mean = (G0 + G30 + G60 + G90 + G120) / 5;
I_mean = (I0 + I30 + I60 + I90 + I120) / 5;

% Matsuda ISI
MATData.Matsuda_ISI = 10000 ./ sqrt((G0 .* I0) .* (G_mean .* I_mean));

% IGI (Insulinogenic Index)
MATData.IGI = (I30 - I0) ./ (G30 - G0);

% Matsuda Disposition Index
MATData.Matsuda_Disposition = MATData.IGI .* MATData.Matsuda_ISI;

% C-peptide Index
MATData.C_peptide_Index = (CP30 - CP0) ./ (G30 - G0); % pmol/L per mg/dL

% IGI / HOMA-IR
MATData.IGI_HOMA_IR = MATData.IGI ./ MATData.HOMA_IR;

%% Incremental areas under the curve (iAUC) and ISSI-2
t = [0, 30, 60, 90, 120];

MATData.AUC_G = 15*(G0 + G30) + 15*(G30 + G60) + 15*(G60 + G90) + 15*(G90 + G120);
MATData.AUC_G_incremental = MATData.AUC_G - (G0 * 120);

MATData.AUC_I = 15*(I0 + I30) + 15*(I30 + I60) + 15*(I60 + I90) + 15*(I90 + I120);
MATData.AUC_I_incremental = MATData.AUC_I - (I0 * 120);

MATData.AUC_CP = 15*(CP0 + CP30) + 15*(CP30 + CP60) + 30*(CP60 + CP120);
MATData.AUC_CP_incremental = MATData.AUC_CP - (CP0 * 120); % pmol/L*min

% ISSI-2
MATData.ISSI_2 = (MATData.AUC_I ./ MATData.AUC_G) .* MATData.Matsuda_ISI;

%% Other glycemic indices
MATData.Glucose_1hr = G60;
MATData.Glucose_2hr = G120;

% Peak Glucose and Time to Peak (row by row)
G_matrix = [G0, G30, G60, G90, G120];
[MATData.Peak_Glucose, max_idx] = max(G_matrix, [], 2);
MATData.Time_to_Peak_Glucose = t(max_idx)';

MATData.Glucose_Range = MATData.Peak_Glucose - min(G_matrix, [], 2);

%% T1D progression indices
% The ORIGINAL formula is used here, as published (Sosenko et al.,
% Diabetes Care 2015; explicitly confirmed - logₑ, i.e. natural log,
% C-peptide in ng/mL - by Redondo et al., Diabetologia 2021),
% NOT the "recast" version with coefficients rescaled for pmol/L used
% internally by AOMM (which is only an equivalent algebraic
% reparameterization, see src/AOMM_gui/DatasetProcessor.cs lines ~452-503).
% CP0_ng/CP60_ng are therefore used (raw ng/mL, no conversion needed:
% the input file is already in ng/mL) with the paper's literal
% coefficients, without any unit conversion constant.

MATData.Index60 = 0.3695 .* log(CP0_ng) + 0.0165 .* G60 - 0.3644 .* CP60_ng;

MATData.DPTRS = 1.571 .* log(CP0_ng + 1) + 0.115 .* G120 - 0.193 .* (CP60_ng ./ CP0_ng) ...
             + 0.0019 .* age + 1.5 .* log(BMI);


%% AOMM results

aommFileName = 'AOMM_results_summary_sample10.csv';
opts = detectImportOptions(aommFileName);
opts.VariableNamingRule = 'preserve';
AOMMData = readtable(aommFileName, opts);

mapping = { ...
    'HOMA_IR',                          'HOMA_IR'; ...
    'HOMA_beta_pct',                    'HOMA_B'; ...
    'QUICKI',                           'QUICKI'; ...

    'Matsuda_ISI',                      'Matsuda_ISI'; ...
    'Insulinogenic_index_mcU_per_mg',   'IGI'; ...
    'Matsuda_Disposition_Index',        'Matsuda_Disposition'; ...
    'C_peptide_index_pmol_per_mg',      'C_peptide_Index'; ...
    'ISSI_2_Retnakaran2009',            'ISSI_2'; ...
    'IGI_over_HOMA_IR_Wareham',         'IGI_HOMA_IR'; ...

    'Glucose_iAUC_0_120_mg_dL_min',     'AUC_G_incremental'; ...
    'Insulin_iAUC_0_120_mcU_mL_min',    'AUC_I_incremental'; ...
    'C_peptide_iAUC_0_120_pmol_L_min',  'AUC_CP_incremental'; ...

    'Glucose_60min_mg_dL',              'Glucose_1hr'; ...
    'Glucose_peak_mg_dL',               'Peak_Glucose'; ...
    'Glucose_time_to_peak_min',         'Time_to_Peak_Glucose'; ...
    'Glucose_range_mg_dL',              'Glucose_Range'; ...
    'Glucose_120min_mg_dL',             'Glucose_2hr'; ...

    'Index60_Sosenko2015',              'Index60'; ...
    'DPTRS_Sosenko2008',                'DPTRS' ...
};

% rename columns
for k = 1:size(mapping, 1)
    oldName = mapping{k, 1};
    newName = mapping{k, 2};
    if ismember(oldName, AOMMData.Properties.VariableNames)
        AOMMData.Properties.VariableNames{oldName} = newName;
    end
end

% align by ID
[in_both, idx_csv] = ismember(MATData.id, AOMMData.Subject_ID);
MATData = MATData(in_both, :);
AOMMData = AOMMData(idx_csv(in_both), :);

variableNames = mapping(:, 2);
figure('Units', 'normalized', 'Position', [0.05, 0.05, 0.9, 0.85], 'Name', 'Index of correlation');
tlo = tiledlayout(5, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

%% correlation

for i = 1:length(variableNames)
    varName = variableNames{i};

    % if the variable exists in both datasets
    if ismember(varName, AOMMData.Properties.VariableNames) && ismember(varName, MATData.Properties.VariableNames)

        aommValue = AOMMData.(varName);
        matlabValue = MATData.(varName);

        validIdx = ~isnan(aommValue) & ~isnan(matlabValue);
        x = aommValue(validIdx);
        y = matlabValue(validIdx);

        nexttile;
        % Scatter plot
        scatter(x, y, 35, 'filled', 'MarkerFaceColor', [0 0.4470 0.7410], 'MarkerFaceAlpha', 0.6);
        hold on;
        minVal = min([x; y]);
        maxVal = max([x; y]);
        plot([minVal, maxVal], [minVal, maxVal], 'r--', 'LineWidth', 1.5); % y = x

        % Pearson correlation coefficient R
        if length(x) > 1
            R_matrix = corrcoef(x, y);
            R = R_matrix(1, 2);
            title(sprintf('%s (R = %.4f)', strrep(varName, '_', '\_'), R), 'FontSize', 9);
        else
            title(strrep(varName, '_', '\_'), 'FontSize', 9);
        end

        grid on;
        xlabel('AOMM', 'FontSize', 8);
        ylabel('MATLAB', 'FontSize', 8);
    end
end
