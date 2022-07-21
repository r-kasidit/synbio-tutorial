%% Synbio-tutorial
% Purpose: Model optimization for synthetic biology design
%
% Written by Nachon Raethong, 19-JUL-2022
%% Initialize constraint-based reconstruction and analysis toolbox
initCobraToolbox;

%% IMPORT THE MODEL
% For .mat model file 
load 'model/iTN656.mat'
matModel = model
clear model;

% For .xml model file 
% Import model by RAVEN
xmlModel_1raven = importModel('model/iTN656.xml',false)

% Import model by COBRA
xmlModel_2cobra = readCbModel('model/iTN656.xml')
xmlModel_2raven = ravenCobraWrapper(xmlModel_2cobra) %Converting COBRA structure to RAVEN

% For .xlsx model file 
xlsxModel = importExcelModel('model/iTN656.xlsx')


%% Define constraints
%  The simulation was constrainted according to the experiment:
%  (i)   Anaerobic condition 
%  (ii)  Complex media composed of a complement of carbon source (e.g. glucose or sucrose),
%        amino acids, vitamins, lipids and ions.


new7maxsOUT = xmlModel_1raven;
new7maxsOUT = setParam(new7maxsOUT,'lb',{'EXC_BOTH_o2_e'},[-0.000000000001]); 

amino_acids = {'EXC_BOTH_arg__L_e', 'EXC_BOTH_cys__L_e', 'EXC_BOTH_ile__L_e', 'EXC_BOTH_leu__L_e', 'EXC_BOTH_lys__L_e', 'EXC_BOTH_met__L_e', 'EXC_BOTH_thr__L_e', 'EXC_BOTH_tyr__L_e', 'EXC_BOTH_val__L_e', 'EXC_BOTH_ala__L_e', 'EXC_BOTH_asn__L_e', 'EXC_BOTH_phe__L_e', 'EXC_BOTH_trp__L_e', 'EXC_BOTH_pro__L_e', 'EXC_BOTH_gln__L_e', 'EXC_BOTH_asp__L_e', 'EXC_BOTH_gly_e', 'EXC_BOTH_ser__L_e', 'EXC_BOTH_glu__L_e', 'EXC_BOTH_his__L_e'};
for i = 1:numel(amino_acids)
    new7maxsOUT=setParam(new7maxsOUT,'lb',amino_acids(i),[-5]);
end
lipid= {'EXC_BOTH_hdcea_e', 'EXC_BOTH_ocdcya_e', 'EXC_BOTH_ocdctr_e'};
for i = 1:numel(lipid)
    new7maxsOUT=setParam(new7maxsOUT,'lb',lipid(i),[-5]);
end

vitamin = {'EXC_BOTH_4abz_e', 'EXC_BOTH_pydam_e', 'EXC_BOTH_thm_e', 'EXC_BOTH_pydam_e', 'EXC_BOTH_pnto__R_e', 'EXC_BOTH_nac_e', 'EXC_BOTH_btn_e'};
for i = 1:numel(vitamin)
    new7maxsOUT=setParam(new7maxsOUT,'lb',vitamin(i),[-0.0001]);
end

carbon_sources = {'EXC_BOTH_melib_e', 'EXC_BOTH_raffin_e',...
    'EXC_BOTH_glc__D_e','EXC_BOTH_13ppd_e',...
    'EXC_BOTH_glcn__D_e', 'EXC_BOTH_glyc_e', 'EXC_BOTH_lcts_e',...
    'EXC_BOTH_sucr_e', 'EXC_BOTH_tre_e', ...
    'EXC_BOTH_rib__D_e',...
    'EXC_BOTH_gal_e', 'EXC_BOTH_malt_e', 'EXC_BOTH_arab__L_e',...
    'EXC_BOTH_tre_e', 'EXC_BOTH_sucr_e',...
    'EXC_BOTH_raffin_e', 'EXC_BOTH_melib_e',...
    'EXC_BOTH_malt_e', 'EXC_BOTH_lcts_e',...
    'EXC_BOTH_fuc__L_e', 'EXC_BOTH_arab__L_e', 'EXC_BOTH_glyc_e',...
    'EXC_BOTH_rib__D_e', 'EXC_BOTH_mnl_e',...
    'EXC_BOTH_glc__D_e', 'EXC_BOTH_glcn__D_e', 'EXC_BOTH_gal_e',...
    'EXC_BOTH_fru_e'};

for i = 1:numel(carbon_sources)
    new7maxsOUT=setParam(new7maxsOUT,'lb',carbon_sources(i),[0]);
end


products = {'EXC_BOTH_hxan_e' 'EXC_BOTH_xan_e',...
    'EXC_BOTH_gcald_e', 'EXC_BOTH_btd__RR_e',...
    'EXC_BOTH_ac_e', 'EXC_BOTH_lac__D_e',...
    'EXC_BOTH_orot_e','EXC_BOTH_etoh_e', ...
    'EXC_BOTH_mal__L_e', 'EXC_BOTH_lac__L_e'};
for i = 1:numel(products)
    new7maxsOUT=setParam(new7maxsOUT,'lb',products(i),[0]);
    new7maxsOUT=setParam(new7maxsOUT,'ub',products(i),[1000]);
end

%% Optimize objective function
%  In constrained condition, the uptake rates of oxygen, glucose, 
%  amino acids, vitamins, lipids and ions are limited according to 
%  the experiment. Biomass reaction is selected as an Objective function 
%  for Optimizing KUB-AC5 Growth.

new7maxsOUT2 = setParam(new7maxsOUT,'obj',{'AC5_GROWTH'},[1]); 
new7maxsOUTg20 = setParam(new7maxsOUT2,'lb',{'EXC_BOTH_glc__D_e'},[-20]);
sol = solveLP(new7maxsOUTg20)
printFluxes(new7maxsOUTg20, sol.x, true, 10^-3)
fprintf(['umax = ' num2str(sol.f*-1) ' per hour' '\n']);

%setRavenSolver('cobra')
%% Model validation
%  The growth simulation was carried out by constrainted the model into 4 conditions
%  accroding to the measured uptake rate of the each carbon substrate 
%  i.e. glucose, sucrose, maltose and lactose

new7maxsOUT3 = setParam(new7maxsOUT2,'lb',{'EXC_BOTH_glc__D_e'},[-6.634]); 
new7maxsOUT4 = setParam(new7maxsOUT2,'lb',{'EXC_BOTH_sucr_e'},[-5.033]); 
new7maxsOUT5 = setParam(new7maxsOUT2,'lb',{'EXC_BOTH_lcts_e'},[-0.940]); 
new7maxsOUT6 = setParam(new7maxsOUT2,'lb',{'EXC_BOTH_malt_e'},[-3.310]); 

sol3 = solveLP(new7maxsOUT3);
sol4 = solveLP(new7maxsOUT4);
sol5 = solveLP(new7maxsOUT5);
sol6 = solveLP(new7maxsOUT6);
fprintf(['umax glucose = ' num2str(sol3.f*-1) ' per hour' '\n']);
fprintf(['umax sucrose = ' num2str(sol4.f*-1) ' per hour' '\n']);
fprintf(['umax lactose = ' num2str(sol5.f*-1) ' per hour' '\n']);
fprintf(['umax maltose = ' num2str(sol6.f*-1) ' per hour' '\n']);


% in vivo
%EXC_BOTH_glc__D_e 6.634±0.684 > umax 0.151±0.004 
%EXC_BOTH_lcts_e 0.940±0.322 > umax 0.078±0.005 
%EXC_BOTH_sucr_e 5.033±0.310 > umax 0.247±0.003 
%EXC_BOTH_malt_e 3.310±0.764 > umax 0.199±0.009 

%in silico
%umax glucose = 0.15237 per hour
%umax sucrose = 0.26709 per hour
%umax lactose = 0.076945 per hour
%umax maltose = 0.18302 per hourc
	%in vivo	in silico	%ERROR
%Glucose	0.151	0.152	0.899
%Lactose	0.078	0.077	1.371
%Sucrose	0.247	0.274	8.133
%Maltose	0.199	0.183	8.731
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Identifying the essential/preferable nutrients
% singleOmission
amino_acids = {'EXC_BOTH_arg__L_e', 'EXC_BOTH_cys__L_e',...
    'EXC_BOTH_ile__L_e', 'EXC_BOTH_leu__L_e', 'EXC_BOTH_lys__L_e',...
    'EXC_BOTH_met__L_e', 'EXC_BOTH_thr__L_e', 'EXC_BOTH_tyr__L_e',...
    'EXC_BOTH_val__L_e', 'EXC_BOTH_ala__L_e', 'EXC_BOTH_asn__L_e',...
    'EXC_BOTH_phe__L_e', 'EXC_BOTH_trp__L_e', 'EXC_BOTH_pro__L_e',...
    'EXC_BOTH_gln__L_e', 'EXC_BOTH_asp__L_e', 'EXC_BOTH_gly_e',...
    'EXC_BOTH_ser__L_e', 'EXC_BOTH_glu__L_e', 'EXC_BOTH_his__L_e'};
vitaminB = {'EXC_BOTH_pydam_e', 'EXC_BOTH_thm_e',...
    'EXC_BOTH_pnto__R_e', 'EXC_BOTH_nac_e', 'EXC_BOTH_btn_e'};

nutrients = union(amino_acids,vitaminB);

for i = 1:numel(nutrients)
    new7maxsOUT33=setParam(new7maxsOUT2,'lb',nutrients(i),[-1]);
    new7maxsOUT33=setParam(new7maxsOUT33,'ub',nutrients(i),[1000]);
end

omiModel=new7maxsOUT33;
omiModelr=setParam(omiModel,'lb',{'EXC_BOTH_pro__L_e'},[0]);
solomiModel = solveLP(omiModelr,1);

cs=cell(numel(nutrients),1);
mVersatile_1 = cell(numel(nutrients),1);
mVersatile_05 = cell(numel(nutrients),1);
mVersatile_025 = cell(numel(nutrients),1);
mVersatile_0125 = cell(numel(nutrients),1);
mVersatile_00625 = cell(numel(nutrients),1);
mVersatile_00312 = cell(numel(nutrients),1);
mVersatile_00156 = cell(numel(nutrients),1);
mVersatile_00078 = cell(numel(nutrients),1);
mVersatile_00039 = cell(numel(nutrients),1);
mVersatile_00019 = cell(numel(nutrients),1);
mVersatile_00001 = cell(numel(nutrients),1);
mVersatile_0 = cell(numel(nutrients),1);

for i = 1:numel(nutrients)
    model = omiModel;
    cs{i}=nutrients{i};
    model1 = setParam(model,'lb',nutrients{i},[-1]);
    model05 = setParam(model,'lb',nutrients{i},[-0.5]);
    model025 = setParam(model,'lb',nutrients{i},[-0.25]);
    model0125 = setParam(model,'lb',nutrients{i},[-0.125]);
    model00625 = setParam(model,'lb',nutrients{i},[-0.0625]);
    model00312 = setParam(model,'lb',nutrients{i},[-0.0312]);
    model00156 = setParam(model,'lb',nutrients{i},[-0.0156]);
    model00078 = setParam(model,'lb',nutrients{i},[-0.0078]);
    model00039 = setParam(model,'lb',nutrients{i},[-0.0039]);
    model00019 = setParam(model,'lb',nutrients{i},[-0.0019]);
    model00001 = setParam(model,'lb',nutrients{i},[-0.0001]);
    model0 = setParam(model,'lb',nutrients{i},[0]);
    
    sol1 = solveLP(model1,1);
    sol05 = solveLP(model05,1);
    sol025 = solveLP(model025,1);
    sol0125 = solveLP(model0125,1);
    sol00625 = solveLP(model00625,1);
    sol00312 = solveLP(model00312,1);
    sol00156 = solveLP(model00156,1);
    sol00078 = solveLP(model00078,1);
    sol00039 = solveLP(model00039,1);
    sol00019 = solveLP(model00019,1);
    sol00001 = solveLP(model00001,1);
    sol0 = solveLP(model0,1);

    mVersatile_1{i} = (sol1.f*-1);
    mVersatile_05{i} = (sol05.f*-1);
    mVersatile_025{i} = (sol025.f*-1);
    mVersatile_0125{i} = (sol0125.f*-1);
    mVersatile_00625{i} = (sol00625.f*-1);
    mVersatile_00312{i} = (sol00312.f*-1);
    mVersatile_00156{i} = (sol00156.f*-1);
    mVersatile_00078{i} = (sol00078.f*-1);
    mVersatile_00039{i} = (sol00039.f*-1);
    mVersatile_00019{i} = (sol00019.f*-1);
    mVersatile_00001{i} = (sol00001.f*-1);
    mVersatile_0{i} = (sol0.f*-1);
end


singleOmissionResult = table(mVersatile_1,mVersatile_05,mVersatile_025,...
    mVersatile_0125,mVersatile_00625,mVersatile_00312,...
    mVersatile_00156,mVersatile_00078,mVersatile_00039,...
    mVersatile_00019,mVersatile_00001,mVersatile_0, 'RowNames', nutrients);
writetable(singleOmissionResult,'sigleOmissionResult.txt',...
    'Delimiter','tab',...
    'WriteRowNames',1);

% doubleOmission
mVersatile = cell(numel(nutrients),numel(nutrients));
row=cell(numel(nutrients),1);
column=cell(numel(nutrients),1);
for i = 1:numel(nutrients)
    model = omiModel;
    row{i} = {nutrients{i}};
    model_i = setParam(model,'lb',nutrients{i},[0]);
    for j = 1:numel(nutrients)
        column{j} = {nutrients{j}};
        model_ij = setParam(model_i,'lb',nutrients{j},[0]);
        solmodel_ij = solveLP(model_ij,1);
        mVersatile(i,j) = {(solmodel_ij.f*-1)};
    end
end

doubleOmissionResult = table(mVersatile, 'RowNames', nutrients);
writetable(doubleOmissionResult,'doubleOmissionResult.txt',...
    'Delimiter','tab',...
    'WriteRowNames',1);

%% get GSC from iTN656 for Integrative transcriptomics
% Purpose: to extract a gene set collection from iTN656
% model=importModel('model/iTN656.xml',false);
extractMetaboliteGSC(xmlModel_1raven,'','GSC.txt');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%