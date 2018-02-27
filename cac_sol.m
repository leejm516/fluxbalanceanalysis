load('iCac802_mod_180227.mat')

% Objective function
iCac802_mod.c(1307) = 0; % Biomass
iCac802_mod.c(1442) = 1; % Butyrate

% Bcd reaction
iCac802_mod = addReaction(iCac802_mod, 'R1017', '2 nadh_c + fdxox_c + butd2coa_c -> 2 nad_c + fdxrd_c + btcoa_c');

% Butyrate production
iCac802_mod.ub(1442)=1000; 
iCac802_mod.lb(1442)=0;

% Acetate production
iCac802_mod.ub(1450)=1000; 
iCac802_mod.lb(1450)=0;

% Solvent production
%iCac802_mod.ub(1459)=0; 
%iCac802_mod.ub(1457)=0; 

% Formate dehydrogenase
iCac802_mod.ub(1156)=0; 
iCac802_mod.lb(1156)=0;

iCac802_mod.ub(1157)=0; % Ferredoxin-NAD+ reductase
iCac802_mod.ub(1176)=0; % Ferredoxin-NADP+ reductase

% Disable the theronine aldolase rxn and add the pser transaminase rxn
iCac802_mod.lb(916)=0; 
iCac802_mod.ub(916)=0;
iCac802_mod = addReaction(iCac802_mod, 'PSERT', '3php_c + glu-L_c <=> akg_c + pser-L_c');

% phosphoglycerate kinase
iCac802_mod.lb(155)=-1000;
iCac802_mod.rev(155)=1;

% Redefine GAP dehydrogenase rxn
iCac802_mod = addReaction(iCac802_mod, 'R0742', 'h2o_c + nad_c + g3p_c -> 3pg_c + h_c + nadh_c');


% Nitrogenase
iCac802_mod.lb(49)=0;

% R1053
iCac802_mod.lb(767)=0;
iCac802_mod.ub(767)=0;


% Modify transhydrogenase rxn (assuming existence of nfnAB)
iCac802_mod = addReaction(iCac802_mod, 'R1787', 'fdxrd_c + 2 nadp_c + nadh_c + h_c -> fdxox_c + 2 nadph_c + nad_c');
% Disable transhydrogenase
%iCac802_mod.lb(1345)=0;
%iCac802_mod.ub(1345)=0;



% CoA Transferase (Acetate)
iCac802_mod.lb(1215)=0;

% Glutamate dehydrogenase
iCac802_mod.lb(158)=-1000;
iCac802_mod.lb(159)=0;
iCac802_mod.ub(159)=0;
iCac802_mod.rev(158)=1;


% Glutamate synthase
iCac802_mod.ub(388)=0;

% Disable the wrong rxn of folate reductase 
iCac802_mod.lb(820)=0;
iCac802_mod.ub(820)=0;
iCac802_mod.lb(821)=0;
iCac802_mod.ub(821)=0;

% Mannitol PTS
iCac802_mod.lb(36)=0;

% 5,10-MethyleneTHF reductase
iCac802_mod = addReaction(iCac802_mod, 'R0085', 'nadp_c + 5mthf_c <=> nadph_c + h_c + mlthf_c');

% Nucleoside kinase
iCac802_mod.lb(783)=0;
iCac802_mod.lb(1144)=0;



sol = optimizeCbModel(iCac802_mod, 'max', 0, 0)
