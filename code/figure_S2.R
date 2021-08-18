
# Read in data
clinical <- '~/Desktop/repos/Klebsiella_2021/data/clinical_maxfit_reps/flux_samples.tsv'
clinical <- read.delim(clinical, sep='\t', header=TRUE, row.names=1)
clinical <- as.data.frame(apply(clinical, 2, as.numeric))
laboratory <- '~/Desktop/repos/Klebsiella_2021/data/laboratory_maxfit_reps/flux_samples.tsv'
laboratory <- read.delim(laboratory, sep='\t', header=TRUE, row.names=1)
laboratory <- as.data.frame(apply(laboratory, 2, as.numeric))

# Subset to exchange reactions
exchanges <- c('EX_cd2_e','EX_cbi_e','EX_cbl1_e','EX_cgly_e','EX_cellb_e','EX_chol_e','EX_ca2_e','EX_12ppd__S_e',
               'EX_13ppd_e','EX_14glucan_e','EX_15dap_e','EX_23camp_e','EX_23ccmp_e','EX_23cgmp_e','EX_cit_e','EX_cl_e',
               'EX_23cump_e','EX_23dappa_e','EX_26dap__M_e','EX_2ddglcn_e','EX_cmp_e','EX_co2_e','EX_34dhpac_e',
               'EX_glyc2p_e','EX_3amp_e','EX_3cmp_e','EX_3gmp_e','EX_3hcinnm_e','EX_3hpppn_e','EX_cobalt2_e',
               'EX_colipa_e','EX_crn_e','EX_csn_e','EX_cu_e','EX_cu2_e','EX_cyan_e','EX_cynt_e','EX_3ump_e','EX_4abut_e',
               'EX_glyc3p_e','EX_4hoxpacd_e','EX_cys__D_e','EX_cys__L_e','EX_cytd_e','EX_glyclt_e','EX_dad_2_e',
               'EX_4hphac_e','EX_gmp_e','EX_5dglcn_e','EX_damp_e','EX_dca_e','EX_LalaDgluMdap_e','EX_dcmp_e','EX_gsn_e',
               'EX_dcyt_e','EX_ddca_e','EX_dgmp_e','EX_dgsn_e','EX_LalaDgluMdapDala_e','EX_abt__D_e','EX_ac_e',
               'EX_acald_e','EX_acgal_e','EX_acgal1p_e','EX_dha_e','EX_gthox_e','EX_diact_e','EX_dimp_e','EX_acgam_e',
               'EX_acgam1p_e','EX_din_e','EX_acmana_e','EX_gthrd_e','EX_acmum_e','EX_dms_e','EX_dmso_e','EX_acolipa_e',
               'EX_dopa_e','EX_gtp_e','EX_dtmp_e','EX_acser_e','EX_dump_e','EX_ade_e','EX_gua_e','EX_adn_e','EX_adocbl_e',
               'EX_duri_e','EX_etha_e','EX_ethso3_e','EX_ag_e','EX_etoh_e','EX_h_e','EX_f6p_e','EX_fald_e',
               'EX_agm_e','EX_akg_e','EX_ala_B_e','EX_ala__D_e','EX_ala__L_e','EX_alaala_e','EX_fe2_e','EX_fe3_e',
               'EX_fe3dhbzs_e','EX_fe3hox_e','EX_fecrm_e','EX_all__D_e','EX_amp_e','EX_feenter_e','EX_anhgm_e',
               'EX_feoxam_e','EX_arab__L_e','EX_h2_e','EX_arg__L_e','EX_for_e','EX_ascb__L_e','EX_fru_e','EX_h2o_e',
               'EX_frulys_e','EX_asn__L_e','EX_aso3_e','EX_h2o2_e','EX_asp__L_e','EX_fruur_e','EX_fuc__L_e','EX_btd_RR_e',
               'EX_h2s_e','EX_but_e','EX_fum_e','EX_g1p_e','EX_butso3_e','EX_hacolipa_e','EX_bz_e','EX_g3pc_e',
               'EX_g3pe_e','EX_halipa_e','EX_g3pg_e','EX_g3pi_e','EX_g3ps_e','EX_pi_e','EX_pnto__R_e','EX_g6p_e',
               'EX_ppa_e','EX_gal_e','EX_hdca_e','EX_gal_bD_e','EX_ppal_e','EX_pppn_e','EX_ppt_e','EX_gal1p_e',
               'EX_hdcea_e','EX_galct__D_e','EX_galctn__D_e','EX_galctn__L_e','EX_galt_e','EX_pro__L_e','EX_progly_e',
               'EX_psclys_e','EX_pser__L_e','EX_galur_e','EX_gam_e','EX_ptrc_e','EX_hg2_e','EX_pyr_e','EX_r5p_e',
               'EX_gam6p_e','EX_gdp_e','EX_rbt_e','EX_his__L_e','EX_rib__D_e','EX_rmn_e','EX_glc__D_e','EX_glcn_e',
               'EX_glcr_e','EX_sbt__D_e','EX_glcur_e','EX_hom__L_e','EX_glcur1p_e','EX_gln__L_e','EX_glu__L_e','EX_gly_e',
               'EX_ser__D_e','EX_ser__L_e','EX_skm_e','EX_so2_e','EX_so3_e','EX_so4_e','EX_glyald_e','EX_glyb_e',
               'EX_spmd_e','EX_hxa_e','EX_succ_e','EX_sucr_e','EX_sulfac_e','EX_glyc_e','EX_glyc__R_e','EX_hxan_e',
               'EX_tartr__L_e','EX_taur_e','EX_tcynt_e','EX_idon__L_e','EX_thm_e','EX_thr__L_e','EX_ile__L_e',
               'EX_thrp_e','EX_thym_e','EX_imp_e','EX_thymd_e','EX_indole_e','EX_tma_e','EX_tmao_e','EX_tre_e',
               'EX_inost_e','EX_trp__L_e','EX_tsul_e','EX_ins_e','EX_ttdca_e','EX_ttdcea_e','EX_isetac_e','EX_tungs_e',
               'EX_tym_e','EX_tyr__L_e','EX_tyrp_e','EX_uacgam_e','EX_udpacgal_e','EX_k_e','EX_udpg_e','EX_udpgal_e',
               'EX_kdo2lipid4_e','EX_udpglcur_e','EX_ump_e','EX_lac__D_e','EX_ura_e','EX_urea_e','EX_lac__L_e',
               'EX_uri_e','EX_val__L_e','EX_lcts_e','EX_xan_e','EX_xmp_e','EX_xtsn_e','EX_leu__L_e','EX_xyl__D_e',
               'EX_xylu__L_e','EX_lipa_e','EX_zn2_e','EX_lipa_cold_e','EX_lys__L_e','EX_lyx__L_e','EX_mal__D_e',
               'EX_mal__L_e','EX_malt_e','EX_malthx_e','EX_maltpt_e','EX_malttr_e','EX_maltttr_e','EX_man_e',
               'EX_man6p_e','EX_manglyc_e','EX_melib_e','EX_met__D_e','EX_met__L_e','EX_metsox_R__L_e','EX_metsox_S__L_e',
               'EX_mg2_e','EX_mmet_e','EX_mn2_e','EX_mnl_e','EX_mobd_e','EX_mso3_e','EX_n2o_e','EX_na1_e',
               'EX_nac_e','EX_nh4_e','EX_ni2_e','EX_nmn_e','EX_no_e','EX_no2_e','EX_no3_e','EX_o2_e','EX_o2s_e',
               'EX_ocdca_e','EX_ocdcea_e','EX_octa_e','EX_orn_e','EX_orot_e','EX_pacald_e','EX_peamn_e','EX_phe__L_e',
               'EX_pheme_e')
clinical_exchanges <- intersect(exchanges, colnames(clinical))
clinical <- clinical[, clinical_exchanges]
laboratory_exchanges <- intersect(exchanges, colnames(laboratory))
laboratory <- laboratory[, laboratory_exchanges]
rm(exchanges, clinical_exchanges, laboratory_exchanges)

# Separate into conserved and unique reactions
core_rxns <- intersect(colnames(clinical), colnames(laboratory))
clinical_core <- clinical[, core_rxns]
clinical_core$BIOMASS_ <- NULL
laboratory_core <- laboratory[, core_rxns]
laboratory_core$BIOMASS_ <- NULL
rm(core_rxns)

# Merge data for supervised learning
clinical_core$condition <- 1
laboratory_core$condition <- 0
flux_samples <- rbind(clinical_core, laboratory_core)
flux_samples$condition <- as.factor(flux_samples$condition)

# Run AUCRF and obtain feature lists
library(AUCRF)
set.seed(906801)
all_aucrf <- AUCRF(condition ~ ., data=flux_samples, pdel=0, k0=15)
flux_samples$condition <- NULL
print(all_aucrf)

# Assemble feature table
top_rxns_importance <- all_aucrf$ranking
rf_rxns <- as.data.frame(cbind(labels(top_rxns_importance), as.vector(top_rxns_importance)))
colnames(rf_rxns) <- c('id','mda')
rf_rxns$mda <- as.numeric(as.character(rf_rxns$mda))
rf_rxns <- subset(rf_rxns, mda > abs(min(rf_rxns$mda)))
#rf_rxns <- rf_rxns[order(-rf_rxns$mda),][1:15,]
rm(all_aucrf, top_rxns_importance)
#write.table(rf_rxns, file='/home/mjenior/Desktop/repos/Klebsiella_2021/data/AA_RF_results.tsv', 
#            quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE)

# Save time with existing RF results
rf_results <- read.delim('~/Desktop/repos/Klebsiella_2021/data/AA_RF_results.tsv', sep='\t', header=TRUE)
rf_results <- rf_results[order(rf_results$mda),] 


pdf(file='~/Desktop/repos/Klebsiella_2021/results/Figure_S4.pdf', width=2.5, height=4)
par(mar=c(3, 1, 1.5, 1), mgp=c(1.4, 0.5, 0), xpd=FALSE, lwd=1.7)
dotchart(rf_results$mda, bg='gray49', xlim=c(0,50), main='Amino Acid Exchanges', cex.main=1.1,
         pch=21, lwd=1.7, pt.cex=1.2, cex=0.8, xaxt='n')
axis(side=1, at=seq(0,50,10), cex.axis=0.8, lwd=1.7)
text(x=-1, y=seq(1.4,14.4,1), labels=rf_results$substrate, pos=4, cex=0.9)
mtext('Mean Decrease Accuracy (%)', side=1, padj=2.5, cex=0.7)
legend('bottomright', legend='OOB-AUCopt = 1.0', bty='n', pt.cex=0, cex=0.7)
dev.off()



clinical_leucine <- clinical_core$EX_leu__L_e
clinical_valine <- clinical_core$EX_val__L_e
clinical_arginine <- clinical_core$EX_arg__L_e
laboratory_leucine <- laboratory_core$EX_leu__L_e
laboratory_valine <- laboratory_core$EX_val__L_e
laboratory_arginine <- laboratory_core$EX_arg__L_e	



