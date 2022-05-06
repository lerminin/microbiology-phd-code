#############################################
## 
## This script describe alpha and beta diversity
## analysis with 16S community profiling data
##
#############################################

# Load all required packages:
pkgs <- c("ggplot2", "ggprism", "here", "dplyr", "ggpubr", "readxl", "purrr", "filesstrings", "readr", "tibble", "broom", "lme4", "optimx", "emmeans", "expss", "phyloseq", "vegan", "SRS", "tidyr", "DHARMa")
lapply(pkgs, library, character = TRUE)

sessionInfo()
# R version 4.1.2 (2021-11-01)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19044)
# 
# Matrix products: default
# 
# locale:
#   [1] LC_COLLATE=English_Canada.1252  LC_CTYPE=English_Canada.1252   
# [3] LC_MONETARY=English_Canada.1252 LC_NUMERIC=C                   
# [5] LC_TIME=English_Canada.1252    
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] DHARMa_0.4.5          tidyr_1.1.4           SRS_0.2.2             shinybusy_0.3.0      
# [5] shinycssloaders_1.0.0 DT_0.20               shiny_1.7.1           vegan_2.5-7          
# [9] lattice_0.20-45       permute_0.9-7         phyloseq_1.38.0       expss_0.11.1         
# [13] maditr_0.8.2          emmeans_1.7.2         optimx_2021-10.12     lme4_1.1-27.1        
# [17] Matrix_1.3-4          broom_0.7.10          tibble_3.1.6          readr_2.1.0          
# [21] filesstrings_3.2.2    stringr_1.4.0         purrr_0.3.4           readxl_1.3.1         
# [25] ggpubr_0.4.0          dplyr_1.0.7           here_1.0.1            ggprism_1.0.3        
# [29] ggplot2_3.3.5        
# 
# loaded via a namespace (and not attached):
#   [1] TH.data_1.1-0          minqa_1.2.4            colorspace_2.0-2       ggsignif_0.6.3        
# [5] ellipsis_0.3.2         rprojroot_2.0.2        estimability_1.3       htmlTable_2.4.0       
# [9] XVector_0.34.0         rstudioapi_0.13        fansi_0.5.0            mvtnorm_1.1-3         
# [13] codetools_0.2-18       splines_4.1.2          knitr_1.37             ade4_1.7-18           
# [17] jsonlite_1.7.3         nloptr_1.2.2.3         cluster_2.1.2          compiler_4.1.2        
# [21] backports_1.3.0        assertthat_0.2.1       fastmap_1.1.0          later_1.3.0           
# [25] htmltools_0.5.2        tools_4.1.2            igraph_1.2.11          gtable_0.3.0          
# [29] glue_1.5.0             GenomeInfoDbData_1.2.7 reshape2_1.4.4         Rcpp_1.0.7            
# [33] carData_3.0-4          Biobase_2.54.0         cellranger_1.1.0       vctrs_0.3.8           
# [37] Biostrings_2.62.0      rhdf5filters_1.6.0     multtest_2.50.0        ape_5.5               
# [41] nlme_3.1-153           iterators_1.0.13       xfun_0.29              mime_0.12             
# [45] lifecycle_1.0.1        rstatix_0.7.0          MASS_7.3-54            zlibbioc_1.40.0       
# [49] zoo_1.8-9              scales_1.1.1           promises_1.2.0.1       hms_1.1.1             
# [53] parallel_4.1.2         biomformat_1.22.0      sandwich_3.0-1         rhdf5_2.38.0          
# [57] strex_1.4.2            stringi_1.7.6          S4Vectors_0.32.3       foreach_1.5.1         
# [61] checkmate_2.0.0        BiocGenerics_0.40.0    boot_1.3-28            GenomeInfoDb_1.30.0   
# [65] rlang_0.4.11           pkgconfig_2.0.3        matrixStats_0.61.0     bitops_1.0-7          
# [69] Rhdf5lib_1.16.0        htmlwidgets_1.5.4      tidyselect_1.1.1       plyr_1.8.6            
# [73] magrittr_2.0.1         R6_2.5.1               IRanges_2.28.0         generics_0.1.2        
# [77] multcomp_1.4-18        DBI_1.1.2              pillar_1.7.0           withr_2.4.3           
# [81] mgcv_1.8-38            survival_3.2-13        abind_1.4-5            RCurl_1.98-1.5        
# [85] crayon_1.5.0           car_3.0-12             utf8_1.2.2             tzdb_0.2.0            
# [89] grid_4.1.2             data.table_1.14.2      digest_0.6.27          xtable_1.8-4          
# [93] httpuv_1.6.3           numDeriv_2016.8-1.1    stats4_4.1.2           munsell_0.5.0         
# [97] languageR_1.5.0   

pal <- c("#808080",
         "#a6cee3",
          "#1f78b4",
          "#fb9a99",
          "#e31a1c",
          "#b2df8a",
          "#33a02c",
          "#fdbf6f",
          "#ff7f00",
          "#cab2d6",
          "#6a3d9a",
          "#ffff99",
          "#b15928",
          "#adf2f2",
          "#009698",
         "#DADADA")
# Load colourblind-friendly palette:
cbPalette <- c("#E69F00", "#56B4E9", "#D55E00", "#0072B2", "#CC79A7", "#009E73", "#F0E442", "#999999")

# Bring in the files (using dada2 assignTaxonomy files):
# OTU file
otu_mat <- read.csv(here("16S", "data", "Counts_no_contam.csv"))
# Taxonomy file
tax_mat <- read.csv(here("16S", "data", "Taxa_no_contam.csv"))
# Metadata file 
samples_df <- read_csv(here("Fieldwork", "output_data", "master.csv"))
head(samples_df)

# Remove non-Bacteria domain ASVs (eukaryotic and archaea)
nonbac <- as_vector(tax_mat[,2] == "Bacteria") # remove non-Bacteria
tax_mat <- subset(tax_mat, nonbac)
otu_mat <- subset(otu_mat, nonbac)
nonbac2 <- as_vector(!is.na(tax_mat[,3])) # remove phylum with NA
tax_mat <- subset(tax_mat, nonbac2)
otu_mat <- subset(otu_mat, nonbac2)
nonbac3 <- as_vector(tax_mat[,5] != "Chloroplast" | is.na(tax_mat[,5])) # remove Cyanobac chloroplasts
tax_mat <- subset(tax_mat, nonbac3)
otu_mat <- subset(otu_mat, nonbac3)

# Make the phyloseq objects by renaming rows in the dataframe and transforming into matrix:
row.names(otu_mat) <- otu_mat$X
otu_mat <- select(otu_mat, -c(X, Mock1, Mock2, Mock3, Mock4, Mock5, Neg2, Neg3, 
                              HG04REP2, HG04REP3, ML06REP2, ML06REP3, 
                              HG03BSAW, ST01, ST06,  RB02, RB03, # samples with low coverage from 2020 run
                              HG03,
                              Mock1, KitNeg1, KitNeg2, KitNeg3, SeqNeg1, SeqNeg2, SeqNeg3, Undetermined))

otu_mat <- otu_mat %>% dplyr::rename(
  #RB02s = TH1,
  #ST04s = TH2,
  #HG03s = TH3,
  #ML05s = TH4,
  #FH04s = TH5,
  ST01 = ST01BSAW,
  ST06 = ST06BSAW, 
  RB02 = RB02BSAW, 
  RB03 = RB03BSAW, 
  HG04 = HG04REP1,
  ML06 = ML06REP1) %>% 
  select(sort(names(.)))
colnames(otu_mat)

# Rarify using SRS: to 6146 reads (to 5090 reads if non-Bacteria and chloroplast ASVs are removed at start)
sort(rowSums(t(otu_mat))) # to check number of counts per sample
specnumber(t(otu_mat))
raremax <- min(rowSums(t(otu_mat)))
otu_srs <- SRS(otu_mat, raremax)
#rarecurve(t(otu_srs), step = 1, sample = raremax)
rownames(otu_srs) <- rownames(otu_mat)
otu_mat <- as.matrix(otu_srs)
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE)

row.names(tax_mat) <- tax_mat$X
tax_mat <- select(tax_mat, -c(X))
tax_mat <- as.matrix(tax_mat)
TAX <- tax_table(tax_mat)

# Remove 2021 samples and fix cow column 
samples_df <- samples_df %>% 
  mutate(source = "water") %>% 
  mutate(year = as.character(year))
# Add sediment samples
samples_sed <- samples_df %>% 
  mutate(sediment_collected = ifelse(is.na(sediment_collected), "yes", sediment_collected)) %>% 
  filter(sediment_collected == "yes") %>% 
  mutate(sampleID = paste(sampleID, "S", sep = "")) %>% 
  mutate(source = "sediment")
# Remove water samples with low number of reads
samples_wat <- samples_df %>% 
#   filter(sampleID != "ST01") %>% 
#   filter(sampleID != "ST06") %>% 
   filter(sampleID != "HG03") %>% 
#   filter(sampleID != "RB02") %>% 
#   filter(sampleID != "RB03") %>% 
    filter(!is.na(sampleID))
# Join water and sediment dfs together
samples_all <- full_join(samples_wat, samples_sed)

samples_df2 <- column_to_rownames(samples_all, var = "sampleID")
samples_df_ordered <- samples_df2[order(row.names(samples_df2)),]
samples <- sample_data(samples_df_ordered)

ps <- phyloseq(OTU, TAX, samples)
ps

# Prune ASVs not present in at least one sample:
ps <- prune_taxa(taxa_sums(ps) > 0, ps)
ps


################################
## Alpha diversity with vegan ##
################################

# Compute alpha diversity metrics
sample_data(ps)$invsimpson <- diversity(otu_srs, MARGIN=2, index = "invsimpson")
sample_data(ps)$simpson <- diversity(otu_srs, MARGIN=2, index = "simpson")
sample_data(ps)$shannon <- diversity(otu_srs, MARGIN=2, index = "shannon")

hist(sample_data(ps)$invsimpson) # gamma distribution

# # Modeling distribution to find which distribution is best (using fitdistrplus package):
# descdist(sample_data(ps)$invsimpson)
# descdist(log(sample_data(ps)$invsimpson))
# plot(fitdist(sample_data(ps)$invsimpson, distr = "gamma"))
# plot(fitdist(log(sample_data(ps)$invsimpson), distr = "gamma"))

# Making glms
df <- as_tibble(sample_data(ps), rownames = "sampleID") %>% 
  filter(!(siteID == "BV")) %>% 
  filter(!(sampleID == "ML04S")) # sample was taken from edge of dugout and is extremely high (>300) compared to the rest
m1 <- glm(invsimpson ~ source*siteID*year, data = df, family = Gamma(link = "identity"))
summary(m1)
pchisq(m1$deviance, df = m1$df.residual, lower.tail = F)# if chi is <0.05, model is bad fit
m2 <- glm(invsimpson ~ source*siteID*year, data = df, family = Gamma(link = "log"))
anova(m1, m2, test = "LRT")# believe this is saying m1 is better than m2

# Checking gamma model assumptions with DHARMa package
qqnorm(residuals(m1))
qqline(residuals(m1))
plot(simulateResiduals(fittedModel = m1, n = 500))
testDispersion(simulateResiduals(fittedModel = m1))
plot(predict(m1), residuals(m1, type = "deviance"))
pchisq(m1$deviance, df = m1$df.residual, lower.tail = F) # if chi is <0.05, model is bad fit
plotResiduals(simulateResiduals(fittedModel = m1, n = 500), form = factor(df$source))
plotResiduals(simulateResiduals(fittedModel = m1, n = 500), form = factor(df$siteID))
plotResiduals(simulateResiduals(fittedModel = m1, n = 500), form = factor(df$year))

# Calculate marginal means with package emmeans - for fixed effects 
ref_grid(m1) @ grid
emm <- emmeans(m1, ~ year*source*siteID, type = "response");emm
set.seed(25)
emmc <- contrast(emm, method = "pairwise", adjust = "tukey", by = c("year", "siteID"));emmc
write.csv(as_tibble(emmc), here("stats", "invsimpson_emmeans.csv"))
emme <- as_tibble(eff_size(emm, sigma = sigma(m1), edf = df.residual(m1), method = "pairwise", by = c("siteID", "year")));emme
emmip(emm, siteID ~ year + source, CIs = TRUE)
# emme$p_value <- "<0.0001"
# emme$p_value[3] <- "0.0001"
# emme$p_value[7] <- "0.0004"
# write.csv(emme, here("stats", "invsimpson_emmeans.csv"))

# Effect size plot
ggplot(emme, aes(siteID, effect.size)) +
  geom_point() +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=0, position = position_dodge(0.5)) +
  geom_hline(yintercept = 0) 

# Plot marginal means with year as separate panel:
ggplot(as_tibble(emm), aes(siteID, emmean, colour = source)) +
  geom_point(position = position_dodge(0.2)) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=0, position = position_dodge(0.2)) +
  geom_point(data = df, aes(x = siteID, y = invsimpson, colour = source), alpha = 0.2, position = position_jitter(h=0.1, w=0.1)) +
  scale_colour_manual(values = cbPalette) +
  labs(y = "Inverse Simpson diversity index", x = "Site") +
  theme_prism(base_line_size = 0.2, base_fontface = "plain") +
  facet_wrap(~ year)

# Plot marginal means with both years on same panel:
df <- df %>% 
  mutate(source_year = paste(source, year, sep = "_"))

# Load colourblind-friendly palette:
cbPalette <- c("#E69F00", "#D55E00", "#56B4E9", "#0072B2", "#CC79A7", "#009E73", "#F0E442", "#999999")

as_tibble(emm) %>% 
  mutate(source_year = paste(source, year, sep = "_")) %>% 
  ggplot(aes(siteID, emmean, colour = source_year, shape = source_year)) +
  geom_point(position = position_dodge(0.4)) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=0, position = position_dodge(0.4)) +
  geom_point(data = df, aes(x = siteID, y = invsimpson, colour = source_year, shape = source_year), alpha = 0.3, position = position_jitter(h=0.1, w=0.1)) +
  #scale_colour_manual(values = cbPalette) +
  labs(y = "Inverse Simpson diversity index", x = "Site") +
  theme_prism(base_line_size = 0.2, base_fontface = "plain") +
  scale_colour_manual(name = "Source & Year",
                      labels = c("Sediment 2020", "Sediment 2021", "Water 2020", "Water 2021"), values = cbPalette) +   
  scale_shape_manual(name = "Source & Year",
                     labels = c("Sediment 2020", "Sediment 2021", "Water 2020", "Water 2021"), values = c(19, 17, 19, 17)) +
  theme(legend.position = "top")

ggsave(here("figures_thesis", "16S_invsimpson.png"), height = 5, width = 9)


######################
## Ordination - RDA ##
######################

# Subset phyloseq object into water and sediment:
ps # main object
ps_wat <- subset_samples(ps, source == "water" & siteID != "BV")
ps_sed <- subset_samples(ps, source == "sediment")

# Transfer phyloseq object into vegan-compatible format 
df <- as_tibble(sample_data(ps_wat), rownames = "sampleID")
spp <- as.matrix(t(otu_table(ps_wat)))

# Hellinger transformation on spp counts (gives low weights to rare species):
#spe.hel <- decostand(chooseTaxa(varespec, n.occ = 5), method = "hellinger") # keeps spp present in at least 5 samples
spph <- decostand(spp, method = "hellinger")

# Make fake depth data for the four samples with NA values (forgot to bring piece of equipment to the field): (based on adjacent sampling and mean for 2021)
if(nrow(df) == 74){
  df$depth_m[40] <- 0.9
  df$depth_m[44] <- 1.0
  df$depth_m[54] <- 2.0
  df$depth_m[55] <- 2.3
}

# Calculate distance matrix:
x <- vegdist(spph, method="bray", binary=FALSE)
# look at an unconstrained ordination first, it is always a good idea to look at both unconstrained and constrained ordinations
# set the seed to reproduce the same result in the fture
set.seed(1234)
mds <- metaMDS(spph, distance = "bray", k = 3, trymax = 50)
# extract x and y coordinates from MDS plot into new dataframe, so you can plot with ggplot 
z <- data.frame(mds$points)
mds$stress # 0.137
# colour by sample type
ggplot(z, aes(MDS1, MDS2, color=df$siteID)) + geom_point() + theme_bw() + ggtitle('Unconstrained ordination, stress:0.15')

# RDA - removed chloride (colinear with TOC), avgchltotal (colinear with chla), DO_percent
simpleRDA <- rda(spph ~ scale(julianday) + scale(pH) + scale(microcystin_ppb) + scale(temp_degC) + scale(SO4) + scale(TIC_ppm) + scale(TOC_ppm) + scale(nitratenitrite_mgperL) + scale(alkalinity) + scale(avgchla_ugperL) + scale(totalnitrogen_ugperL) + scale(totalphosphorous_mgperL) + scale(ammonium_mgperL) + scale(solublereactivephosphate_mgperL) + scale(depth_m) + scale(cond_uS) + scale(DO_mgperL) + Condition(siteID), data=df)
summary(simpleRDA)
screeplot(simpleRDA) #bstick not available for constrained ordinations
anova.cca(simpleRDA)
anova.cca(simpleRDA, by = "axis") # first four are significant
anova.cca(simpleRDA, by = "margin")

# unadjusted & adjusted R^2 retrieved from the rda result
RsquareAdj(simpleRDA)$r.squared
RsquareAdj(simpleRDA)$adj.r.squared

# Proportion of variance explained by first two axes 
sum((as.vector(m4$CCA$eig)/sum(m4$CCA$eig))[1:2])

# Testing RDA for significance:
# for 2020 samples only:
h <- how(nperm=499, within = Within(type = "series", mirror = FALSE), blocks = df$siteID)
# for all samples:
h <- how(within = Within(type = "series", mirror = FALSE), plots = with(df, Plots(strata = year, type = "none")), blocks = df$siteID)


# Testing variables in RDA for collinearity - checking variance inflation factors (linear dependencies between enviro variables) - VIF>20 = strong collinearity, VIF>10 = potential concern
df_num <- select_if(df, is.numeric)  
df_num <- df_num[-c(2, 4, 10:12, 14)]
pairs(df_num[1:19])
vif.cca(simpleRDA)
sqrt(vif.cca(simpleRDA)) # if >2, multicollinearity is high
head(goodness(simpleRDA))
head(inertcomp(ord, proportional = TRUE))
languageR::pairscor.fnc(df) # only for numeric variables

# Variable selection in RDA:
m1 <- ordiR2step(rda(spph ~ 1, data = df), m2)
m1 <- ordistep(simpleRDA, rda(spph ~ 1, data = df), direction = "backward")
# best backwards model: spph ~ scale(julianday) + scale(TIC_ppm) + scale(nitratenitrite_mgperL) +      scale(totalnitrogen_ugperL) + scale(DO_mgperL) + Condition(siteID) 
# best model: spph ~ siteID + scale(TOC_ppm) + scale(julianday) + scale(microcystin_ppb) + scale(pH) + scale(DO_mgperL) + scale(ammonium_mgperL) + scale(chloride), data = df)
# best model without site & year: spph ~ scale(microcystin_ppb) + scale(avgchla_ugperL) + scale(chloride) +      scale(pH) + scale(julianday) + scale(TOC_ppm)
# best model with site as condition: scale(microcystin_ppb) + scale(avgchla_ugperL) + scale(TOC_ppm) + scale(pH) + scale(julianday) + scale(alkalinity)


m2 <- rda(spph ~ scale(microcystin_ppb) + scale(avgchla_ugperL) + scale(pH) + scale(julianday) + scale(TOC_ppm) + scale(alkalinity) + Condition(siteID), data = df)
m3 <- rda(spph ~ scale(depth_m) + scale(temp_degC) + scale(DO_mgperL) + scale(pH) + scale(julianday) + scale(avgchla_ugperL) + scale(microcystin_ppb) + scale(TIC_ppm) + scale(SO4) + scale(ammonium_mgperL) + scale(solublereactivephosphate_mgperL) + scale(nitratenitrite_mgperL), data = df)
# All water samples:
m4 <- rda(spph ~ scale(microcystin_ppb) + depth_m + scale(avgchla_ugperL) + scale(pH) + scale(julianday) + scale(cond_uS) + scale(DO_mgperL) + scale(temp_degC) + Condition(siteID + year), data = df)
m1 <- ordiR2step(rda(spph ~ 1, data = df), m4, direction = "forward")
# best model: 


# Plotting triplot:
library(ggord)
# Triplot: three different entities in the plot: sites, response variables and explanatory variables (arrowheads are on the explanatory variables)
# Scaling 1
plot(simpleRDA, scaling=1, main="Triplot RDA matrix ~ env - scaling 1 - wa scores")
spe.sc <- scores(simpleRDA, choices=1:2, scaling=1, display="sp")
arrows(0,0,spe.sc[,1], spe.sc[,2], length=0, lty=1, col='red')
# Scaling 2
plot(simpleRDA, main="Triplot RDA matrix ~ env - scaling 2 - wa scores")
spe2.sc <- scores(simpleRDA, choices=1:2, display="sp") # scores() choices= indicates which axes are to be selected, make sure to specify the scaling if its different than 2 
arrows(0,0,spe2.sc[,1], spe2.sc[,2], length=0, lty=1, col='red')

# Alt plot:
plot(m4, display = "sites", type = "text")
plot(m4, display = "species", type = "text")
# extract scores:
scrs <- scores(vegan_rda, display = "species", choices = c(1,2), scaling = "species", correlation = TRUE)
head (scrs)

library(ggord)
ggord(m1, df$siteID, vec_ext=0.7, ellipse=TRUE, size=1, addsize=-1) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 


############################
## Astrobio Mike tutorial ##
##  Taxonomic summaries   ##
############################

# To make a bar chart by proportion of phylum and break down Proteobacteria by class:

# Tax: use phyloseq to make a count table that has summed all ASVs that were in the same phylum
phyla_counts_tab <- otu_table(tax_glom(ps, taxrank="phylum"))

# Tax: making a vector of phyla names to set as row names
phyla_tax_vec <- as.vector(tax_table(tax_glom(ps, taxrank="phylum"))[,2]) 
rownames(phyla_counts_tab) <- as.vector(phyla_tax_vec)

# Tax: make a new category for NA phylum
unclassified_tax_counts <- colSums(otu_table(ps)) - colSums(phyla_counts_tab)
phyla_and_unidentified_counts_tab <- rbind(phyla_counts_tab, "Unclassified phyla"=unclassified_tax_counts)

# Tax: remove the Proteobacteria, so they can be added back broken down by class
temp_major_taxa_counts_tab <- phyla_and_unidentified_counts_tab[!row.names(phyla_and_unidentified_counts_tab) %in% "Proteobacteria", ]

# Tax: make a count table broken down by class (contains classes beyond the Proteobacteria too at this point)
class_counts_tab <- otu_table(tax_glom(ps, taxrank="class"))

# Tax: make a table that holds the phylum & class info
class_tax_phy_tab <- tax_table(tax_glom(ps, taxrank="class"))
phy_tmp_vec <- class_tax_phy_tab[,2]
class_tmp_vec <- class_tax_phy_tab[,3]
rows_tmp <- row.names(class_tax_phy_tab)
class_tax_tab <- data.frame("phylum"=phy_tmp_vec, "class"=class_tmp_vec, row.names = rows_tmp)

# Tax: make a vector of just the Proteobacteria classes
proteo_classes_vec <- as.vector(class_tax_tab[class_tax_tab$phylum == "Proteobacteria", "class"])

# Tax: change the row names like above so that they correspond to the taxonomy, rather than an ASV identifier
rownames(class_counts_tab) <- as.vector(class_tax_tab$class) 

# Tax: make a table of the counts of the Proteobacterial classes
proteo_class_counts_tab <- class_counts_tab[row.names(class_counts_tab) %in% proteo_classes_vec, ] 

# Tax: create new category for Proteobacteria with NA class 
proteo_no_class_annotated_counts <- phyla_and_unidentified_counts_tab[row.names(phyla_and_unidentified_counts_tab) %in% "Proteobacteria", ] - colSums(proteo_class_counts_tab)

# Tax: combine Proteobacteria class and phylum tables
major_taxa_counts_tab <- rbind(temp_major_taxa_counts_tab, proteo_class_counts_tab, "Unresolved_Proteobacteria"=proteo_no_class_annotated_counts)
identical(colSums(major_taxa_counts_tab), colSums(otu_table(ps)))

# Tax: generate a proportions table for summarizing
major_taxa_proportions_tab <- apply(major_taxa_counts_tab, 2, function(x) x/sum(x)*100)
dim(major_taxa_proportions_tab)

# Tax: only keep taxa that make up more than 7% of counts in any sample
temp_filt_major_taxa_proportions_tab <- data.frame(major_taxa_proportions_tab[apply(major_taxa_proportions_tab, 1, max) > 7, ])
dim(temp_filt_major_taxa_proportions_tab) 
filtered_proportions <- colSums(major_taxa_proportions_tab) - colSums(temp_filt_major_taxa_proportions_tab)
filt_major_taxa_proportions_tab <- rbind(temp_filt_major_taxa_proportions_tab, "Other phlya"=filtered_proportions)

# Tax: copy the table so it's safe to manipulate
filt_major_taxa_proportions_tab_for_plot <- filt_major_taxa_proportions_tab

# Tax: add taxa name column so it's not just as R row names 
filt_major_taxa_proportions_tab_for_plot$Major_Taxa <- row.names(filt_major_taxa_proportions_tab_for_plot)

# Tax: transform in narrow/long format and compare to old table
filt_major_taxa_proportions_tab_for_plot.g <- pivot_longer(data= filt_major_taxa_proportions_tab_for_plot, names_to = "sampleID", values_to = "proportion", c(-Major_Taxa))
#filt_major_taxa_proportions_tab_for_plot.g <-  gather(filt_major_taxa_proportions_tab_for_plot, Sample, Proportion, -Major_Taxa)
head(filt_major_taxa_proportions_tab_for_plot.g)
head(filt_major_taxa_proportions_tab_for_plot)

# Tax: make aesthetic table for "color" and "characteristics" of each sample and merge into our plotting table
sample_info_for_merge<-data.frame("sampleID"=row.names(sample_data(ps)), "siteID"=sample_data(ps)$siteID, "source"=sample_data(ps)$source, "year" = sample_data(ps)$year, stringsAsFactors=F)
filt_major_taxa_proportions_tab_for_plot.g2 <- merge(filt_major_taxa_proportions_tab_for_plot.g, sample_info_for_merge)
filt_major_taxa_proportions_tab_for_plot.g2 <- filt_major_taxa_proportions_tab_for_plot.g2  %>% 
mutate(Major_Taxa = factor(Major_Taxa, levels = c("Acidobacteriota","Actinobacteriota", "Alphaproteobacteria","Bacteroidota","Bdellovibrionota","Campylobacterota","Chloroflexi",  "Cyanobacteria","Desulfobacterota","Firmicutes", "Gammaproteobacteria", "Gemmatimonadota",   "Planctomycetota", "Spirochaetota", "Verrucomicrobiota", "Other phlya")))

# Tax: plot individual samples - FH 2021
a <- filt_major_taxa_proportions_tab_for_plot.g2 %>% 
   filter(siteID == "FH" & year == "2021") %>% 
   mutate(source = factor(source, levels = c("water", "sediment"))) %>% 
   mutate(sampleID = str_remove(sampleID, pattern = "S$")) %>% 
  # filter(sampleID == "FH02" | sampleID == "FH03" | sampleID == "FH04" | sampleID == "FH05" | sampleID == "FH06" | sampleID == "FH07" | sampleID == "FH08" | sampleID == "FH02S" | sampleID == "FH03S" | sampleID == "FH04S" | sampleID == "FH05S" | sampleID == "FH06S" | sampleID == "FH07S" | sampleID == "FH08S") %>% 
   # mutate(sampleID = factor(sampleID, levels = c("FH02", "FH03", "FH04", "FH05", "FH06", "FH07" , "FH08" , "FH02S", "FH03S" , "FH04S", "FH05S", "FH06S", "FH07S", "FH08S"))) %>% 
   ggplot(aes(x=sampleID, y=proportion, fill=Major_Taxa)) +
  geom_bar(stat="identity") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1), legend.title=element_blank()) +
  labs(x="Sample", y="% of 16S rRNA gene copies recovered") +
  scale_fill_manual(values = pal) +
   theme_prism(base_line_size = 0.2, base_fontface = "plain", base_size = 12) +
   theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1, size = 9), strip.background = element_rect(fill = "white", colour = "grey"), strip.text.x = element_text(margin = margin (t = 3, b = 5))) +
   labs(x="Samples", y="% of 16S rRNA gene\ncopies recovered") + 
   facet_grid(~ source);a
 
# Tax: plot individual samples - FH 2020
b <- filt_major_taxa_proportions_tab_for_plot.g2 %>% 
  filter(siteID == "FH" & year == "2020") %>% 
  mutate(source = factor(source, levels = c("water", "sediment"))) %>% 
  mutate(sampleID = str_remove(sampleID, pattern = "S$")) %>% 
  ggplot(aes(x=sampleID, y=proportion, fill=Major_Taxa)) +
  geom_bar(stat="identity") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1), legend.title=element_blank()) +
  labs(x="Sample", y="% of 16S rRNA gene copies recovered") +
  scale_fill_manual(values = pal) +
  theme_prism(base_line_size = 0.2, base_fontface = "plain", base_size = 12) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1, size = 9), strip.background = element_rect(fill = "white", colour = "grey"), strip.text.x = element_text(margin = margin (t = 3, b = 5))) +
  labs(x="Samples", y="% of 16S rRNA gene\ncopies recovered") + 
  facet_grid(~ source);b

ggarrange(b, a, labels = c("A", "B"), nrow = 2, ncol = 1, heights = c(1,1), common.legend = TRUE, legend = "right", align = "h")
ggsave(here("figures_thesis", "16S_FHbarchart.png"), height = 8, width = 9)
 
# Tax: plot individual samples - HG 2021 Bdellovibrio
temp <- subset_samples(ps, row.names(sample_data(ps)) == "HG16")
head(sort(otu_table(temp), decreasing = TRUE), n = 50)
temp2 <- row.names(sort(otu_table(temp), decreasing = TRUE))
temp3 <- temp2[1:50]
TAX[rownames(tax_table(TAX)) %in% temp3,]
filt_major_taxa_proportions_tab_for_plot.g2 %>% 
  filter((sampleID == "HG16" | sampleID == "RB04") & Major_Taxa == "Bdellovibrionota")

c <- filt_major_taxa_proportions_tab_for_plot.g2 %>% 
  filter(siteID == "HG" & year == "2021" & source == "water" & sampleID != "HG11" & sampleID != "HG12") %>% 
  mutate(source = factor(source, levels = c("water", "sediment"))) %>% 
  mutate(sampleID = str_remove(sampleID, pattern = "S$")) %>% 
  ggplot(aes(x=sampleID, y=proportion, fill=Major_Taxa)) +
  geom_bar(stat="identity") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1), legend.title=element_blank()) +
  labs(x="Sample", y="% of 16S rRNA gene copies recovered") +
  scale_fill_manual(values = pal) +
  theme_prism(base_line_size = 0.2, base_fontface = "plain", base_size = 12) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1, size = 9), strip.background = element_rect(fill = "white", colour = "grey"), strip.text.x = element_text(margin = margin (t = 3, b = 5))) +
  labs(x="Samples", y="% of 16S rRNA gene\ncopies recovered");c

# Tax: plot individual samples - RB 2020 Bdellovibrio
d <- filt_major_taxa_proportions_tab_for_plot.g2 %>% 
  filter(siteID == "RB" & year == "2020" & source == "water") %>% 
  mutate(source = factor(source, levels = c("water", "sediment"))) %>% 
  mutate(sampleID = str_remove(sampleID, pattern = "S$")) %>% 
  ggplot(aes(x=sampleID, y=proportion, fill=Major_Taxa)) +
  geom_bar(stat="identity") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1), legend.title=element_blank()) +
  labs(x="Sample", y="% of 16S rRNA gene copies recovered") +
  scale_fill_manual(values = pal) +
  theme_prism(base_line_size = 0.2, base_fontface = "plain", base_size = 12) +
  theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1, size = 9), strip.background = element_rect(fill = "white", colour = "grey"), strip.text.x = element_text(margin = margin (t = 3, b = 5))) +
  labs(x="Samples", y="% of 16S rRNA gene\ncopies recovered");d

ggarrange(c, d, labels = c("A", "B"), nrow = 1, ncol = 2, widths = c(1,1), common.legend = TRUE, legend = "right", align = "h")
ggsave(here("figures_thesis", "16S_Bdellobarchart.png"), height = 4.2, width = 9)

# Tax: facet wrap samples by site for sediment
 filt_major_taxa_proportions_tab_for_plot.g2 %>% 
  filter(source == "sediment") %>% 
   filter(!(siteID == "BV")) %>% 
   ggplot(aes(x=sampleID, y=proportion, fill=Major_Taxa)) +
   geom_bar(width=0.8, stat="identity") +
   facet_grid(~ siteID, scales="free_x") +
   scale_fill_manual(values=as.vector(pal)) +
   geom_vline(data= filter(filt_major_taxa_proportions_tab_for_plot.g2, siteID == "FH"), aes(xintercept = 7.5)) +
   geom_vline(data= filter(filt_major_taxa_proportions_tab_for_plot.g2, siteID == "HG"), aes(xintercept = 7.5)) +
   geom_vline(data= filter(filt_major_taxa_proportions_tab_for_plot.g2, siteID == "ML"), aes(xintercept = 7.5)) +
   geom_vline(data= filter(filt_major_taxa_proportions_tab_for_plot.g2, siteID == "RB"), aes(xintercept = 6.5)) +
   geom_vline(data= filter(filt_major_taxa_proportions_tab_for_plot.g2, siteID == "ST"), aes(xintercept = 7.5)) +
   theme_prism(base_line_size = 0.2, base_fontface = "plain", base_size = 12) +
   theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1, size = 9), strip.background = element_rect(fill = "white", colour = "grey"), strip.text.x = element_text(margin = margin (t = 3, b = 5))) +
   labs(x="Sediment samples", y="% of 16S rRNA gene copies recovered") 
 
 ggsave(here("figures_thesis", "16S_sedbarchart.png"), height = 4.2, width = 12)
 
 # Tax: facet wrap samples by site for water
 filt_major_taxa_proportions_tab_for_plot.g2 %>% 
   filter(source == "water") %>% 
   filter(!(siteID == "BV"))  %>% 
   ggplot(aes(x=sampleID, y=proportion, fill=Major_Taxa)) +
   geom_bar(width=0.8, stat="identity") +
   facet_grid(~ siteID, scales="free_x") +
   scale_fill_manual(values=as.vector(pal)) +
   geom_vline(data= filter(filt_major_taxa_proportions_tab_for_plot.g2, siteID == "FH"), aes(xintercept = 8.5)) +
   geom_vline(data= filter(filt_major_taxa_proportions_tab_for_plot.g2, siteID == "HG"), aes(xintercept = 6.5)) +
   geom_vline(data= filter(filt_major_taxa_proportions_tab_for_plot.g2, siteID == "ML"), aes(xintercept = 7.5)) +
   geom_vline(data= filter(filt_major_taxa_proportions_tab_for_plot.g2, siteID == "RB"), aes(xintercept = 6.5)) +
   geom_vline(data= filter(filt_major_taxa_proportions_tab_for_plot.g2, siteID == "ST"), aes(xintercept = 8.5)) +
   theme_prism(base_line_size = 0.2, base_fontface = "plain", base_size = 12) +
   theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1, size = 9), strip.background = element_rect(fill = "white", colour = "grey"), strip.text.x = element_text(margin = margin (t = 3, b = 5))) +
 labs(x="Water samples", y="% of 16S rRNA gene copies recovered") 

 ggsave(here("figures_thesis", "16S_waterbarchart.png"), height = 4.2, width = 12)
 

###############
## PERMANOVA ##
###############
 
# PERMANOVA of water samples, analysis of dissimilarities, but be careful of dispersion problem (if one group is more variable than others, then that can be interpreted as difference in means, rather than difference in variance)
# testing for differences in group means/centroids
 
 # Subset phyloseq object into water and sediment:
ps # main object
ps_wat <- subset_samples(ps, source == "water" & siteID != "BV" & year == "2021" & julianday > 149 & siteID != "FH")
# Prune ASVs not present in at least one sample:
ps_wat <- prune_taxa(taxa_sums(ps_wat) > 0, ps_wat)
 
# Transfer phyloseq object into vegan-compatible format 
df <- as_tibble(sample_data(ps_wat), rownames = "sampleID")
spp <- as.matrix(t(otu_table(ps_wat)))
 
# Hellinger transformation on spp counts (gives low weights to rare species):
#spe.hel <- decostand(chooseTaxa(varespec, n.occ = 5), method = "hellinger") # keeps spp present in at least 5 samples
spph <- decostand(spp, method = "hellinger")
 
# Calculate dissimilarity matrix for betadisper:
dis <- vegdist(spph)
mod <- betadisper(dis, df$cows_present);mod
boxplot(mod)
set.seed(25)
h <- how(within = Within(type = "series", mirror = FALSE), 
          #plots = with(df, Plots(strata = siteID, type = "none")) 
          blocks = df$siteID
 )
permutest(mod, permutations = h, pairwise = TRUE)
plot(mod)
# significant differences in variance are present 
 
# If variable fails betadisper test, cannot do the following adonis (permanova):
permanova <- adonis2(spph ~ cows_present, data = df, permutations = h, by = "margin");permanova
# significant difference in means are also present, along with significant differences in variance

# Sediment samples:
ps_sed <- subset_samples(ps, source == "sediment" & siteID != "BV" & year == "2021" & julianday > 149 & siteID != "FH")
# Prune ASVs not present in at least one sample:
ps_wat <- prune_taxa(taxa_sums(ps_sed) > 0, ps_sed)

# Transfer phyloseq object into vegan-compatible format 
df <- as_tibble(sample_data(ps_sed), rownames = "sampleID")
spp <- as.matrix(t(otu_table(ps_sed)))

# Hellinger transformation on spp counts (gives low weights to rare species):
spph <- decostand(spp, method = "hellinger")

# Calculate dissimilarity matrix for betadisper:
dis <- vegdist(spph)
mod <- betadisper(dis, df$cows_present);mod
boxplot(mod)
set.seed(25)
h <- how(within = Within(type = "series", mirror = FALSE), 
         #plots = with(df, Plots(strata = siteID, type = "none")) 
         blocks = df$siteID
)
permutest(mod, permutations = h, pairwise = TRUE)
plot(mod)
# significant differences in variance are present 

# If variable fails betadisper test, cannot do the following adonis (permanova):
set.seed(25)
permanova <- adonis2(spph ~ cows_present, data = df, permutations = h, by = "margin");permanova
# significant difference in means are also present, along with significant differences in variance

# Plot NMDS
mds <- metaMDS(spph, distance = "bray", k = 3, trymax = 50)
# extract x and y coordinates from MDS plot into new dataframe, so you can plot with ggplot 
z <- data.frame(mds$points)
mds$stress # 0.137
# colour by sample type
ggplot(z, aes(MDS1, MDS2, color=df$siteID)) + geom_point() + theme_bw() + ggtitle('Unconstrained ordination, stress:0.15')

