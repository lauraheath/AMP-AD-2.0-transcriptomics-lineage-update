# This script runs GSEA enrichment using branch-specific pseudotime effect sizes to rank genes and plots the results


# setup -------------------------------------------------------------------

BiocManager::install("fgsea")

# libraries
library(synapser)
library(fgsea)
library(tidyverse)



# plotting theme
theme_set(theme_bw())


#This function quantifies the number of significantly enriched GO terms by AD biological domain
bd.tally <- function( enrVct, biodomDefTbl){
  bdt <- bind_cols(
    domain = unique(biodomDefTbl$Biodomain),
    n_term = map_dbl( unique(biodomDefTbl$Biodomain),
                      ~ biodomDefTbl %>% filter(Biodomain == .x) %>% 
                        dplyr::select(GOterm_Name) %>% distinct() %>% nrow()),
    n_sig_term = map_dbl( unique(biodomDefTbl$Biodomain),
                          ~ enrVct[ enrVct %in% 
                                      biodomDefTbl$GOterm_Name[
                                        biodomDefTbl$Biodomain == .x]] %>% 
                            length()) ) %>% 
    mutate(domain = fct_reorder(domain, n_sig_term, .desc = F)) %>% 
    arrange(domain)
  bdt$prop <- bdt$n_sig_term / bdt$n_term
  return(bdt)
}


# data: biodomain definitions, biodomain plotting colors, and braak state-gene table
biodom <- readRDS( synapser::synGet('syn25428992')$path  )
biodom.annotated <- biodom %>% filter(!is.na(n_hgncSymbol)) %>% pull(hgnc_symbol, name=GOterm_Name)
dom.cols <- read_csv( synGet('syn26856828')$path )

#female degs by state
pt_data <- read_csv( synGet('syn39989123')$path )
#pt_data <- read_csv(file="female_DEanova_statsALLGENES.csv")
#OR
#male degs by state
pt_data <- read_csv( synGet('syn39990047')$path )
#pt_data <- read_csv(file="male_DEanova_statsALLGENES.csv")


# Run enrichments by pseudotime state -------------------------------------


# enrichment for each state; all significant genes (pval < 0.05), arranged by effect size
enr <- map_dfr(
  unique(pt_data$state),
  ~ pt_data %>% 
    filter(state == .x
           , pvalue < 0.05
           # , !(gene_names %in% (dlpfc.pt %>% filter(pvalue <= 0.05) %>% pull(gene_names))) # use this to consider genes unique to the neuropath
    ) %>% #
    arrange(desc(effect)) %>% 
    pull(effect, name = gene_names) %>%
    fgseaMultilevel(biodom.annotated, ., eps = 0, scoreType = 'std') %>% 
    mutate(state = .x)
)


# count up the number of significant terms per biodomain
bdt <- bd.tally(enr$pathway[enr$pval <= 0.05], biodom) %>% 
  mutate(domain = fct_reorder(domain, n_sig_term, .desc = T)) %>% 
  arrange(domain)


# add biodomain annotations and plotting colors to enrichment results
enr <- enr %>% 
  left_join(., biodom %>% dplyr::select(pathway=GOterm_Name, Biodomain), by='pathway') %>%
  full_join(., dom.cols, by = c('Biodomain'='domain')) %>%
  mutate( Biodomain = fct_relevel(Biodomain, as.character(bdt$domain) )) %>% 
  mutate( state = fct_reorder( as.factor(state), as.numeric(state), .desc = T ))


# plot! -------------------------------------------------------------------

#tiff(file='FEMALE_biodomains.tiff',height=200,width=200,units='mm',res=300)
tiff(file='MALE_biodomains.tiff',height=200,width=200,units='mm',res=300)

enr %>%   
  ggplot(aes( factor(state), NES ))+ 
  scale_color_identity()+ scale_fill_identity()+
  scale_size_continuous(
    limits = -log10(enr$pval) %>% range(),
    range = c(.7,5))+
  geom_jitter(
    data = subset(enr, pval > 0.01),
    color = 'grey85', alpha = .2, size = .5, show.legend = F
  )+
  geom_jitter(
    data = subset(enr, pval < 0.01 ),
    aes(color = color, size = -log10(pval) ), #
    alpha = .5, show.legend = T
  )+
  geom_violin(
    data = subset(enr, pval < 0.01 & NES < 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_violin(
    data = subset(enr, pval < 0.01 & NES > 0),
    scale='width', aes(fill = color), color = 'grey50',#, color = col
    alpha = .3, show.legend = F
  )+
  geom_hline(yintercept = 0, lty = 3)+
  facet_wrap(~Biodomain)+
  labs(x='')+  coord_flip() + theme(legend.position = 'right') +
  #ggtitle('Female transcriptomics, biodomains in each state vs State1; state pval < 0.01; term pval < 0.01')
  ggtitle('Male transcriptomics, biodomains in each state vs State1; state pval < 0.01; term pval < 0.01')
dev.off()


enr$leadingEdge <- as.character(enr$leadingEdge)
enr <- as.data.frame(enr)
enr2 <- subset(enr, enr$padj<0.1)
#move leadingEdge column to end of dataframe
leadingEdge <- as.data.frame(enr2$leadingEdge)
enr3 <- subset(enr2, select = -c(leadingEdge))
enr3 <- cbind(enr3, leadingEdge)

# enr2$leadingEdge2 <- gsub("\"", "", enr2$leadingEdge)
# enr2$leadingEdge<-NULL
# #find extra tabs
# library(stringr)
# enr2$leadingEdge <- str_replace(enr2$leadingEdge2, "\",\r\n", "\",") 


write.table(enr3, file="data_objects/F_GOenrichment.txt", sep = "\t", row.names=FALSE)
file <- synapser::File(path='data_objects/F_GOenrichment.txt', parentId='syn38349639')
file <- synapser::synStore(file)

file <- synapser::File(path='FEMALE_biodomains.tiff', parentId='syn38347670')
file <- synapser::synStore(file)



write.table(enr3, file="data_objects/M_GOenrichment.txt", sep = "\t", row.names=FALSE)
file <- synapser::File(path='data_objects/M_GOenrichment.txt', parentId='syn38349639')
file <- synapser::synStore(file)

file <- synapser::File(path='MALE_biodomains.tiff', parentId='syn38347670')
file <- synapser::synStore(file)


