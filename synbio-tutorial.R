library(readr)

GLS <- read_delim("/Users/nachonase/Documents/GitHub/Synbio-tutorial/GLS.txt", 
                  "\t", escape_double = FALSE, trim_ws = TRUE, show_col_types = FALSE)

GSC <- read_delim("/Users/nachonase/Documents/GitHub/Synbio-tutorial/GSC.txt", 
                  "\t", escape_double = FALSE, trim_ws = TRUE, show_col_types = FALSE)

gsc_nachon <- GSC
names(gsc_nachon) <- c("s", "g")
gsc_nachon <- gsc_nachon[, c(2,1)]

gls_nachon <- GLS
gsc_nachon$g <- as.character(gsc_nachon$g)
gsc_nachon$s <- as.character(gsc_nachon$s)
library(piano)
gsc_nachon_obj <- loadGSC(gsc_nachon)
pvals <- gls_nachon$p
names(pvals) <- gls_nachon$g
directions <- gls_nachon$FC
names(directions) <- gls_nachon$g
nachon_data <- list(gsc = gsc_nachon_obj, pvals = pvals, directions = directions)
gsares <- runGSA(geneLevelStats= pvals , directions=directions,
                 gsc=gsc_nachon_obj, nPerm=500,
                 geneSetStat = "reporter")

networkPlot2(gsares, class="distinct", direction="both", adjusted = TRUE,
             significance = 0.3, geneSets = NULL, lay = "visNetwork")
dev.new()
hm <- GSAheatmap(gsares)
#dev.off()
#writeFilesForKiwi(gsares,"kiwiSynbio")
#3 files are located in Documents!!

#View(GSAsummaryTable(gsares))
GSAsummaryTable(gsares, save=T, file="gsaresSynbio.txt")

