#QUALITY CONTROL
library(affyPLM) 
PT <- ReadAffy()
PLM <- fitPLM(PT) #This does RMA


#TO ACCESS NUMBER OF GENES: 
expr_matrix <- coefs(PLM)
probeset_names <- rownames(expr_matrix)
n_probesets <- length(unique(probeset_names))

#NUMBER OF CONTROL PROBES
keep_idx <- !grepl("^AFFX", probeset_names)
control_probes <- probeset_names[!keep_idx]
n_control <- length(unique(control_probes))

#CHIPS PSEUDO-IMAGES
par(mfrow = c(6,4))
image(PLM, type = "weights") #If this doesn't work, reinstall affyPLM


#RLE 
par(mfrow = c(1,1))
NUSE(PLM, main = "NUSE of dataset")
abline(h = 1.05, col = "red") #bound for acceptable NUSE
NUSE(PLM, type = "stats")


#PLOTTING THE POSSIBLY BAD CHIPS, AND COMPARING THEM WITH  “”GOOD” CHIPS
par(mfrow = c(3,2))
for (x in c(9,14,15,19,40,42)) {
  image(PLM, which =x, type = "weights")
}
for (x in c(4,12,23,31,34,45)) {
  image(PLM, which =x, type = "weights")
}
image(PLM, which =14, type = "weights", add.legend=TRUE)

colors <- rep("lightgrey", 45)
colors[14] <- "red" # we conclude that the 14th chip in our dataset (Hyperdip-50-C28) has bad quality, so we label it in red
for (i in c(9,15,19,40,42)){
  colors[i] <- "orange" #Suspicious chips

}

#ADDING COLORED LABELS TO THE FAULTY/SUSPICIOUS CHIPS IN THE NUSE
NUSE(PLM, main = "NUSE of dataset", col = colors, xaxt="n")
abline(h = 1.05, col = "red")







#DIFFERENTIAL EXPRESSION
#The following is done without the file Hyperdip-50-C28.CEL in the working directory
library(limma)
PT <- ReadAffy()
PLM <- fitPLM(PT)


#CONSTRUCTING AND FITTING THE LINEAR MODEL
origin <- substr(sampleNames(PLM),1,1) |> factor(levels = c("N","H"))
design <- model.matrix(~origin)  
LM <- lmFit(PLM,design)

#EMPIRICAL BAYES
LM <- eBayes(LM, trend=TRUE, robust=TRUE)

#ADJUSTING P-VALUES (done in toptable()) AND SORTING BY log FC
results <- topTable(LM, number = Inf, sort.by = "logFC")
de_genes <- results[results$adj.P.Val < 0.05, ]

#NUMBER OF  DE GENES 
nrow(de_genes)

#VISUALISATION (VOLCANO PLOT)
M <- LM$coefficients[,2]
neglog10adjP <- -log10(results$adj.P.Val) #-log10 to improve readability
#M <- LM$coefficients[,2]
#B <- LM$lods[,2]

plot(M,neglog10adjP,pch=".", ylab = "-log10(adjP)")
abline(h=-log10(0.05),lty=3)
abline(v=c(-1,1),lty=5)
bigM <- (abs(M) > 1)
smalladjP <- (neglog10adjP > -log10(0.05))

points(M[bigM&!smalladjP],neglog10adjP[bigM&!smalladjP],pch=".",col="green", cex = 3)
points(M[smalladjP&!bigM],neglog10adjP[smalladjP&!bigM],pch=".",col="magenta")
points(M[bigM&smalladjP],neglog10adjP[bigM&smalladjP],pch=".",col="blue", cex = 5)

legend("topright",
       legend = c("adjP = 0.05",
                  "|M| = 1"),
       lty = c(3, 5),
       col = "black",
       bty = "o", #n is possible as well, maybe use ?
       cex = 0.8) 

#TO GET THE 6 BLUE PROBE SETS
probeset_names[bigM&smalladjP]

#TOP50 DE GENES
Top50 <- de_genes[1:50,]
library(xtable)
xtable_res <- xtable(Top50, caption="Top 50 DE genes")
print(xtable_res)




#CLUSTERING (REMOVE NORMAL AND BAD QUALITY CHIPS BEFORE DOING THIS)
library(cluster)
library(affyPLM) 
PT <- ReadAffy()
PLM <- fitPLM(PT) #This does RMA



expr_matrix <- coefs(PLM)

#DISTINGUISHING CONTROL PROBES
probeset_names <- rownames(expr_matrix)
keep_idx <- !grepl("^AFFX", probeset_names)
control_probes <- probeset_names[!keep_idx]


#RANKING PROBES BY VARIANCE
top_var_indices <-  apply(expr_matrix, 1, var) |> order( decreasing = TRUE) 


#CHECKING THE RANK OF CONTROL PROBES
sorted_probes <- rownames(expr_matrix)[top_var_indices]
is_control <- sorted_probes %in% control_probes
control_pos <- which(is_control)
as.table(setNames(control_pos, sorted_probes[control_pos]))

#MDS PLOT TO CHECK FOR BATCH EFFECTS
plotMDS(PLM)

#TAKING OUT CONTROL PROBES
expr_matrix <- expr_matrix[keep_idx, ]


#REDUCING NOISE BY LIMITING Nº OF GENES IN CONSIDERATION 
top_var_indices <-  apply(expr_matrix, 1, var) |> order( decreasing = TRUE) #again, as we took out control probes
n_genes <- 1000 
Tvar_matrix <- expr_matrix[top_var_indices[1:n_genes], ]


#CLUSTERING
samples.cor.ward <-hclust(as.dist(1-cor(Tvar_matrix)), method = "ward.D2")
genes.cor.ward <- hclust (as.dist(1-cor(t(Tvar_matrix))), method = "ward.D2") #needed for later


# SILHOUETTE ANALYSIS
n_samples <- length(samples.cor.ward$labels)
avg_sil <- numeric(n_samples-2)
for (k in 2:(n_samples-1)) {
  clusters <- cutree(samples.cor.ward, k)
  sil <- silhouette(clusters, as.dist(1-cor(Tvar_matrix)))
  avg_sil[k-1] <- mean(sil[, 3], na.rm=TRUE)
}
par(mfrow = c(1,1))
plot(2:(n_samples-1), avg_sil[1:(n_samples-2)], type = "b", xlab = "Number of clusters", 
     ylab = "Average silhouette width", 
     main = "Silhouette Analysis of Hyperdiploid>50 clusterings")

best_three_avsil <- sort(avg_sil, decreasing = TRUE)[1:3]
top3_k <- order(avg_sil, decreasing = TRUE)[1:3] + 1
par(mfrow = c(1,3))
for (k in top3_k) {
  clusters <- cutree(samples.cor.ward, k)
  sil <- silhouette(clusters, as.dist(1-cor(Tvar_matrix)))
  plot(sil, main = paste("Silhouette plot",k-1))
}
data <- matrix(best_three_avsil, ncol=1) 
rownames(data) <- top3_k
colnames(data) <- c('Average Sil.') 
data <- as.table(t(data))

#HEATMAPS FOR CLUSTERS WITH MAX AVG SIL
library(gplots)
for (n_c in top3_k) {
  #Coloring
  P <-rainbow(n_c)
  #coloring <- c("1" = P[1], "2" = P[2], "3" = P[3], "4" = P[4])
  clusters_n_c <- cutree(samples.cor.ward, k = n_c)
  cluster_colors <- P[clusters_n_c]
  heatmap.2(Tvar_matrix,
            Rowv = as.dendrogram(genes.cor.ward),
            Colv = as.dendrogram(samples.cor.ward),
            ColSideColors = cluster_colors,
            trace = "none",           # cleaner look
            density.info = "none",
            col = heat.colors(100),   # or another palette
            key = TRUE,               # show color legend
            key.title = "Legend",
            cexCol = 0.75,
            srtCol = 35,
            labRow = NA,
            key.xlab = "Expression value")

}



