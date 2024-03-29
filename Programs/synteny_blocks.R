library(pacman)
p_load(data.table, rtracklayer, ggplot2, ggthemes, plyranges, ggpubr, plotly,
    BRGenomics, reshape2, dplyr, gplots, Biostrings, scales, gtools, annotables,
    ChIPseeker, TxDb.Mmusculus.UCSC.mm10.knownGene, org.Mm.eg.db, org.Dr.eg.db,
    TxDb.Drerio.UCSC.danRer11.refGene, BSgenome.Drerio.UCSC.danRer11,
    BSgenome.Mmusculus.UCSC.mm10)

################################################################################

# Declaring paths:
PROJECT <- paste0("/home/daniele/Desktop/IV_course/I_semester/",
                    "Kursinis_projektas/Tbx5_analysis_II/")
INPUTS <- paste0(PROJECT, "Inputs/")
INTER_FILES <- paste0(PROJECT, "Intermediate_data_files/")
FIGURES <- paste0(PROJECT, "Figures/")
RESULTS <- paste0(PROJECT, "Peak_annotations/")
GTF_FASTA <- paste0(PROJECT, "Seq_Annot/")

################################################################################

# Reading a file that describes synteny between organisms:
synteny_blocks <- read.csv(paste0(INTER_FILES, "synteny_blocks.csv"), skip = 1)
samples <- read.csv(file = paste0(INTER_FILES, "sample_key.csv"))

# Mouse vs. Zebrafish
# Getting synteny block positions and chromosome number for mouse:
synteny_chr_pos <<- synteny_blocks[, c(7:8, 6, 3:4, 2)]
colnames(synteny_chr_pos) <- c("startpos_MM", "endpos_MM", "chr_MM",
                                "startpos_DR", "endpos_DR", "chr_DR")

# Creating a dataframe that stores data:
df_mm <- data.frame(chr = paste0("chr", c(synteny_chr_pos[, 3])),
                    start = c(synteny_chr_pos[, 1]),
                    end = c(synteny_chr_pos[, 2]), strand = "+")

# Checking whether the starting region position is bigger than end position:
for (i in 1:length(rownames(df_mm))) {
    if (df_mm[i, 2] > df_mm[i, 3]) {
        reverse1 = df_mm[i, 2]
        df_mm[i, 2] = df_mm[i, 3]
        df_mm[i, 3] = reverse1
        df_mm[i, 4] <- "-"
    }
}

# Creating a GRanges object for created dataframe (GRange object for
# synteny blocks):
syn_grange <- makeGRangesFromDataFrame(df_mm)

# Reading bigBed files:
bbfiles <- list.files(path = paste0(INPUTS, "BigBed/"), "*bb")
grl <- GRangesList()

# Declaring what chromosomes are analyzed:
chr_abr <- c(paste0("chr", 1:19), "chrX", "chrY")

# Creating GRanges objects and extracting certain chromosomes:
for (i in 1:length(bbfiles)) {
    grl[[i]] <- import(paste0(paste0(INPUTS, "BigBed/"), bbfiles[i])) %>%
                       dplyr::filter(seqnames %in% chr_abr)
    names(grl)[i] <- samples$Graph_names[samples$Filename == bbfiles[i]]
}

################################################################################

chr_syn <- as.data.frame(table(synteny_blocks[, 6]))
colnames(chr_syn) <- c("Chr", "Freq")
chr_syn$Chr <- sub("^", "Chr", chr_syn$Chr)
chr_order <- c(paste0("Chr", 1:19), "ChrX")
chr_syn$Chr <- factor(chr_syn$Chr, levels = chr_order)

chr_syn <- chr_syn[order(chr_syn$Chr), ]
colours <- ifelse(chr_syn$Freq >= 100, "#b61f1f","black")
sizes <- ifelse(chr_syn$Freq >= 100, 6, 5)
bold <- ifelse(chr_syn$Freq >= 100, 2, 0)

################################################################################

plot0 <- ggplot(chr_syn, aes(x = Chr, y = Freq)) + 
    geom_bar(width = 0.8, size = 0.2, fill = "#be800d", color = "#523a0f",
             stat = 'identity', position = position_dodge(0.5)) +
    geom_text(aes(label = Freq), colour = colours, vjust = -0.2, hjust = 0.5,
              size = sizes, fontface = bold) +
    labs(title = "", x = "", y = "") + 
    ylim(0, 2.5) +
    scale_y_continuous(labels = label_number(suffix = " K", scale = 1e-2)) +
    theme(panel.background = element_rect(fill = "#eeeef1",
                                          colour = "#4c0001"),
          panel.grid.major.y = element_line(colour = "#cab5b5", size = 0.3,
                                            linetype = "dashed"),
          panel.grid.minor.y = element_line(colour = "#cab5b5", size = 0.3,
                                            linetype = "dashed"),
          panel.grid.major.x = element_line(colour = "#cab5b5", size = 0.2,
                                            linetype = "longdash"),
          panel.grid.minor.x = element_line(colour = "#cab5b5", size = 0.2,
                                            linetype = "longdash"),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 11,
                                     face = "bold", color = "black"),
          axis.text.y = element_text(size = 11, face = "bold", color = "black"),
          axis.title.x = element_text(size = 2),
          axis.title.y = element_text(size = 16),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.title = element_text(face = "bold", size = 10),
          legend.text = element_text(size = 10))

png(file = paste0(FIGURES, "Synteny_Chr_MM.png"))
plot0
dev.off()
plot0

################################################################################

chr_syn <- as.data.frame(table(synteny_blocks[, 2]))
colnames(chr_syn) <- c("Chr", "Freq")
chr_syn$Chr <- sub("^", "Chr", chr_syn$Chr)
chr_order <- c(paste0("Chr", 1:25))
chr_syn$Chr <- factor(chr_syn$Chr, levels = chr_order)

chr_syn <- chr_syn[order(chr_syn$Chr), ]
colours <- ifelse(chr_syn$Freq >= 70, "#b61f1f","black")
sizes <- ifelse(chr_syn$Freq >= 70, 6, 5)
bold <- ifelse(chr_syn$Freq >= 70, 2, 0)

plot1 <- ggplot(chr_syn, aes(x = Chr, y = Freq)) + 
    geom_bar(width = 0.8, size = 0.2, fill = "#611d02", color = "#310f01",
             stat = 'identity', position = position_dodge(0.5)) +
    geom_text(aes(label = Freq), colour = colours, vjust = -0.2, hjust = 0.5,
              size = sizes, fontface = bold) +
    labs(title = "", x = "", y = "") + 
    ylim(0, 2.5) +
    scale_y_continuous(labels = label_number(suffix = " K", scale = 1e-2)) +
    theme(panel.background = element_rect(fill = "#eeeef1",
                                          colour = "#4c0001"),
          panel.grid.major.y = element_line(colour = "#cab5b5", size = 0.3,
                                            linetype = "dashed"),
          panel.grid.minor.y = element_line(colour = "#cab5b5", size = 0.3,
                                            linetype = "dashed"),
          panel.grid.major.x = element_line(colour = "#cab5b5", size = 0.2,
                                            linetype = "longdash"),
          panel.grid.minor.x = element_line(colour = "#cab5b5", size = 0.2,
                                            linetype = "longdash"),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,
                                     size = 11, face = "bold",
                                     color = "black"),
          axis.text.y = element_text(size = 11, face = "bold",
                                     color = "black"),
          axis.title.x = element_text(size = 2),
          axis.title.y = element_text(size = 16),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.title = element_text(face = "bold", size = 10),
          legend.text = element_text(size = 10))

png(file = paste0(FIGURES, "Synteny_Chr_DR.png"))
plot1
dev.off()
plot1

################################################################################

plot1 <- ggplot(chr_syn, aes(x = Chr, y = Freq)) + 
    geom_bar(width = 0.8, size = 0.2, fill = "#611d02", color = "#310f01",
             stat = 'identity', position = position_dodge(0.5)) +
    geom_text(aes(label = Freq), colour = colours, vjust = -0.2, hjust = 0.5,
              size = sizes, fontface = bold) +
    labs(title = "", x = "", y = "") + 
    ylim(0, 2.5) +
    scale_y_continuous(labels = label_number(suffix = " K", scale = 1e-2)) +
    theme(panel.background = element_rect(fill = "#eeeef1",
                                          colour = "#4c0001"),
          panel.grid.major.y = element_line(colour = "#cab5b5", size = 0.3,
                                            linetype = "dashed"),
          panel.grid.minor.y = element_line(colour = "#cab5b5", size = 0.3,
                                            linetype = "dashed"),
          panel.grid.major.x = element_line(colour = "#cab5b5", size = 0.2,
                                            linetype = "longdash"),
          panel.grid.minor.x = element_line(colour = "#cab5b5", size = 0.2,
                                            linetype = "longdash"),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,
                                     size = 11, face = "bold",
                                     color = "black"),
          axis.text.y = element_text(size = 11, face = "bold",
                                     color = "black"),
          axis.title.x = element_text(size = 2),
          axis.title.y = element_text(size = 16),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.title = element_text(face = "bold", size = 10),
          legend.text = element_text(size = 10))

png(file = paste0(FIGURES, "Synteny_Chr_DR.png"))
plot1
dev.off()
plot1

################################################################################

indicators <- c("synteny_overlaps", "peak_count")
percentages <- c()
grl_reduced <- GRangesList()

# Creating an empty dataframe and naming it's columns:
df_final <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(df_final) <- c("names", "indic", "value", "percentage")

# Determining whether overlaps between synteny blocks and peaks exist:
for (i in 1:length(grl)) {
    total_peaks = length(grl[[i]])
    overlap_count = length(subsetByOverlaps(grl[[i]], syn_grange))
    values <- c(overlap_count, total_peaks)

    for (y in 1:length(indicators)) {
        row = c(samples$Graph_names[i], indicators[y], values[y],
                (overlap_count / total_peaks) * 100)
        df_final[nrow(df_final) + 1, ] <- row
    }

    grl_reduced[[i]] <- subsetByOverlaps(grl[[i]], syn_grange)
    names(grl_reduced)[i] <- names(grl[i])
}

################################################################################

plot2 <- ggplot(df_final, aes(fill = indic, y = as.numeric(value), x = names)) + 
    geom_bar(width = 1, size = 0.2, colour = "#3f2704", stat = 'identity',
             position = position_dodge(0.8)) +
    geom_text(aes(label = ifelse(indic == "synteny_overlaps",
                                 paste0(round(as.numeric(percentage),
                                              digits = 2), "%"), ""),
              fontface = 2), vjust = -1.2, hjust = -0.5, size = 5) +
    labs(title = "", x = "", y = "") + 
    scale_fill_manual(values = c("#d37f02", "#f3bc47"),
                      labels = c("Pikų skaičius", "Synteny blokų skaičius")) +
    scale_y_continuous(labels = label_number(suffix = " K", scale = 1e-3)) +
    guides(fill = guide_legend(title = "Spalvų paaiškinimas")) +
    theme(panel.background = element_rect(fill = "#eeeef1",
                                          colour = "#4c0001"),
          panel.grid.major.y = element_line(colour = "#cab5b5", size = 0.3,
                                            linetype = "dashed"),
          panel.grid.minor.y = element_line(colour = "#cab5b5", size = 0.3,
                                            linetype = "dashed"),
          panel.grid.major.x = element_line(colour = "#cab5b5", size = 0.2,
                                            linetype = "longdash"),
          panel.grid.minor.x = element_line(colour = "#cab5b5", size = 0.2,
                                            linetype = "longdash"),
          axis.text.x = element_text(angle = 0, hjust = 1, vjust = 1, size = 11,
                                     face = "bold", color = "black"),
          axis.text.y = element_text(size = 11, face = "bold",
                                     color = "black"),
          axis.title.x = element_text(size = 2),
          axis.title.y = element_text(size = 16),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.title = element_text(face = "bold", size = 10),
          legend.text = element_text(size = 10)) +
    coord_flip()

png(file = paste0(FIGURES, "Synteny_blocks_peaks.png"))
plot2
dev.off()
plot2

################################################################################

mm_known_genes <- TxDb.Mmusculus.UCSC.mm10.knownGene
grl_synteny <- GRangesList()

# Annotating Mus musculus peaks that belong to synteny block:
for (object in 1:length(grl_reduced)) {
    name <- gsub(" ", "_", names(grl_reduced[object]))
    peak <- grl_reduced[[object]]
    seqlengths(peak) <- seqlengths(peak) - 10001
    peak_annotation <- annotatePeak(peak, tssRegion = c(-3000, 3000),
                         TxDb = mm_known_genes, annoDb = "org.Mm.eg.db")

    print(peak_annotation)
    
    mm_annot <- as.data.frame(peak_annotation@anno)
    print(length(rownames(mm_annot)))
    entrezids <- unique(mm_annot$geneId)
    entrez2gene <- grcm38 %>% filter(entrez %in% entrezids) %>%
                                dplyr::select(entrez, symbol)

    m <- match(mm_annot$geneId, entrez2gene$entrez)
    mm_annot <- cbind(mm_annot[, 1:14], gene_symbol = entrez2gene$symbol[m],
                        mm_annot[, 15:16])

    # Defining the same grl_smlr object that has two extra columns with
    # gene id and gene symbol:
    grl_synteny[[object]] <- mm_annot
    names(grl_synteny)[object] <- names(grl)[object]
}

write.table(grl_synteny, file = paste0(INTER_FILES, "/grl_synteny.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

################################################################################

genes_samples <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(genes_samples) <- c("Sample", "Gene_count", "Percentage")

for (smpl in 1:length(grl_synteny)) {
    genes_tot <- length(unique(grl_synteny[[smpl]]$gene_symbol))

    row = c(names(grl_synteny[smpl]),
            genes_tot, (genes_tot / length(genes(mm_known_genes))) * 100)
    genes_samples[nrow(genes_samples) + 1, ] <- row
}

# Diagram that shows how many unique genes were identified for each
# sample belonging to the synteny block:
plot3 <- ggplot(genes_samples, aes(x = Sample, y = as.numeric(Percentage))) +
    geom_bar(stat = "identity", fill = "#cc8b12",
             width = .5, color = "#8d5c00") +
    labs(title = "", x = "", y = "Anotuotų pikų % M. musculus") +
    geom_text(aes(label = paste0(round(as.numeric(Percentage), digit = 2), "%")),
              color = "#030101", size = 5, vjust = -1, fontface = 2) +
    coord_cartesian(ylim = c(0, as.numeric(max(genes_samples$Percentage)) + 20)) +
    theme(axis.text = element_text(size = 10, colour = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 12,
                                     vjust = 1, face = "bold"),
          axis.text.y = element_text(size = 12, face = "bold"),
          axis.title.x = element_text(size = 14, colour = "black"),
          axis.title.y = element_text(size = 14, colour = "black"),
          panel.grid.major = element_line(color = "#eeeeee"),
          plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
          panel.background = element_rect(fill = "#eeeef1", colour = "#3a1010"),
          panel.grid.major.y = element_line(colour = "#cab5b5",
                                            size = 0.3, linetype = "dashed"),
          panel.grid.minor.y = element_line(colour = "#cab5b5",
                                            size = 0.3, linetype = "dashed"),
          panel.grid.major.x = element_line(colour = "#cab5b5",
                                            size = 0.2, linetype = "longdash"),
          panel.grid.minor.x = element_line(colour = "#cab5b5",
                                            size = 0.2, linetype = "longdash"),
          legend.position = c(0.777, 0.788))

png(file = paste0(FIGURES, "Unique_genes_MM.png"))
plot3
dev.off()
plot3

################################################################################

# Getting positions of Danio rerio genes:
gene_ranges <- genes(TxDb.Drerio.UCSC.danRer11.refGene)

# Adding 'gene_symbol' column to gene_ranges GRange:
gene_ranges$gene_symbol <- mapIds(org.Dr.eg.db, keys = gene_ranges$gene_id,
                column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")

# Getting synteny block positions and chromosome number for zebrafish:
fish_chr_pos <- synteny_blocks[, c(3:4, 2)]
colnames(fish_chr_pos) <- c("startpos", "endpos", "chr")

# Creating a dataframe that stores data:
df_dr <- data.frame(chr = paste0("chr", c(fish_chr_pos[, 3])),
                    start = c(fish_chr_pos[, 1]),
                    end = c(fish_chr_pos[, 2]), strand = "+")

for (i in 1:length(rownames(df_dr))) {
    if (df_dr[i, 2] > df_dr[i, 3]) {
        reverse1 = df_dr[i, 2]
        df_dr[i, 2] = df_dr[i, 3]
        df_dr[i, 3] = reverse1
        df_dr[i, 4] <- "-"
    }
}

# Creating a GRanges object for created dataframe:
syn_grange_dr <- makeGRangesFromDataFrame(df_dr)

dr_genes <- GRangesList()

# Getting the genes that overlap with Mus musculus genes:
for (smpl2 in 1:length(grl_synteny)) {
    selected_genes <- tolower(gene_ranges$gene_symbol) %in%
                            tolower(unique(grl_synteny[[smpl2]]$gene_symbol))
    genes <- gene_ranges[selected_genes]

    dr_genes[[smpl2]] <- genes
    names(dr_genes)[smpl2] <- names(grl_synteny[smpl2])
}

# Creating an empty dataframe to store data about gene counts:
genes_synteny <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(genes_synteny) <- c("Sample", "Gene_count", "Percentage")

# Finding genes that belong to synteny block:
for (smpl3 in 1:length(grl_synteny)) {
    common_genes <- subsetByOverlaps(dr_genes[[smpl3]], syn_grange_dr)
    total_gene_length <- length(tolower(gene_ranges$gene_symbol))
    row = c(names(grl_synteny[smpl3]), length(common_genes),
            (length(common_genes) / total_gene_length) * 100)
    genes_synteny[nrow(genes_synteny) + 1, ] <- row
}

################################################################################

plot4 <- ggplot(genes_synteny, aes(x = Sample, y = as.numeric(Percentage))) +
    geom_bar(stat = "identity", fill = "#611d02",
            width = .5, color = "#310f01") +
    labs(title = "", x = "", y = "Genų % D. rerio") +
    geom_text(aes(label = paste0(round(as.numeric(Percentage), digit = 2), "%")),
              color = "#030101", size = 5, vjust = -1, fontface = 2) +
    coord_cartesian(ylim = c(0, as.numeric(max(genes_synteny$Percentage)) + 3)) +
    theme(axis.text = element_text(size = 10, colour = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 12,
                                     vjust = 1, face = "bold"),
          axis.text.y = element_text(size = 12, face = "bold"),
          axis.title.x = element_text(size = 14, colour = "black"),
          axis.title.y = element_text(size = 14, colour = "black"),
          panel.grid.major = element_line(color = "#eeeeee"),
          plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
          panel.background = element_rect(fill = "#eeeef1", colour = "#3a1010"),
          panel.grid.major.y = element_line(colour = "#cab5b5",
                                            size = 0.3, linetype = "dashed"),
          panel.grid.minor.y = element_line(colour = "#cab5b5",
                                            size = 0.3, linetype = "dashed"),
          panel.grid.major.x = element_line(colour = "#cab5b5",
                                            size = 0.2, linetype = "longdash"),
          panel.grid.minor.x = element_line(colour = "#cab5b5",
                                            size = 0.2, linetype = "longdash"),
          legend.position = c(0.777, 0.788))

png(file = paste0(FIGURES, "Unique_genes_DR.png"))
plot4
dev.off()
plot4

################################################################################

grl_mm <- GRangesList()
for (object in 1:length(grl)) {
    name <- gsub(" ", "_", names(grl[object]))
    peak <- grl[[object]]
    seqlengths(peak) <- seqlengths(peak) - 10001
    peak_annotation <- annotatePeak(peak, tssRegion = c(-3000, 3000),
                         TxDb = mm_known_genes, annoDb = "org.Mm.eg.db")

    mm_annot <- as.data.frame(peak_annotation@anno)
    entrezids <- unique(mm_annot$geneId)
    entrez2gene <- grcm38 %>% filter(entrez %in% entrezids) %>%
                                dplyr::select(entrez, symbol)

    m <- match(mm_annot$geneId, entrez2gene$entrez)
    mm_annot <- cbind(mm_annot[, 1:14], gene_symbol = entrez2gene$symbol[m],
                      mm_annot[, 15:16])

    # Defining the same grl_smlr object that has two extra columns with
    # gene id and gene symbol:
    grl_mm[[object]] <- mm_annot
    names(grl_mm)[object] <- names(grl)[object]
}

################################################################################

# Defining a dataframe that stores values describing the percentage of
# genes that belong to synteny blocks:
belong_block <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(belong_block) <- c("Sample", "Block_count", "Percentage")

find_common_genes <- function(sample, name_title) {
    count <- 0
    sample <- grl_mm[[1]]
    name_title <- names(grl_mm[1])
    mm <- sample       # Mus musculus
    dr <- gene_ranges       # Danio rerio
    # name_title <- names(grl_mm[y])
    synteny_block_n <- length(rownames(synteny_chr_pos))

    df_mm <- as.data.frame(na.omit(mm)) # d
    df_dr <- as.data.frame(na.omit(dr)) # c

    # Retrieving data from certain Mus musculus GRange (as dataframe) columns:
    df_mm <- df_mm[, c("seqnames", "start", "end", "width", "strand",
                       "annotation", "gene_symbol")]
    
    # Creating a dataframe that will store data about common genes (chromosome,
    # start and end positions, organism name (either Mus musculus
    # or Danio rerio)):
    common_genes <- data.frame(matrix(ncol = 5, nrow = 0))
    colnames(common_genes) <- c("Gene", "Chr", "Start", "End", "Organism")

    ########################## FINDING COMMON GENES ############################

    # Loop that iterates every synteny block (every row from dataset that
    # describes syntenic regions between M.musculus and D.rerio) and determines
    # common genes:
    for (i in 1:synteny_block_n) {
        syn_row <- synteny_chr_pos[i, ]

        # Extracing M.musculus rows that match the conditions (match
        # chromosome numbers and belong to synteny block interval):
        genes_mm <- na.omit(df_mm[df_mm$seqnames == paste0("chr",
                                                           syn_row$"chr_MM") &
                                  df_mm$start >= syn_row$"startpos_MM" &
                                  df_mm$end <= syn_row$"endpos_MM", ])

        # Extracing D.rerio rows that match the conditions (match chromosome
        # numbers and belong to synteny block interval):
        genes_dr <- na.omit(df_dr[df_dr$seqnames == paste0("chr",
                                                           syn_row$"chr_DR") &
                                  df_dr$start >= syn_row$"startpos_DR" &
                                  df_dr$end <= syn_row$"endpos_DR", ])
        
        # Finding common genes:
        common_g <- intersect(tolower(genes_mm[,"gene_symbol"]),
                              tolower(genes_dr[,"gene_symbol"]))

        # Condition that forms a new dataframe if genes belong to synteny
        # block and they are common to both M.musculus and D.rerio organisms:
        if (length(rownames(genes_dr)) != 0 & length(rownames(genes_mm)) != 0 &
            length(common_g) != 0) {
            print(paste0("i value = ", i))
            count = count + 1 
            print(paste0("count = ", count))

            # Data about common genes in D.rerio:
            cg_data_dr <- genes_dr[genes_dr$"gene_symbol" %in% common_g,
                                   c("seqnames", "start", "end", "gene_symbol")]
            cg_data_dr[, ncol(cg_data_dr) + 1] <- rep("Danio_rerio",
                                                      length(nrow(cg_data_dr)))
            colnames(cg_data_dr)[5] <- "organism"

            # Data about common genes in M.musculus:
            cg_data_mm <- genes_mm[tolower(genes_mm$"gene_symbol") %in% common_g,
                                   c("seqnames", "start", "end", "annotation",
                                     "gene_symbol")]
            cg_data_mm$"gene_symbol" <- tolower(cg_data_mm$"gene_symbol")
            cg_data_mm[, ncol(cg_data_mm) + 1] <- rep("Mus_musculus",
                                                      length(nrow(cg_data_mm)))
            colnames(cg_data_mm)[6] <- "organism"

            genes <- rbind(cg_data_dr, cg_data_mm[, c("seqnames", "start", "end",
                                                    "gene_symbol", "organism")])
            common_genes <- rbind(common_genes, genes)
        }
    }

    # write.table(belong_block,
    #             file = paste0(INTER_FILES, "/Belong_to_block.txt"),
    #             sep = "\t", quote = FALSE, row.names = FALSE)

    ################### GETTING SEQUENCES FOR COMMON GENES #####################

    # Creating a sequence dataset that stores M.musculus gene sequences:
    gene_mm <- common_genes[common_genes$"organism" == "Mus_musculus", ]
    mm_gr <- df_mm[df_mm$start %in% gene_mm$"start" &
                   df_mm$end %in% gene_mm$"end" &
                   df_mm$seqnames %in% gene_mm$"seqnames",
                   c("seqnames", "start", "end", "gene_symbol")]

    mm_gr <- makeGRangesFromDataFrame(mm_gr, keep.extra.columns = TRUE)
    seq_mm <- getSeq(BSgenome.Mmusculus.UCSC.mm10, mm_gr)
    names(seq_mm) <- mm_gr$gene_symbol

    # Creating a sequence dataset that stores D.rerio gene sequences:
    gene_dr <- common_genes[common_genes$"organism" == "Danio_rerio", ]
    dr_gr <<- df_dr[df_dr$start %in% gene_dr$"start" &
                   df_dr$end %in% gene_dr$"end" &
                   df_dr$seqnames %in% gene_dr$"seqnames",
                   c("seqnames", "start", "end", "gene_symbol", "gene_id")]

    dr_gr <<- makeGRangesFromDataFrame(dr_gr, keep.extra.columns = TRUE)
    seq_dr <- getSeq(BSgenome.Drerio.UCSC.danRer11, dr_gr)
    names(seq_dr) <- dr_gr$gene_symbol

    ####################### COUNTING HITS OF TBX5 PWM ##########################

    mpwm <- read.table(paste0(INTER_FILES, "TBX5_MOUSE.H11MO.0.D.pwm"))
    mpwm <- t(mpwm)
    rownames(mpwm) <- c("A", "C", "G", "T")

    # Counting PWM matches for Danio rerio:
    gene_hits_dr <- data.frame(matrix(ncol = 3, nrow = 0))
    colnames(gene_hits_dr) <- c("Gene", "Hits", "Length")

    # Determining the number of Tbx5 binding hits:
    for (gene in 1:length(seq_dr)) {
        hits <- countPWM(as.matrix(mpwm), as.character(seq_dr[gene]),
                         min.score = "75%")
        row <- c(names(seq_dr[gene]), hits, nchar(as.character(seq_dr[gene])))
        gene_hits_dr[nrow(gene_hits_dr) + 1, ] <- row
    }
    gene_hits_dr <- as.data.frame(gene_hits_dr[order(gene_hits_dr$Gene), ])
    gene_hits_dr[, ncol(gene_hits_dr) + 1] <- rep("Danio rerio",
                                                  length(nrow(gene_hits_dr)))

    # Counting PWM matches for Mus musculus:
    gene_hits_mm <- data.frame(matrix(ncol = 3, nrow = 0))
    colnames(gene_hits_mm) <- c("Gene", "Hits", "Length")

    # Determining the number of Tbx5 binding hits:
    for (gene in 1:length(seq_mm)) {
        hits <- countPWM(as.matrix(mpwm), as.character(seq_mm[gene]),
                         min.score = "75%")
        row <- c(names(seq_mm[gene]), hits, nchar(as.character(seq_mm[gene])))
        gene_hits_mm[nrow(gene_hits_mm) + 1, ] <- row
    }
    gene_hits_mm <- gene_hits_mm[order(gene_hits_mm$Gene), ]

    # Summing the values of the repeating genes in Mus musculus:
    gene_hits_mm <- data.table(gene_hits_mm)
    gene_hits_mm <- as.data.frame(
                        gene_hits_mm[, list(Hits = sum(as.numeric(Hits)),
                                            Length = sum(as.numeric(Length))),
                                            by = 'Gene'])

    gene_hits_mm[, ncol(gene_hits_mm) + 1] <- rep("Mus musculus",
                                                  length(nrow(gene_hits_mm)))

    # Creating a combined dataframe that stores Tbx5 hit values of
    # M.musculus and D.rerio organisms:
    combined_df <- rbind(gene_hits_dr, gene_hits_mm)
    colnames(combined_df)[4] = "Organism"

    plot5 <- ggplot(combined_df, aes(x = as.numeric(Hits),
                                     y = as.numeric(Length),
                                     colour = Organism)) +
        geom_point(alpha = 0.7, size = 5, shape = 20) +
        scale_colour_manual(values = c("#301c05", "#b96517")) +
        labs(colour = "Organizmas", title = name_title, x = "", y = "Ilgis") +
        scale_x_continuous(limits = c(0, 80)) +
        coord_cartesian(ylim = c(0, 400000)) +
        scale_y_continuous(labels = label_number(suffix = " K", scale = 1e-3)) +
        geom_vline(xintercept = 20, linetype = "dashed") + 
        geom_hline(yintercept = 100000, linetype = "dashed") +
        theme(panel.background = element_rect(fill = "#eeeef1",
                                              colour = "#4c0001"),
            panel.grid.major.y = element_line(colour = "#cab5b5", size = 0.3,
                                              linetype = "dashed"),
            panel.grid.minor.y = element_line(colour = "#cab5b5", size = 0.3,
                                              linetype = "dashed"),
            panel.grid.major.x = element_line(colour = "#cab5b5", size = 0.2,
                                              linetype = "longdash"),
            panel.grid.minor.x = element_line(colour = "#cab5b5", size = 0.2,
                                              linetype = "longdash"),
            axis.text.x = element_text(angle = 0, hjust = 1, vjust = 1,
                                       size = 12, face = "bold",
                                       color = "black"),
            axis.text.y = element_text(size = 12, face = "bold",
                                       color = "black"),
            axis.title.x = element_text(size = 16),
            axis.title.y = element_text(size = 16),
            plot.title = element_text(hjust = 0.5, face = "bold"),
            legend.title = element_text(face = "bold", size = 10),
            legend.text = element_text(size = 10, color = "black"))
    # png(file = paste0(FIGURES, "PWM_matches.png"))
    # plot5
    # dev.off()
    plot5
}

################################################################################

# Creating a list that stores ggplot objects:
plots <- list()
for (smpl in 1:length(grl_mm)) {

    remove_y <- theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank())
    
        plots[[smpl]] <- find_common_genes(grl_mm[[smpl]],
                                           names(grl_mm[smpl])) + remove_y
}

png(file = paste0(FIGURES, "PWM_matches_all.png"), width = 1200)
ggarrange(plotlist = plots, nrow = 1, common.legend = TRUE, widths = 5,
          legend = "bottom")
dev.off()

################################################################################

# library(ReactomePA)
# library(clusterProfiler)

# ggo <- groupGO(gene = (dr_gr)$gene_id, "org.Dr.eg.db", ont = "BP", level = 2,
#                readable = TRUE)
# head(summary(ggo))
# barplot(ggo, drop = TRUE, showCategory = 20)