ddw <- DEWSeq::DESeqDataSetFromSlidingWindows(countData  = as.data.frame(countData),
                                      colData    = colData,
                                      annotObj   = as.data.frame(annotationData),
                                      tidy       = TRUE,
                                      design     = ~type+experiment+protein+treatment)

ddw <- DESeq2::estimateSizeFactors(ddw)
DESeq2::sizeFactors(ddw)



countData[which(greater_10>0),]

colSums(countData[which(greater_10>0),][,-1])



count_thresholds <- c(5, 10, 15, 20)
scaling_mat <- list()
N_vec <- c()

for(i in seq_along(count_thresholds)){

    greater_10 <- apply(counts(ddw) >= count_thresholds[i], 1, prod)

    ddw_greater10 <- ddw[which(greater_10>0), ]
    ddw_greater10_mRNAs <- ddw_greater10#[ rowData(ddw_greater10)[,"gene_type"] == "protein_coding", ]


    iqr_samples <- apply(counts(ddw_greater10_mRNAs), 2, function(col_elem){
            quantile(col_elem, c(0.25, 0.75))
        })

    iqr_within <- apply(counts(ddw_greater10_mRNAs), 1, function(row_elem){

           all(row_elem >= iqr_samples[1,] & row_elem <= iqr_samples[2,])
      })


    N_vec[i] <- sum(iqr_within)

    print(head(counts(ddw_greater10_mRNAs[iqr_within, ])))

    ddw_greater10_mRNAs <- ddw_greater10_mRNAs[iqr_within, ]
    ddw_greater10_mRNAs <- estimateSizeFactors(ddw_greater10_mRNAs)

    scaling_mat[[i]] <- as.matrix(DESeq2::sizeFactors(ddw_greater10_mRNAs))

}

i <- length(count_thresholds) + 1

countData_bu <- countData


countData[,-1] <- apply(countData[,-1], 2, function(col_el)col_el + 1)

ddw <- DEWSeq::DESeqDataSetFromSlidingWindows(countData  = as.data.frame(countData[,-1]),
                                      colData    = colData,
                                      annotObj   = as.data.frame(annotationData),
                                      tidy       = TRUE,
                                      design     = ~type+experiment+protein+treatment)

ddw_mRNAs <- ddw[ rowData(ddw)[,"gene_type"] == "protein_coding", ]


iqr_samples <- apply(counts(ddw_mRNAs), 2, function(col_elem){
            quantile(col_elem, c(0.25, 0.75))
        })

iqr_within <- apply(counts(ddw_mRNAs), 1, function(row_elem){

           all(row_elem >= iqr_samples[1,] & row_elem <= iqr_samples[2,])
      })


N_vec[i] <- sum(iqr_within)

ddw_mRNAs <- ddw_mRNAs[iqr_within, ]
ddw_mRNAs <- estimateSizeFactors(ddw_mRNAs, )

scaling_mat[[i]] <- as.matrix(DESeq2::sizeFactors(ddw_mRNAs))

