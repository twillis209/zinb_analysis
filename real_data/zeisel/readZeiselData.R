# Using Risso code to read in Zeisel data which is stored in a rather complicated text file
data <- read.table("expression_mRNA_17-Aug-2014.txt", sep='\t', stringsAsFactors = FALSE, comment.char = '%')

tissue <- as.factor(as.matrix(data)[1,-(1:2)])
group <- as.factor(as.matrix(data)[2,-(1:2)])
nmolecule <- as.numeric(as.matrix(data)[3,-(1:2)])
well <- as.factor(as.matrix(data)[4,-(1:2)])
sex <- as.factor(as.matrix(data)[5,-(1:2)])
age <- as.numeric(as.matrix(data)[6,-(1:2)])
diameter <- as.numeric(as.matrix(data)[7,-(1:2)])
cell_id <- as.matrix(data)[8,-(1:2)]
batch <- as.factor(sapply(strsplit(cell_id, "_"), function(x) x[1]))
position <- as.factor(sapply(strsplit(cell_id, "_"), function(x) x[2]))

level1 <- as.factor(as.matrix(data)[9,-(1:2)])
level2 <- as.factor(as.matrix(data)[10,-(1:2)])

counts <- as.matrix(data[12:NROW(data),-(1:2)])
counts <- matrix(as.numeric(counts), ncol=ncol(counts), nrow=nrow(counts))
rownames(counts) <- data[12:NROW(data),1]
colnames(counts) <- data[8, -(1:2)]
write.csv(counts, 'zeiselCounts.csv', row.names=T)
