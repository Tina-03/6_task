
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("BiocStyle", "dada2", "phyloseq", "DECIPHER"))
install.packages(c("knitr","ggplot2","gridExtra","phangorn"))
install.packages("gridExtra")
install.packages("phangorn")
library("knitr")
library("BiocStyle")
library("ggplot2")
library("gridExtra")
library("dada2")
library("phyloseq")
library("DECIPHER")
library("phangorn")

# устанавливаем рабочую директорию с fastq файлами
miseq_path <- "C:\\Users\\tanya\\Downloads\\Таня\\NGS\\мит\\MiSeq_SOP\\MiSeq_SOP"

list.files(miseq_path)

# сортируем по названию файлы с прямыми и обратными ридами
fnFs <- sort(list.files(miseq_path, pattern="_R1_001.fastq"))
fnRs <- sort(list.files(miseq_path, pattern="_R2_001.fastq"))
fnFs[1:3]
# извлекаем из названий файлов названия образцов
sampleNames <- sapply(strsplit(fnFs, "_"), `[`, 1)

# сохраняем в переменные полные пути к файлам с ридами
fnFs <- file.path(miseq_path, fnFs)
fnRs <- file.path(miseq_path, fnRs)
fnFs[1:3]
fnRs[1:3]

# визуализация качества секвенирования в зависимости от позиции
plotQualityProfile(fnFs[1:4]) # для прямых ридов : зеленые средние значения, обрезаем 10 с конца
plotQualityProfile(fnRs[1:4]) # для обратных ридов :конец хуже из-за технологии - обрезаем более 90

# создаем директорию для отфильтрованных данных
filt_path <- file.path(miseq_path, "filtered")
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sampleNames, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sampleNames, "_R_filt.fastq.gz"))

# фильтруем и обрезаем риды - аналог триммоматика, возьми прямые риды и сохрани в прямые отфильтрованные, транклен - то, что хотим оставить
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE
head(out)

# находим уникальные риды - дедупликация
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

names(derepFs) <- sampleNames
names(derepRs) <- sampleNames

# оценка параметров модели для коррекции ошибок секвенирования
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF)
plotErrors(errR)

# коррекция ошибок секвенирования
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
dadaFs[[1]]

# объединяем прямые и обратные риды
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs)

# сводная таблица с количеством одинаковых последовательностей в образцах
seqtabAll <- makeSequenceTable(mergers[!grepl("Mock", names(mergers))])

# сводная таблица с длинами последовательностей - выводятся самые часто встречающиеся
table(nchar(getSequences(seqtabAll)))

# идентификация и удаление химерных последовательностей
seqtabNoC <- removeBimeraDenovo(seqtabAll)



#### далее оцениваем принадлежность последовательностей известным организмам

# делаем прогноз для последовательностей с использованием наивного байеса и обучающей выборки
fastaRef <- "C:\\Users\\tanya\\Downloads\\Таня\\NGS\\мит\\rdp_train_set_16.fa.gz"

taxTab <- assignTaxonomy(seqtabNoC, refFasta = fastaRef, multithread=TRUE)
unname(head(taxTab))

# выполняем множественное выравнивание
seqs <- getSequences(seqtabNoC)
names(seqs) <- seqs
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA,verbose=FALSE)

# строим филогенетическое дерево (neighbor-joining tree)
phangAlign <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phangAlign)
treeNJ <- NJ(dm)
plot(treeNJ,cex=0.1)

# строим филогенетическое дерево (GTR+G+I maximum likelihood tree) оч долго
fit <- pml(treeNJ, data=phangAlign)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))

plot(fitGTR,cex=0.1)

# импортируем информацию по образцам - СЛЕД ЭТАП ВЫПОЛНЯТЬ НЕ ТРЕБУЕТСЯ
#НУЖЕН только список бактерий, далее интерпретация на геномах
samdf <- read.csv("MIMARKS_Data_combined.csv",row.names=1)

# создаем объект класса phyloseq
ps <- phyloseq(otu_table(seqtabNoC, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxTab),phy_tree(fitGTR$tree))
save(ps,file = "ps.RData")
