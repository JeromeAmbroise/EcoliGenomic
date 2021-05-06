#' Screen within a reference sequence the presence of a querry sequence
#'
#' With this function, the user provides a reference sequence and a querry sequence and the function compute the percentage of the reference which is covered by the querry. This computation is computed using the GenomicRanges Bioconductor objects.
#'
#'
#' @param reference The path of the reference sequence that you want to screen. The sequence should be a .fasta file including one or several sequences
#' @param querry The path of the querry sequence. he sequence should be a .fasta file including one or several sequences
#' @param dir.out the directory of the output
#'
#' @return numeric value of the percentage of the reference sequence which is covered by the querry sequence.
#' @import GenomicRanges IRanges Biostrings
#'
#' @export

screen_Blast <- function (reference, querry,min.pc.ident,dir.out)
{
  library(Biostrings)
  try(unlink("temp", recursive = TRUE))
  dir.create("temp")
  dir.create("temp/dbblast")
  myarg <- paste0("-in ", reference, " -out temp/dbblast/db -dbtype nucl")
  system2(command = "makeblastdb", args = myarg, stdout = F)
  myarg <- paste0("-query ", querry, " -db temp/dbblast/db -out temp/blast.txt -num_threads 8 -num_alignments 10 -outfmt \"7 qacc bitscore qlen length pident qstart qend sacc sstart send \"")
  system2(command = "blastn", args = myarg)
  blast <- try(read.table("temp/blast.txt", comment.char = "#"), silent = T)
  if (class(blast) == "data.frame")
  {
    colnames(blast) <- c("querry.access", "bitscore", "querry.length", "alignment.lenght", "pc.ident.",
                         "querry.start", "querry.end", "subject.access", "subject.start", "subject.end")
    blast <- blast[blast$pc.ident.>min.pc.ident,]
    start <- blast$subject.start
    end <- blast$subject.end
    new.start <- start
    new.end <- end
    new.start[start>end] <- end[start>end]
    new.end[start>end] <- start[start>end]
    data.frame(start,end,new.start,new.end)
    GR <- GRanges(seqnames = blast$subject.access,ranges = IRanges(start =new.start ,end = new.end))
    GR.disjoin <- disjoin(GR)

    #sequences.levels <- levels(factor(seqnames(GR.disjoin)))
    sequences <- readDNAStringSet(reference)
    sequences.levels <- names(sequences)

    hitlength <- numeric()
    seqlength <- numeric()
    pc.identity <- numeric()
    for(i in seq_along(sequences.levels))
    {
      hitlength[i] <- sum(width(GR.disjoin[seqnames(GR.disjoin)==sequences.levels[i]]))
      pc.identity[i] <- mean(blast$pc.ident.[blast$subject.access==sequences.levels[i]])
      seqlength[i] <- width(sequences[names(sequences)==sequences.levels[i]])
    }

    pc.covered <- round(100*(hitlength/seqlength),1)
    result <- data.frame(sequences.levels,seqlength,hitlength,pc.identity,pc.covered)
    result <- na.omit(result)

    resname <- basename(reference)
    resname <- unlist(strsplit(resname,split='_'))
    resname <- resname[length(resname)]
    resname <- gsub(resname,pattern = '.fasta',replacement = '')
    write.csv(result,paste0(dir.out,'/',resname,'_',gsub(basename(querry),pattern = '.fasta',replacement = ''),'.csv'))
    percentage <- max(pc.covered)


  }
  else{percentage <- 0}
  unlink('temp',recursive = T)
  return(percentage)
}
