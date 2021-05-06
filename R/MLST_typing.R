#' type MLST of a draft genome
#'
#' You provide the path of the draft genome and
#'
#' @param reference the reference sequence that you want to screen. Fasta file in one or severa sequences
#' @param querry the querry sequence. Fasta file in one or severa sequences
#' @param dir.out the directory of the output
#'
#' @return numeric value of the percentage of the reference sequence which is covered by the querry
#' @import GenomicRanges IRanges Biostrings
#'
#' @export





typemlst <- function(genomePath,type)
{

  if(type=='ecoli-warwick')
    {
    allelesPath <- list.files(system.file("extdata/typing/MSLT-warwick/alleles/",package = "EcoliGenomics"),full.names=T)
    STpath      <- list.files(system.file("extdata/typing/MSLT-warwick/profiles/",package = "EcoliGenomics"),full.names=T)
    }
  if(type=='plasmid-IncACcgPMLST')
    {
    allelesPath <- list.files(system.file("extdata/typing/plasmid/IncACcgPMLST/alleles/",package = "EcoliGenomics"),full.names=T)
    STpath      <- list.files(system.file("extdata/typing/plasmid/IncACcgPMLST/profiles/",package = "EcoliGenomics"),full.names=T)
    }
  if(type=='plasmid-IncACPMLST')
    {
    allelesPath <- list.files(system.file("extdata/typing/plasmid/IncACPMLST/alleles/",package = "EcoliGenomics"),full.names=T)
    STpath      <- list.files(system.file("extdata/typing/plasmid/IncACPMLST/profiles/",package = "EcoliGenomics"),full.names=T)
    }
  if(type=='plasmid-IncFRST')
    {
    allelesPath <- list.files(system.file("extdata/typing/plasmid/IncFRST/alleles/",package = "EcoliGenomics"),full.names=T)
    STpath      <- list.files(system.file("extdata/typing/plasmid/IncFRST/profiles/",package = "EcoliGenomics"),full.names=T)
    }
  if(type=='plasmid-IncHI1MLST')
    {
    allelesPath <- list.files(system.file("extdata/typing/plasmid/IncHI1MLST/alleles/",package = "EcoliGenomics"),full.names=T)
    STpath      <- list.files(system.file("extdata/typing/plasmid/IncHI1MLST/profiles/",package = "EcoliGenomics"),full.names=T)
    }

  if(type=='plasmid-IncHI2DLST')
    {
    allelesPath <- list.files(system.file("extdata/typing/plasmid/IncHI2DLST/alleles/",package = "EcoliGenomics"),full.names=T)
    STpath      <- list.files(system.file("extdata/typing/plasmid/IncHI2DLST/profiles/",package = "EcoliGenomics"),full.names=T)
    }

  if(type=='plasmid-IncI1MLST')
    {
    allelesPath <- list.files(system.file("extdata/typing/plasmid/IncI1MLST/alleles/",package = "EcoliGenomics"),full.names=T)
    STpath      <- list.files(system.file("extdata/typing/plasmid/IncI1MLST/profiles/",package = "EcoliGenomics"),full.names=T)
    }

  if(type=='plasmid-IncNMLST')
    {
    allelesPath <- list.files(system.file("extdata/typing/plasmid/IncNMLST/alleles/",package = "EcoliGenomics"),full.names=T)
    STpath      <- list.files(system.file("extdata/typing/plasmid/IncNMLST/profiles/",package = "EcoliGenomics"),full.names=T)
    }


  if(exists('allelesPath')==F)
  {cat('specify one of the following type: ecoli-warwick plasmid-IncACcgPMLST plasmid-IncACPMLST plasmid-IncFRST plasmid-IncHI1MLST plasmid-IncHI2DLST plasmid-IncI1MLST plasmid-IncNMLST')  }

  else
  {
    try(unlink("temp", recursive=TRUE))
    dir.create('temp')
    dir.create('temp/dbblast')
    myarg <-paste0('-in ',genomePath,' -out temp/dbblast/db -dbtype nucl')
    system2(command = 'makeblastdb', args = myarg,stdout=F)

    Nalleles <- length(allelesPath)
    result <- character()

    for (i in 1:Nalleles)
    {
      myarg <-  paste0('-query ',allelesPath[i],' -db temp/dbblast/db -out temp/blast.txt -num_threads 8 -num_alignments 10 -outfmt "7 qacc qlen length pident qstart qend sacc sstart send "' )
      system2(command = 'blastn', args = myarg)

      blast <- try(read.table('temp/blast.txt', comment.char = '#'),silent=T)
      if(class(blast)=='data.frame')
      {
        colnames(blast) <- c('querry.access','querry.length','alignment.lenght','pc.ident.','querry.start','querry.end','subject.access','subject.start','subject.end')
        blast <- blast[(blast$alignment.lenght==blast$querry.length)&(blast$pc.ident.==100),]
        result[i] <- paste0('',as.character(blast$querry.access[which.max(blast$querry.length)]))
      }
      else{result[i] <- ''}
    }

    ST.combi <- read.table(STpath,sep='\t',header=T)
    ST.combi <- ST.combi[duplicated(ST.combi)==F,]
    resultnum <- substr(result,start=unlist(lapply(gregexpr("[0123456789]",result), function(l) l[[1]])),stop=nchar(result))


    names(resultnum) <- gsub(pattern='.fasta',replacement='',x=basename(allelesPath))

    if(min(nchar(resultnum))>0)
    {
      idx <- numeric()
      for(i in 1:Nalleles){idx[i] <- grep(toupper(gsub(pattern='.fasta',replacement='',x=basename(allelesPath))[i]),toupper(colnames(ST.combi)))}
      ST.combiReduced <- ST.combi[,idx]
      rownames(ST.combiReduced) <- ST.combi$ST
      ST <- paste0('',rownames(ST.combiReduced[apply(mapply("==",as.numeric(resultnum),ST.combiReduced),1,sum)==Nalleles,]))
    }
    else{ST<-''}
    names(ST) <- 'ST'

    result <- c(resultnum,ST)
    return(result)

  }
}



