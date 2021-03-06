% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MLST_typing.R
\name{typemlst}
\alias{typemlst}
\title{MLST typing of a draft genome}
\usage{
typemlst(genomePath, type)
}
\arguments{
\item{genomePath}{The path of the genome (complete or WGS) of the reference sequence that you want to screen. The sequence should be a .fasta file including one or several sequences}

\item{type}{The name of the MLST scheme to be applied. List of possible values includes: ecoli-warwick plasmid-IncACcgPMLST plasmid-IncACPMLST plasmid-IncFRST plasmid-IncHI1MLST plasmid-IncHI2DLST plasmid-IncI1MLST plasmid-IncNMLST}
}
\value{
numeric vector including the Allele Type and Sequence Type of the genome
}
\description{
With this function, the user provides a genome (complete or WGS) as a fasta file and the function returns the MLST profile of the genome. Currently, severala MLST scheme are provided for MLST typing of E. coli (e.g Warwick scheme and plasmid incompatibility group (IncFRST,..)
}
