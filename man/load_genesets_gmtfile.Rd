% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/load_genesets_gmtfile.R
\name{load_genesets_gmtfile}
\alias{load_genesets_gmtfile}
\title{parse genesets in GMT format where gene identifiers are numeric Entrez gene IDs}
\usage{
load_genesets_gmtfile(filename, label)
}
\arguments{
\item{filename}{input file for this function should be the full path to genesets defined in GMT format}

\item{label}{a shortname for the genesets in this file, for example "GO_CC", "KEGG", "MY_DB_V1". This will be stored in the 'source' column of the resulting table. Importantly, multiple testing correction in GOAT is grouped by this 'source' column so you probably want to use a different label for each collection-of-genesets that you load. Must not be empty, only allowed characters are; upper/lower-case letter, numbers 0-9 and underscore}
}
\value{
tibble with columns; source (character), source_version (character), id (character), name (character), genes (list), ngenes (int)
}
\description{
parse genesets in GMT format where gene identifiers are numeric Entrez gene IDs
}
\section{Example data;}{
URL: https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp#C5
download this data: KEGG subset of curated pathways -->> NCBI (Entrez) Gene IDs
filename should be something like "c2.cp.kegg.v2023.1.Hs.entrez.gmt"
}

\examples{
  # TODO: update the filename to your downloaded file
  f = "C:/DATA/c2.cp.kegg.v2023.1.Hs.entrez.gmt"
  if(file.exists(f)) genesets_asis = load_genesets_gmtfile(f, label = "KEGG")
}
