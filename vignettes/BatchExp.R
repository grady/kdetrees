library(kdetrees)
##library(hgm)
library(BatchExperiments)
##library(devtools)
library(ape)


pmcoa <- function(dynamic,k=1.5,...){
  i <- which(detect.complete.outliers(pMCOA(dynamic,...),k=k)$TFgn)
  nout <- attr(dynamic,"nout")
  hit <- sum(seq_len(nout) %in% i)
  list(hit=hit/nout, type1=(length(i)-hit)/(length(dynamic)-nout))
}


## a "problem"
makeDataset <- function(static,neff,ncoal=100,sp.file){
    coaltrees <- read.tree(text=system2("sim/treesim.py",
                               c("-n",neff,"-s",sp.file,"-N",ncoal),stdout=TRUE),
                           keep.multi=TRUE)
    nout <- 1
    structure(c(.uncompressTipLabel(sample(static,nout)),
              .uncompressTipLabel(coaltrees)),nout=nout)
}

## algorithm
kdesim <- function(dynamic,...){
    kde.obj <- kdetrees(dynamic,...)
    i <- kde.obj$i
    nout <- attr(dynamic,"nout")
    hit <- sum(seq_len(nout) %in% i)
    c(hit=hit/nout,type1=(length(i)-hit)/(length(kde.obj$density)-nout))
}

out.trees <- unname(read.nexus("sim/species1k.nex"))

## reg <- makeExperimentRegistry("rocExp", seed=42, packages=c("kdetrees","ape"))
## addRegistrySourceFiles(reg, "sim/my-pmcoa.R", src.now = TRUE)

addProblem(reg, "coalGen", static=out.trees, dynamic=makeDataset)
addAlgorithm(reg, "kdetrees", kdesim)
addAlgorithm(reg, "pmcoa", pmcoa)

## removeProblem(reg, "coalGen")
## removeAlgorithm(reg,"kdetrees")
## removeAlgorithm(reg,"pmcoa")

## 'X design' = X + parameter sets
## 'experiment' = problem designn + algorithm design + number of reps
sp.files <- c("sim/species.nex","sim/species2.nex","sim/species5.nex")
neffs <- round(exp(seq(log(500),log(3000),len=10)))
neffDesign <- makeDesign("coalGen", exhaustive=list(sp.file=sp.files, neff=neffs))

ks <- seq(-2,3,by=0.125)
methods <- c("geodesic","dissimilarity")
kdeDesign <- makeDesign("kdetrees", exhaustive=list(k=ks,distance=methods))
addExperiments(reg, neffDesign, kdeDesign, repls=200)

pmcDesign <- makeDesign("pmcoa", exhaustive=list(k=ks, distance="patristic"))
addExperiments(reg, neffDesign, pmcDesign, repls=200, skip.defined=TRUE)

## removeExperiments(reg,findExperiments(reg))

testJob(reg)

## restart here
reg <- loadRegistry("rocExp-files")
options(BBmisc.ProgressBar.style="off")
submitJobs(reg, chunk(findNotDone(reg), chunk.size=500))
##submitJobs(reg, chunk(findNotStarted(reg), chunk.size=100))


res <- reduceResultsExperiments(reg,findDone(reg), fun=function(job,res) as.list(res))
save(res,file="sim/res.Rdata")


library(plyr)
res2 <- ddply(res,c("prob","neff","sp.file","algo","distance","k"),summarise,hit=mean(hit),type1=mean(type1))
