library(kdetrees)
##library(hgm)
library(BatchExperiments)
##library(devtools)
library(distory)
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
    kde.obj <- kdetrees(dynamic,...,bw=list(prop=0.19))
    i <- kde.obj$i
    nout <- attr(dynamic,"nout")
    hit <- sum(seq_len(nout) %in% i)
    c(hit=hit/nout,type1=(length(i)-hit)/(length(kde.obj$density)-nout))
}

out.trees <- unname(read.nexus("sim/species1k.nex"))

reg <- makeExperimentRegistry("rocExp", seed=42, packages=c("kdetrees","ape"))
addRegistrySourceFiles(reg, "sim/my-pmcoa.R", src.now = TRUE)

addProblem(reg, "coalescent", static=out.trees, dynamic=makeDataset)
addAlgorithm(reg, "kdetrees", kdesim)
addAlgorithm(reg, "pmcoa", pmcoa)

## removeProblem(reg, "coalGen")
## removeAlgorithm(reg,"kdetrees")
## removeAlgorithm(reg,"pmcoa")

## 'X design' = X + parameter sets
## 'experiment' = problem designn + algorithm design + number of reps
sp.files <- c("sim/species.nex","sim/species2.nex","sim/species5.nex")
neffs <- round(exp(seq(log(500),log(3000),len=4)))
rocDesign <- makeDesign("coalescent", exhaustive=list(sp.file=sp.files, neff=neffs))

ks <- seq(-2,3,by=0.125)
methods <- c("geodesic","dissimilarity")
kdeDesign <- makeDesign("kdetrees", exhaustive=list(k=ks,distance=methods))
addExperiments(reg, rocDesign, kdeDesign, repls=1)#000)

pmcDists <- c("patristic","nodal")
pmcDesign <- makeDesign("pmcoa", exhaustive=list(k=ks[ks<2 & ks>-1.3], distance=pmcDists))
addExperiments(reg, rocDesign, pmcDesign, repls=1)#000)


## removeExperiments(reg,findExperiments(reg))

tt <- makeDataset(out.trees, 500,sp.file="sim/species.nex")
kdetrees(tt)
dm <- as.matrix(dist.multiPhylo(tt))
bw <- kdetrees:::bw.nn(dm,.18)
kdetrees:::bhv.consts(tt,bw)
kdetrees:::bhv.orthant.lb(tt[[1]],bw[1])
testJob(reg)
## restart here
reg <- loadRegistry("rocExp-files")
options(BBmisc.ProgressBar.style="off")
submitJobs(reg, chunk(findNotDone(reg), chunk.size=100))

submitJobs(reg, chunk(findNotStarted(reg), chunk.size=300))
waitForJobs(reg)
showStatus(reg)
res <- reduceResultsExperiments(reg, fun=function(job,res) as.list(res))
save(res,file="sim/res2.Rdata")
