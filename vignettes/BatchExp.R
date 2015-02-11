library(kdetrees)
library(hgm)
library(BatchExperiments)
library(devtools)
library(ape)



## a "problem"
makeDataset <- function(static,neff,ncoal=100,sp.file){
    coaltrees <- read.tree(text=system2("sim/treesim.py",
                               c("-n",neff,"-s",sp.file,"-N",ncoal),stdout=TRUE),
                           keep.multi=TRUE)
    nout <- 1
    structure(c(sample(static,nout),coaltrees),nout=nout)
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

reg <- makeExperimentRegistry("rocExp", seed=42, packages=c("kdetrees","ape"))
addProblem(reg, "coalGen", static=out.trees, dynamic=makeDataset)
addAlgorithm(reg, "kdetrees", kdesim)

removeProblem(reg, "coalGen")
removeAlgorithm(reg,"kdetrees")

## 'X design' = X + parameter sets
## 'experiment' = problem designn + algorithm design + number of reps
sp.files <- c("sim/species.nex","sim/species2.nex","sim/species5.nex")
neffs <- round(exp(seq(log(500),log(3000),len=10)))
neffDesign <- makeDesign("coalGen", exhaustive=list(sp.file=sp.files, neff=neffs))

ks <- seq(-2,3,by=0.125)
methods <- c("geodesic","dissimilarity")
rocDesign <- makeDesign("kdetrees", exhaustive=list(k=ks,distance=methods))

addExperiments(reg, neffDesign, rocDesign, repls=100)

## removeExperiments(reg,findExperiments(reg))

testJob(reg)

dd <- makeDataset(out.trees,neff=1000,ncoal=20,"sim/species.nex")

kdesim(dd)
kdetrees(dd)

