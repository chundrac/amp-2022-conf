require(phytools)
require(phangorn)
require(rstan)
require(cmdstanr)

set.seed(1234)

model.type <- commandArgs(trailingOnly = TRUE)[1]
#tree.index <- as.numeric(commandArgs(trailingOnly = TRUE)[2])

mod <- cmdstan_model(paste('models/',model.type,'.stan',sep=''))

data.df <- read.csv('character_data.csv',row.names=1)

data.df <- data.df[,which(apply(data.df,2,function(x){length(unique(x[x!='absent']))})==4)]

#data.df <- data.df[,which(apply(data.df,2,function(x){length(x[x!='absent'])})>30|apply(data.df,2,function(x){length(unique(x[x!='absent']))})==4)]

trees <- read.nexus('austronesian-mapped.nex')

tree <- maxCladeCred(trees)

#data.lame <- data.frame(language=names(data.df.i),feature=data.df.i)
#ggtree(tree.i) %<+% data.lame + geom_tippoint(aes(col=feature)) + geom_tiplab(size=2) + theme(legend.position = c(0.2,0.9))

fit.list <- list()
for (feat.index in 1:ncol(data.df)) {
  data.df.i <- data.df[,feat.index]
  names(data.df.i) <- row.names(data.df)
  #data.df.i <- data.df.i[data.df.i!='absent']
  tree.i <- keep.tip(tree,names(data.df.i))
  tree.i <- reorder.phylo(tree.i,'pruningwise')
  data.binarized <- to.matrix(data.df.i,seq=c("absent","FALSE 0","FALSE 1","TRUE 0","TRUE 1"))
  J <- ncol(data.binarized)
  bin.states <- data.binarized[tree.i$tip.label,]
  bin.states <- rbind(as.matrix(bin.states),matrix(1,nrow=tree.i$Nnode,ncol=ncol(bin.states)))
  parent <- tree.i$edge[,1]
  child <- tree.i$edge[,2]
  b.lens <- tree.i$edge.length/1000
  N <- length(unique(c(parent,child)))
  T <- length(child[which(!child %in% parent)])
  tip.lik <- bin.states
  data.list <- list(N=N,
                    T=T,
                    B=length(parent),
                    brlen=b.lens,
                    child=child,
                    parent=parent,
                    tiplik=tip.lik,
                    J=J)
  #stanfit <- stan(file=paste('models/',model.type,'.stan',sep=''),data=data.list,chains=4)
  fit <- mod$sample(data=data.list,chains=4,parallel_chains=4)
  stanfit <- rstan::read_stan_csv(fit$output_files())
  fit.list[[feat.index]] <- stanfit
  save.image(paste('output/',model.type,'_','mcc','_fit.Rdata',sep=''))
  print(feat.index/ncol(data.df))
}