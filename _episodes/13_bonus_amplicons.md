---
source: Rmd
title: "Basic plotting,	clustering and mapping amplicon data"
teaching: 0
exercises: 40
questions:
- "How to visualise alpha- and bet-diversity from a set of amplicon sequences?"
- "How to determine which genotypes (OTUs or zOTUs) are contributing to the differences?"
objectives:
- "Organise abundance and classification data in a useful way"
- "Begin to explore the diversity and biogeography of a set of environmental samples"
keypoints:
- "This lesson shows several useful methods from baseR. This process is being ported to Tidyverse and [Phyloseq](https://joey711.github.io/phyloseq/)"
---

## Getting Sarted


~~~
library(vegan)
library(clustsig)
library(mvabund)

#To Do
library(tidyverse)
library(DESeq2)
~~~
{: .language-r}


Load the otu table, environmental data and taxonomy (and tidy the data using gather, join and spread functions)

~~~
otu.table<-read_csv()
taxonomy<-read_csv()
env.data<-read_csv()

~~~
{: .language-r}

> ## Get to know your data first.
>
> Use your cheatsheets to explore and summarise the data
>
> ~~~
> dim(otu.table)
> str(otu.table)
> ~~~
> {: .language-r}
>
> > ## Solution
> >
> >
> > ~~~
> > [1] "This new template looks good"
> > ~~~
> > {: .output}
> {: .solution}
{: .challenge}

## nMDS clustering and simprof

Once the data is cleaned the first *classic* step to consider is to rarefy the data so that there is equal sampling effort at each site. This approach is now recognised as *old school* because multivariate approaches, including model-based multivariate statistics, are now possible with the available computational power. For some background explanation visit Professor David Warton's [Youtube video](https://www.youtube.com/watch?v=KnPkH6d89l4) on his R package [mvabund](http://environmentalcomputing.net/introduction-to-mvabund/). There is a treat for everyone who makes it to the end of the video.

> ## Transpose the table before rarefying
>  because vegan expects sites in rows, species as columns
{: .callout}


Oldshool rarefication of your data to the lowest number of observations in one sample


1. Show lowest 100 number of obs per site in matrix. In practice this is the point where some samples with poor sequence returns can be dropped from the analyses

~~~
sort(colSums(ss2012agg))[1:10]

ss2012aggt<-t(ss2012agg)

ss2012_rare<-rrarefy(ss2012aggt, min(colSums(ss2012agg))
~~~
{: .language-r}

## The rarefied matrix can now be used for simprof and metaMDS

~~~
ss2012nmds<-metaMDS(ss2012_rare, distance="bray")
~~~
{: .language-r}

~~~
ss2012_simprof<-simprof(ss2012_rare, num.expected=100, num.simulated=99, method.cluster="ward.D", method.distance="braycurtis",method.transform="identity", alpha=0.001, sample.orientation="row", const=0, silent=TRUE, increment=100, undef.zero=TRUE, warn.braycurtis=TRUE)
~~~
{: .language-r}

## plot the tree and check significant clusters)

~~~
summary(ss2012_simprof)
~~~
{: .language-r}

~~~
simprof.plot(ss2012_simprof)
~~~
{: .language-r}

### draw a line through branches to pick a rational number of clutsers


~~~
abline(h=220,b=NULL)
~~~
{: .language-r}

### cut the tree at the desired level to define clusters

~~~
mycl <- cutree(simprof_v1$hclust, h=220)
length(unique(mycl))
~~~
{: .language-r}


### add clusters to names d.f

~~~
names$cluster<-mycl
~~~
{: .language-r}


you can give these clusters a defined name by recoding the number


~~~
names$province<-mycl
names$province[names_2$province==2]<-'ASTS1'
~~~
{: .language-r}

Now plot the nmds

~~~
ss2012nmds<-metaMDS(ss2012r, distance="bray")

nmds_cols<-col.cat(mycl_2)
symbols<-c(16,17,18,19) ## code for different symbols when there are, e.g 5 clusters
~~~
{: .language-r}


~~~
ordiplot(example_NMDS,type="n")
orditorp(bobnmds,display="species",col="red",air=0.01)
orditorp(bobnmds,display="sites",cex=1.25,air=0.01)
mtext("2D Stress 0.", 3)
~~~
{: .language-r}


~~~
plot(ss2012nmds, type="n")
  points(ss2012nmds, display="sites", col=nmds_cols, pch=symbols[mycl])
  legend("topright", legend=unique(names_2$province), col=unique(nmds_cols), pch=c(15,16,17,18,19), bty="n")
~~~
{: .language-r}


###now do envfit and add the significant vectors. Need to redo the simprof

##subset the env data and matrix to those sites present in the both

env_sub<-env[(names_2$shortname  %in% row.names(env)),]
ss2012_sub<-(ss2012r[rownames(ss2012r) %in% row.names(env),])

ss2012_sub_nmds<-metaMDS(ss2012_sub, distance='bray')

ss2012_env<-envfit(ss2012_sub_nmds, env_sub)

plot(ss2012_sub_nmds, type="n")
  points(ss2012_sub_nmds, display="sites", col=nmds_cols[env_sub_merge$mycl], pch=symbols[mycl])
  legend("topright", legend=unique(env_sub_merge$province), col=unique(nmds_cols)[env_sub_merge$mycl], pch=c(15,16,17,18,19), bty="n")



###plot on the map after merging codes_namesa and env_sub


plot( env_sub_merge$lon, env_sub_merge$lat, pch=symbols[env_sub_merge$mycl], col=unique(nmds_cols)[env_sub_merge$mycl], cex=2, ylim=c(17,9))
legend("bottomleft", legend=unique(env_sub_merge$province), col=unique(nmds_cols), pch=c(16,17,18,19), bty="n")


#############

##boxplots to check differences in values


 par(mfrow=c(2,3))
  boxplot(tapply(env_sub_merge$Syn, env_sub_merge$province, summary), col=col_o, main="Syn Counts")
  boxplot(tapply(env_sub_merge$PPE, env_sub_merge$province, summary), col=col_o, main= "Plastidic Protists")
  boxplot(tapply(env_sub_merge$NITRATE, env_sub_merge$province, summary), col=col_o, main="Nitrate")
  boxplot(tapply(env_sub_merge$AMMONIA, env_sub_merge$province, summary), col=col_o, main="Ammonia")
  boxplot(tapply(env_sub_merge$TEMP, env_sub_merge$province, summary), col=col_o, main="˚C")
  boxplot(tapply(env_sub_merge$SILICATE, env_sub_merge$province, summary), col=col_o, main="Silicate")

##script for colour categories

col.cat <- function(x=NULL){
   if(is.null(x)) {
     return(NULL)
   }
   v <- unique(x)
   ret <- niceColors(length(v))
   names(ret) <- as.character(v)
   ret[as.character(x)]

##nice Colours 1

niceColors<-function(n=80){
cols<-rep(c("#75954F","#D455E9","#E34423","#4CAAE1","#451431","#5DE737","#DC9B94","#DC3788","#E0A732","#67D4C1","#5F75E2","#1A3125","#65E689","#A8313C","#8D6F96","#5F3819","#D8CFE4","#BDE640","#DAD799","#D981DD","#61AD34","#B8784B","#892870","#445662","#493670","#3CA374","#E56C7F","#5F978F","#BAE684","#DB732A","#7148A8","#867927","#918C68","#98A730","#DDA5D2","#456C9C","#2B5024","#E4D742","#D3CAB6","#946661","#9B66E3","#AA3BA2","#A98FE1","#9AD3E8","#5F8FE0","#DF3565","#D5AC81","#6AE4AE","#652326","#575640","#2D6659","#26294A","#DA66AB","#E24849","#4A58A3","#9F3A59","#71E764","#CF7A99","#3B7A24","#AA9FA9","#DD39C0","#604458","#C7C568","#98A6DA","#DDAB5F","#96341B","#AED9A8","#55DBE7","#57B15C","#B9E0D5","#638294","#D16F5E","#504E1A","#342724","#64916A","#975EA8","#9D641E","#59A2BB","#7A3660","#64C32A"),
ceiling(n/80))
cols[1:n]
}

##nice colours 2

niceColors<-function(n=80){
  cols<-rep(c("#1F497D", "#438AE0", "#A3C5F0", "#991515", "#B81C1C", "#D61E1E", "#F74D4D", "#F56C6C", "#F98B8B", "#FCB1B1", "#FACFCF", "#FD4626", "#FCEA00", "#FFF681", "#FFFAB9", "#4BACC6", "#98E4EC", "#CCC9FF", "#A8007C", "#E200A7", "#FF5DD5", "#6E009C", "#9433BC", "#AB60CB", "#885AFF", "#AD8DFF", "#08A07C", "#08A07C", "#08A07C", "#FF9900", "#FFD058", "#3A3A3A", "#878787", "#BCBCBC", "#EBEBEB", "#0BE8C2", "#C7CF00", "#DCE171", "#A98C78", "#60BA22", "#99D359", "#582900", "#FDF5AF", "#FF9F00", "#fa94a3", "#f44f4f", "#99E990", "#3EB43E", "#008E00"),
  ceiling(n/80))
  cols[1:n]
}

my_env<-c("SS2012t07CTD1d5_MS247", "SS2012t07UW1d5_MS247", "SS2012t07UW2d5_MS247", "SS2012t07UW3d5_MS247", "SS2012t07UW4d5_MS247", "SS2012t07UW5d5_MS247", "SS2012t07UW14d5_MS247", "SS2012t07CTD2d5_MS247", "SS2012t07UW15d5_MS247", "SS2012t07UW16d5_MS247", "SS2012t07UW18d5_MS247", "SS2012t07UW20d5_MS247", "SS2012t07UW22d5_MS247", "SS2012t07UW23d5_MS247", "SS2012t07UW24d5_MS247", "SS2012t07UW25d5_MS247", "SS2012t07UW26d5_MS247", "SS2012t07CTD3d5_MS247", "SS2012t07UW27d5_MS247", "SS2012t07UW29d5_MS247", "SS2012t07CTD4d5_MS247", "SS2012t07CTD5d5_MS247", "SS2012t07UW31d5_MS247", "SS2012t07UW32d5_MS247", "SS2012t07UW33d5_MS247", "SS2012t07UW34d5_MS247", "SS2012t07CTD6d5_MS247", "SS2012t07UW35d5_MS247", "SS2012t07CTD7d5_MS247", "SS2012t07UW36d5_MS247", "SS2012t07CTD8d5_MS247", "SS2012t07UW37d5_MS247", "SS2012t07UW38d5_MS247", "SS2012t07UW39d5_MS247", "SS2012t07UW41d5_MS247", "SS2012t07UW42d5_MS247", "SS2012t07UW43d5_MS247", "SS2012t07CTD9d5_MS247", "SS2012t07UW44d5_MS247", "SS2012t07CTD10d5_MS247", "SS2012t07UW45d5_MS247", "SS2012t07UW46d5_MS247", "SS2012t07UW47d5_MS247", "SS2012t07UW49d5_MS247", "SS2012t07UW50d5_MS247", "SS2012t07UW51d5_MS247", "SS2012t07CTD11d5_MS247", "SS2012t07CTD12d5_MS247", "SS2012t07CTD13d5_MS247", "SS2012t07UW55d5_MS247", "SS2012t07UW56d5_MS247", "SS2012t07CTD14d5_MS247", "SS2012t07UW58d5_MS247", "SS2012t07UW60d5_MS247", "SS2012t07CTD15d5_MS247", "SS2012t07CTD16d5_MS247", "SS2012t07UW61d5_MS247", "SS2012t07UW62d5_MS247", "SS2012t07UW63d5_MS247", "SS2012t07UW64d5_MS247", "SS2012t07CTD17d5_MS247")




short_names<-c("CTD1d5", "UW1d5", "UW2d5", "UW3d5", "UW4d5", "UW5d5", "UW14d5", "CTD2d5", "UW15d5", "UW16d5", "UW18d5", "UW20d5", "UW22d5", "UW23d5", "UW24d5", "UW25d5", "UW26d5", "CTD3d5", "UW27d5", "UW29d5", "CTD4d5", "CTD5d5", "UW31d5", "UW32d5", "UW33d5", "UW34d5", "CTD6d5", "UW35d5", "CTD7d5", "UW36d5", "CTD8d5", "UW37d5", "UW38d5", "UW39d5", "UW41d5", "UW42d5", "UW43d5", "CTD9d5", "UW44d5", "CTD10d5", "UW45d5", "UW46d5", "UW47d5", "UW49d5", "UW50d5", "UW51d5", "CTD11d5", "CTD12d5", "CTD13d5", "UW55d5", "UW56d5", "CTD14d5", "UW58d5", "UW60d5", "CTD15d5", "CTD16d5", "UW61d5", "UW62d5", "UW63d5", "UW64d5", "CTD17d5")
