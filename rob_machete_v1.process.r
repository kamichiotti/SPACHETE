## simulated data anlysis
par(mfrow=c(2,4))
i=0

## eventually want one of these loops
#for (name in c("bladder","CML_test","CML_UConn","Engstrom","EwingKD","normal_breast","normal_fetal","ov_RNase_Qatar","ov_salzman","prostate","SEQC","tcga")){
#for (name in c("normal_fetal","ov_RNaseR_Qatar","Engstrom","Ewing","normal_breast" ,"CML_test","CML_UConn")){

# to test normal fetal-- look at matrix which is 'test'
## this is the thing to type in even though the script breaks-- test[!(junction %like% "no_fusion")]

spork_dir_stem = "/scratch/PI/horence/rob/spachete_outputs/"
#spork_dir_stem = "/scratch/PI/horence/rob/SPACHETE/"
#for (name in "CML_with_R2s"){ #Unfiltered
#for (name in "CML_new_gtf_test"){ #Filtered
#for (name in "test_new_pos_fix"){ #Fixed positions
#for (name in "normal_fetal_8_22"){
#for (name in "normal_fetal_9_8"){
for (name in "engstrom_9_9"){
#for (name in "CML_uconn_9_5_nup214"){
#for (name in "ewing_9_1"){
#for (name in "rerun_glms"){
#for (name in "normal_breast_9_6"){
for (spork in c(0,1)){
#for (name in "engstrom"){

dir=paste("/scratch/PI/horence/gillian/All_AppendedReports_Jun20/",name,"/",sep="")
# for spachete
if (spork==1){
dir=paste(spork_dir_stem,name,"/",sep="")
}

if (name == "tcga"){
dir=paste("/scratch/PI/horence/machoutput/")
dir="/scratch/PI/horence/sb/sbdata/";#jun27d0alpacxx1gcaagg/reports/AppendedReports/"
#dir2="/scratch/PI/horence/sb/data_dir/";#glmReports
}

print ("dir")
print(dir)
require(data.table)

for (tmyfile in list.files(dir)){
myfile=tmyfile

if ((spork==1) | (name =="tcga")){
myfiledir1=paste(dir,tmyfile,"/reports/AppendedReports/",sep="")
#myfiledir2=paste(dir,tmyfile,"/glmReports/",sep="")
print(myfile)
myfile=c(list.files(myfiledir1))
}
print(paste("myfile is ",myfile))

if (length(myfile)>0){
if ((myfile %like% "Appended" & myfile %like% "naive" )| (name=="tcga" & myfile %like% "circJuncProbs.txt_cdf") ){# & myfile %like% "SRR1594022"){

print (paste(myfile))

## diff directory structure for spork
if (!(spork==1)|(name == "tcga")){
if (!is.null(tryCatch(read.delim(paste(dir,myfile,sep=""),sep="\t"), error=function(e) NULL))){
print (myfile)
m=data.table(read.delim(paste(dir,myfile,sep=""),sep="\t"))
}
}
## diff directory structure for spork
if (name == "tcga"|spork==1){
if (!is.null(tryCatch(read.delim(paste(myfiledir1,myfile,sep=""),sep="\t"), error=function(e) NULL))){
print (myfile)
m=data.table(read.delim(paste(myfiledir1,myfile,sep=""),sep="\t"))
}
}

m[,sample:=paste("spork_",spork,"_",myfile,sep="")]
m[,sampleType:=name ]
g=data.frame(m)
names(m)[1]="junction"
if (dim(m)[2]>17){
names(m)[19]="badfj1is1"
names(m)[20]="badfj2is1"

if ((is.null(match("productPhat.y",names(m))))){
setnames(m,"p_predicted.x","productPhat.x")
setnames(m,"p_predicted.y","productPhat.y")
setnames(m,"p_value.x","junction_cdf.x")
setnames(m,"p_value.y","junction_cdf.y")


}
g=data.frame(m[ productPhat.y !="-" & productPhat.x!="-",])
print (head(g))
g$junction_cdf.y=as.numeric(as.vector(g$junction_cdf.y))
g$productPhat.y=as.numeric(as.vector(g$productPhat.y))
g$numReads.y=as.numeric(as.vector(g$numReads.y))
g$junction_cdf.x=as.numeric(as.vector(g$junction_cdf.x))
g$productPhat.x=as.numeric(as.vector(g$productPhat.x))
}

print ("finished ")
print (name)
if (i==0){allg=g}
if (i>0){allg=rbind(g,allg)}
i=i+1
print (head(allg))
}
}
}
}
}

print("finishedloop")
allg=data.table(allg)

allg[,junction:=gsub("([|])", ":", paste(junction))]
# example:     chr19:DAZAP1:1432689:+:chr22:SEPT5:19708072:+:fusion      0

allg[,chr1:= as.character(lapply(strsplit(paste(allg$junction), split=":"), "[", 1))]
allg[,gene1:= as.character(lapply(strsplit(paste(allg$junction), split=":"), "[", 2))]
allg[,pos1:= as.numeric(as.character(lapply(strsplit(paste(allg$junction), split=":"), "[", 3)))]
allg[,strand1:= as.character(lapply(strsplit(paste(allg$junction), split=":"), "[", 4))]

allg[,chr2:= as.character(lapply(strsplit(paste(allg$junction), split=":"), "[", 5))]
allg[,gene2:= as.character(lapply(strsplit(paste(allg$junction), split=":"), "[", 6))]
allg[,pos2:= as.numeric(as.character(lapply(strsplit(paste(allg$junction), split=":"), "[", 7)))]
allg[,strand2:= as.character(lapply(strsplit(paste(allg$junction), split=":"), "[", 8))]

allg[,type:= as.character(lapply(strsplit(paste(allg$junction), split=":"), "[", 9))]

## star used:
#g= g[SpliceType %like% "ONLY" & (chr1!=chr2 | abs(pos1-pos2)>1000000 | strand1!=strand2 )]
allg= allg[ (chr1!=chr2 | abs(pos1-pos2)>1000000 | strand1!=strand2 )]

p_pred_thresh=0
p_val_thresh=.7

##MOD
#pass.thresh= allg[  junction_cdf.y > p_val_thresh][order(-as.numeric(as.vector(numReads.y))) ]
#if (sum(names(allg) %like% "bad")>1){
#pass.thresh= allg[ (badfj1is1+badfj2is1)==0 & junction_cdf.y > p_val_thresh][order(-as.numeric(as.vector(numReads.y))) ]
# was .4
#}
#pass.thresh$numReads.y=as.numeric(as.vector(pass.thresh$numReads.y))
#pass.thresh$numReads.x=as.numeric(as.vector(pass.thresh$numReads.x))

## informative plots and calling thresholds
i=0
#pdf (file="machete_plots.pdf")
par(mfrow=c(2,2))

for (my.sample in unique (allg$sample)){

#suff. many counts: >3
train.good=allg [ paste(sample)==my.sample & (badfj1is1+badfj2is1)==0]

# these should be junctions that are bad whereas fj2is bad are more likley to just be linear artifacts

train.bad=allg[(sample %like% "spork_0") &(!(junction %like% "no_fusion")) & paste(sample)==my.sample & badfj1is1>0 ]

if (  dim(train.bad) [1]>0 ){
plot(train.bad[sample==my.sample]$junction_cdf.y, train.bad[sample==my.sample]$numReads.y, main=paste(unique(train.bad[sample==my.sample,sampleType]),"badFJs"))
}
if (  dim(train.bad) [1] == 0 ){
plot(c(0,0), main=paste("NODATA", unique(allg[sample==my.sample,sampleType]),"badFJs"))
}

if (  dim(train.good) [1] == 0 ){
plot(c(0,0), main=paste("NODATA", unique(allg[sample==my.sample,sampleType]),"goodFJs"))
}
if (  dim(train.good) [1] > 0 ){

plot(train.good[sample==my.sample]$junction_cdf.y, train.good[sample==my.sample]$numReads.y, main=paste(unique(train.good[sample==my.sample,sampleType]),"goodFJs", my.sample))
}
INC=100

train.good[, discrete_p:=round(INC*junction_cdf.y)/INC]
for (th in c(0:INC)/INC){
## assign a threshold
tot.bad=length(train.bad$junction_cdf.y)
tot.good=length(train.good$junction_cdf.y)
prop=tot.bad/tot.good
print ("adding proportion of good and bad-- CHANGE")
current.p=sum(train.bad$junction_cdf.y>th)/tot.bad

print(paste("theshold is",th,"current ",current.p))			      
train.good[discrete_p==th,emp_p:=current.p]
print(train.good[discrete_p==th,])
print (paste("current thresh and fdr",th,current.p))
}

if (i>0 & dim(train.good)[1]>0){all.passed=rbind(all.passed,train.good)}
if (i==0){
all.passed=train.good
}
i=i+1
}

print ("COMPLETED")
#

## falsely called over all clled:
## find the good threshold, then report posterior and fdr
all.passed[emp_p=="NaN", emp_p:=0]

# misnomer!

all.passed[,sinfo:=paste(sample,sampleType)]
all.passed[,fusionInfo:=junction]

min.reads=1
min.reads.high.cdf=10
p.thresh=.1
lower.jun.cdf=.2
#test= all.passed[(emp_p<p.thresh &  junction_cdf.y >lower.jun.cdf & numReads.y>min.reads) | (junction_cdf.y >(1-p.thresh) & numReads.y==1 & productPhat.y>.5)]

test= all.passed[(emp_p<p.thresh &  junction_cdf.y >lower.jun.cdf & numReads.y>min.reads) | (junction_cdf_lower.y >.5 & numReads.y==1 & productPhat.y>.5)]

sam=unique(test$sample)
test[,sampleid:=match(sample,sam)]

use.in=data.frame(test[,list(junction,sampleType,sample, numReads.y)])
use.in=cbind(use.in,paste(use.in$sample,use.in$sampleType))
colnames(use.in)[5]="sampleid"
use.in=data.table(use.in)
dr=data.frame(reshape(use.in[,list(junction,sampleid,numReads.y)], timevar="sampleid",idvar="junction",direction="wide"))
dr[is.na(dr)]=0
mxplot=as.matrix(dr[,2:dim(dr)[2]])
rownames(mxplot)=dr[,1]
mxplot[(mxplot>0)]=1
# heatmap(mxplot[order(rownames(mxplot)),], Rowv=NA, scale="none")

## CML replicates
#or all.passed
machete.info=test

test[,sinfo:=paste(sample,sampleType)]
test[,fusionInfo:=junction]

print ("undo loop when generating figures")
if (0==1){
for (testtype in c("Ewing","CML")){ #prostate 

cml.info=test[   sampleType %like% testtype,list(sinfo,fusionInfo, numReads.y)]
cml.wide=reshape(data.frame(cml.info),timevar="sinfo", idvar="fusionInfo", direction="wide")
rownames(cml.wide)=cml.wide[,1]
cml.wide[is.na(cml.wide)]=0
print(cor(cml.wide[,2:dim(cml.wide)[2]], method="spearman"))

#X11()
cml.wide[cml.wide>0]=1
hist(cor(cml.wide[,2:dim(cml.wide)[2]], method="spearman"), main=paste(testtype," hist Machete"))


#X11()

cml.wide= as.matrix(cml.wide[,2:dim(cml.wide)[2]])
cml.wide=cml.wide[,order(colnames(cml.wide))]
cml.wide=cml.wide[order(rownames(cml.wide)),]
heatmap(cml.wide, Rowv=NA, Colv=NA, scale="none", main=paste("machete heatmap min reads", min.reads, testtype))
}


## used to test whether the fusions are 'real'-- test is what we call 'good'

machete.allg=test

## used to test if there are artifacts -- COMMENT if we want 'real' fusions

prostate_mach= allg[sampleType %like% "prostate"  ] 
# 'any detected fusion'
#prostate_mach= allg[sampleType %like%"prostate" & badfj1is1==1 ]

# 'passed our filter"
#prostate_mach=machete.allg[sampleType %like% "prostate"]
# train.bad
prostate_mach[,fusion:=paste(gene1,gene2,sep="-")]
val=data.table(read.delim("prostate_validated_fusions.txt"))
val[,fusion:=paste(gene1,gene2,sep="-")]
# will count the number of validated fusions 
mergeprostate_mach=merge(prostate_mach,val,by="fusion")
print(unique(mergeprostate_mach[, list(fusion, validated)]))
}
#dev.off()
write.table(file="August_2016_V2_supp_machet_processed.tab",test,sep="\t",quote=F)
