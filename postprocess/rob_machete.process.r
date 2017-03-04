## simulated data anlysis
par(mfrow=c(2,4))
rm(list=ls())
i=0

require(data.table)

## eventually want one of these loops
#for (name in c("bladder","CML_test","CML_UConn","Engstrom","Ewing","normal_breast","normal_fetal","ov_RNase_Qatar","ov_salzman","prostate","SEQC","tcga")){

for (name in c("normal_fetal","ov_RNaseR_Qatar","Engstrom","Ewing","normal_breast" ,"CML_test","CML_UConn")){
  # to test normal fetal-- look at matrix which is 'test'
  ## this is the thing to type in even though the script breaks-- test[!(junction %like% "no_fusion")]
  
  #for (name in "Ewing"){ #Unfiltered
  #for (name in "CML_UConn"){
  #for (name in c("Engstrom","normal_breast","normal_fetal","Ewing")){
  #for (name in "normal_breast"){
  #for (name in "normal_fetal"){
  #for (name in "Engstrom"){
  #for (name in "Ewing"){
  #for (name in "ascor"){
  #for (name in "salzman_ov_2016"){
  for (spork in c(0,1)){# should be 0,1
    
    #RB This just redirects file names to the SPORK equivalents
    if (spork==1){
      if (name %like% "salzman"){name="runs_not_for_machete_paper/peter_out_9_19"}
      if (name %like% "CML_U"){ name ="CML_uconn_9_16"}
      if (name %like% "CML_test"){ name ="cml_test_9_21"}
      if (name %like% "Ewing"){ name ="ewing_9_14"}
      if (name %like% "fetal"){ name ="normal_fetal_9_16"}
      if (name %like% "breast"){ name ="normal_breast_9_21"}
      if (name %like% "Engstrom"){name="engstrom_9_9"}
      if (name %like% "ov_RNaseR_Qatar"){name="ovcar3_9_15"}
      
    }
    #CML_filtered"}
    #for (name in "engstrom"){
    #RB This is the newest appended reports for MACHETE from Gillian (all of them)
    dir=paste("/scratch/PI/horence/gillian/All_AppendedReports_Sep18/",name,"/",sep="")
    #dir=paste("/scratch/PI/horence/gillian/All_AppendedReports_Jun20/",name,"/",sep="")
    #dir=paste("/scratch/PI/horence/gillian/AllGLM_May1/",name,"/",sep="")
    # for spachete
    #RB directory path is running SPACHETE processing
    if (spork==1){
      dir=paste("/scratch/PI/horence/rob/spachete_outputs/",name,"/",sep="")
    }
    #RB For Julia's own testing of SPORK
    if (spork==2){
      dir=paste("/scratch/PI/horence/julia/",name,"/",sep="")
    }
    
    #RB This is for 7B, not important for MACHETE paper 1.0
    if (name == "tcga"){
      dir=paste("/scratch/PI/horence/machoutput/")
      dir="/scratch/PI/horence/sb/sbdata/";#jun27d0alpacxx1gcaagg/reports/AppendedReports/"
      #dir2="/scratch/PI/horence/sb/data_dir/";#glmReports
    }
    
    print ("dir")
    print(dir)
    require(data.table)
    #RB Not important for MACHETE paper 1.0
    if (name=="ov_salzman"){
      dir="/scratch/PI/horence/gillian/AllGLM_May1/ov_salzman/AppendedReports/"
    }
    
    #RB This is the AppendedReports read-in loop
    for (tmyfile in list.files(dir)){
      myfile=tmyfile
      
      if ((spork==1) | (name =="tcga")){
        myfiledir1=paste(dir,tmyfile,"/reports/AppendedReports/",sep="")
        
        myfile=c(list.files(myfiledir1))
      }
      if ((spork==2)){
        myfiledir1=paste(dir,tmyfile,"/",sep="")
        myfile=c(list.files(myfiledir1))
      }
      print(paste("myfile is |",myfile,"|",sep=""))
      
      #RB These two if-checks make sure that the file exists and has Appended and
      #RB naive in the name so that the circ and badfj reports aren't used
      if (length(myfile)>0){
        if ((myfile %like% "Appended" & myfile %like% "naive" )| (name=="tcga" & myfile %like% "circJuncProbs.txt_cdf") ){# & myfile %like% "SRR1594022"){
          
          print (paste(myfile))
          #RB Putting my own file output here to check that everything is run
          write(paste(spork,myfile,sep="  "),file="names_file.txt",append=TRUE)
          
          
          ## diff directory structure for spork
          #RB Read in MACHETE styled outputs
          if (!(spork==1)|(name == "tcga")){
            if (!is.null(tryCatch(read.delim(paste(dir,myfile,sep=""),sep="\t"), error=function(e) NULL))){
              print (myfile)
              m=data.table(read.delim(paste(dir,myfile,sep=""),sep="\t"))
            }
          }
          ## diff directory structure for spork
          #RB Read in SPACHETE styled outputs
          if (name == "tcga"|spork==1){
            if (!is.null(tryCatch(read.delim(paste(myfiledir1,myfile,sep=""),sep="\t"), error=function(e) NULL))){
              print (myfile)
              m=data.table(read.delim(paste(myfiledir1,myfile,sep=""),sep="\t"))
            }
          }
          if (spork==2){
            if (!is.null(tryCatch(read.delim(paste(myfiledir1,myfile,sep=""),sep="\t"), error=function(e) NULL))){
              print (myfile)
              m=data.table(read.delim(paste(myfiledir1,myfile,sep=""),sep="\t"))
            }
          }
          
          #RB Assigning spork number to the front of the sample field
          m[,sample:=paste("spork_",spork,"_",myfile,sep="")]
          m[,sampleType:=name ]
          #RB Create a data frame to store this AppendedReport's junctions
          g=data.frame(m)
          names(m)[1]="junction"
          #RB Assigning a badjf1 or badfj2
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
          #RB This is where the current data frame, g, gets added vertically
          #RB into a single large dataframe called allg
          if (i==0){allg=g}
          if (i>0){allg=rbind(g,allg)}
          i=i+1
          print (head(allg))
        }
      }
    }
  }
}

#RB End of Read-In loop
print("finishedloop")
allg=data.table(allg)

#RB All of this is parsing out the junction header as fields for the allg object
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
#RB set allg to be filled with only junctions as defined by distances
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

#RB looping through all samples only once (removed duplicates w/ Unique)
for (my.sample in unique (allg$sample)){
  
  #RB now train.good is full of jcts that are not badfj1 or badjf2
  #suff. many counts: >3
  train.good=allg [ paste(sample)==my.sample & (badfj1is1+badfj2is1)==0]
  
  # these should be junctions that are bad whereas fj2is bad are more likley to just be linear artifacts
  
  if (spork==1){
    # cryptic exons may be detected by spork and also badfj1
    #RB Should the (sample %like% spork) be (sample %like% "spork")
    #RB This is just looking for badfj1's only?
    #RB This line is instantly overwritten by the next line though right?
    train.bad=allg[(sample %like% spork) &(!(junction %like% "no_fusion")) & paste(sample)==my.sample & badfj1is1>0 & badfj2is1==0 ]
    
    #RB Will this line ever work since we only get in this if if spork==1?
    train.bad=allg[(sample %like% "spork_0") & paste(sample)==my.sample & badfj1is1>0 ]
  }
  if (spork==0){
    train.bad=allg[(sample %like% spork) & paste(sample)==my.sample & badfj1is1>0 ]
  }
  
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
  
  #RB The threshold between good and bad jct_cdf is averaged to get a p-value
  #RB The train-goods are already going to be added to the all.passed, they are just
  #RB getting their p-values updated
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
#RB End of assigning p-values to all.passed
print ("COMPLETED")
#

## falsely called over all clled:
## find the good threshold, then report posterior and fdr
all.passed[emp_p=="NaN", emp_p:=0]

# misnomer!

all.passed[,sinfo:=paste(sample,sampleType)]
all.passed[,fusionInfo:=junction]

#RB These are some of the hard thresholds being created
min.reads=1
#min.reads.high.cdf=10
p.thresh=.1

lower.jun.cdf=.4 # 
## sqrt motivation by scaling of standard error for a junction cdf measuement as 1/sqrt(n)

#!!!
#RB This is filtering out the all.passed that pass the thresholds, main filtering step
#RB Has a separate test for is the numreads == 1 or >= 2 (more stringent in the first case)
#!!!
test= all.passed[(emp_p<p.thresh &  junction_cdf.y >.2 & numReads.y>min.reads) | ( numReads.y==1 & junction_cdf_lower.y >.5 & productPhat.y>.5)]

#RB This is a list of all the unique samples
sam=unique(test$sample)
test[,sampleid:=match(sample,sam)]

test[,sinfo:=paste(sample,sampleType)]
test[,fusionInfo:=junction]

#Finally get rid of all the no_fusion annotations (first time looking at this annotation)
stest=test[!(fusionInfo %like% "no_fusion")]

write.table(file="September_2016_V2_supp_machete_processed.tab",stest,sep="\t",quote=F)
