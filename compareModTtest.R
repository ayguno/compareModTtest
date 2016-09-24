
#This function is designed to help us to interpret similarities between various treatments based on the proteins
#being co-effected by the individual treatments
# Multiple TMT10 experiments can be compared
#"ModT test table results" printed by the Shiny server will be provided into the same working directory, along with the class vector files specific for each experiment (both in csv format)
#The function will take a customly prepared csv file as input, which matches the data tables with class vector files
#The function will generate heatmaps of top 10 up/down regulated proteins for each experiment
#The function will also attempt to plot the entire set of experiments into a single scatterplot
#ModT_ClassV_match (chr) = csv table that matches ModT table file names (first column), with the respective Class vector files (second column) 
# keratins= T/F? (logical) False by default. If filtering for human Keratins, choose TRUE and also include "human_keratins_as_downloaded_from_HGNC_06242016.csv" in the working directory

# Functional as of 06/24/2016 (last update)


compareModTtest<-function(ModT_ClassV_match, keratins=F) { 
  
  require(gplots)
  require(RColorBrewer)
  
  
  mt<-read.csv(ModT_ClassV_match,header=TRUE) # Read the ModT_ClassV match table
  
  colnames(mt)<-c("ModT","ClassV","Experiment")
  
  Tindx<-NULL #Treatment index table that will be filled below with treatment names
  
  UP10<-c(1:10)
  DN10<-c(1:10)
  
  UnionUD<-NULL
  UnionUD<-as.list(UnionUD)
  
  UnionUD_gs_id<-NULL
  UnionUD_gs_id<-as.list(UnionUD_gs_id)
  
  geneU <-c(1:10)
  geneD <-c(1:10)
  pvalU <-c(1:10)
  pvalD <-c(1:10)
  AvExpU<-c(1:10)
  AvExpD<-c(1:10)
  
  n<-nrow(mt)
  
  for (i in 1:n) { #Loop over the experiments specified in the match table
    
    MT<-read.csv(file=as.character(mt$ModT[i]),header=TRUE) #Read the ModT results table for the experiment
    ClV<-read.csv(file=as.character(mt$ClassV[i]),header=TRUE) #Read the Class Vector table for the experiment
    
              if(keratins == T) {
                                KRT<-read.csv("human_keratins_as_downloaded_from_HGNC_06242016.csv",header =T)
                                KRT<-KRT$Approved
                                
                                w<-which(MT$geneSymbol %in%  KRT)
                                
                                MT<-MT[-w,]
                                }
    
                                
                                
    MT$geneSymbol<-as.character(MT$geneSymbol)
    
    mG<-which(MT$geneSymbol==c("")) #any rows that has missing gene symbol?
    
    if(length(mG)>0) {
      
      MT$geneSymbol[mG]<-as.character(MT$accession_number[mG]) #Replace any missing gene symbols with Uniprot IDs
      
    }
    
    
    #Exract individual treatments from this experiment
    
    a<-1
    n1<-nrow(ClV)/2
    
    for(i in 1:n1){ #Loop over the individual treatments within each experiment
      
      
      adjp<-which(colnames(MT)==paste("adj.P.Val",ClV[a,2],sep=".")) #Extract relevant column indexes for the specific treatment
      AExp<-which(colnames(MT)==paste("AveExpr",ClV[a,2],sep="."))
      p<-which(colnames(MT)==paste("P.Value",ClV[a,2],sep="."))
      id<-which(colnames(MT)=="accession_number") #finds the uniprot ID column
      
      gs_merge_id<-paste(as.character(MT[,"geneSymbol"]),as.character(MT[,id]),sep="____") #creates a new character vector with the same size as MT where gene symbols and uniprot ids are conctatenated to create a unique column
      
      E<-paste(ClV[a,2],sep="") #Create a new variable on-the-fly (named with the individual treatment) to store these information
      Tindx<-c(Tindx,E) #Collect the names of the treatment data frames for future access
      
      Edf<-data.frame(MT[,"geneSymbol"],MT[,adjp],MT[,AExp],MT[,p],gs_merge_id)#Extract the relevant columns for a specified treatment and coerce them into a data frame
      
      colnames(Edf)<-c(paste("geneSymbol",ClV[a,2],sep="."), 
                       paste("adj.P.Val",ClV[a,2],sep="."),
                       paste("AveExpr",ClV[a,2],sep="."),
                       paste("P.Value",ClV[a,2],sep="."),"gs_merge_id")
      
      #Report Top10 upregulated and downregulated gene symbols
      
      UP<-which(Edf[,3]>0) #First seperate up/down regulated for each treatment 
      DN<-which(Edf[,3]<0)
      
      EdfUP<-Edf[UP,]
      EdfDN<-Edf[DN,]
      
      EdfUPs<-EdfUP[order(EdfUP[,4]),] #then rank based on nominal p-value, get the top 10 genes
      EdfDNs<-EdfDN[order(EdfDN[,4]),]
      
      U<-head(EdfUPs,10)
      D<-head(EdfDNs,10)
      
      geneU<-data.frame(geneU,U[,1]) #Collecting the variables to be plotted into data frames 
      geneD<-data.frame(geneD,D[,1])
      pvalU<-data.frame(pvalU,U[,4])
      pvalD<-data.frame(pvalD,D[,4])
      AvExpU<-data.frame(AvExpU,U[,3])
      AvExpD<-data.frame(AvExpD,D[,3])
      
      
      UnionUD<-union(UnionUD,U[,1]) #Get the concatenated union of Top10 up and down-regulated proteins
      UnionUD<-union(UnionUD,D[,1]) #Finally UnionUD contains the genesymbols of the unions 
      
      UnionUD_gs_id<-union(UnionUD_gs_id,U[,5]) # Same a UnionUD but contains the union of merged gs_id identifiers
      UnionUD_gs_id<-union(UnionUD_gs_id,D[,5])
      
      UP10<-data.frame(UP10,U) #Now this information is concatenated into two tables which will be used below for reports 
      DN10<-data.frame(DN10,D)
      
      
      assign(E,Edf) #Assign the extracted data.frame into the on-the-fly variable
      
      
      
      a<-a+2 #Move to the next treatment within the experiment, until no treatment left
    }
    
    
    #Move to the next experiment in the match table, until all of them completed 
  }
  
  write.csv(file="Top10_Upregulated.csv",UP10)
  write.csv(file="Top10_Downregulated.csv",DN10)
  
  geneU<-geneU[,-1] #Dropping the first columns
  geneD<-geneD[,-1]
  pvalU<-pvalU[,-1]
  pvalD<-pvalD[,-1]
  AvExpU<-AvExpU[,-1]
  AvExpD<-AvExpD[,-1]
  
  colnames(geneU)<-Tindx
  colnames(geneD)<-Tindx  
  colnames(pvalU)<-Tindx
  colnames(pvalD)<-Tindx
  colnames(AvExpU)<-Tindx
  colnames(AvExpD)<-Tindx
  
  
  geneD<-as.matrix(geneD) #Converting the data.frames into matrices
  geneU<-as.matrix(geneU)
  
  pdf(file="Top10_up_down_regulated_proteins.pdf",width = 11) ##png will be better!!
  
  par(mfrow=c(1,2))
  
  colfuncUP <- colorRampPalette(c("white","red"))
  heatmap.2(as.matrix(AvExpU), Rowv = F, Colv = F, trace='n',  
            dendrogram = 'none', cellnote = geneU, notecol='black',
            margins=c(10,5), density.info = 'none',col=colfuncUP(100), 
            symbreaks=F,symkey = F,keysize = 1,
            main="Top10 upregulated proteins ranked by p-value",srtCol = 45,
            ylab = "nominal p-value rank (Smallest-to-Largest)",
            key.title=F, key.xlab = "Average Expression (log2)", 
  )
  
  
  colfuncDN <- colorRampPalette(c("blue","white"))
  heatmap.2(as.matrix(AvExpD), Rowv = F, Colv = F, trace='n', 
            dendrogram = 'none', cellnote = geneD, notecol='black',
            margins=c(10,5), density.info = 'none',col=colfuncDN(100),
            symbreaks=F,symkey = F, keysize = 1,
            main="Top10 downregulated proteins ranked by p-value",srtCol = 45,
            ylab = "nominal p-value rank (Smallest-to-Largest)",
            key.title=F, key.xlab = "Average Expression (log2)")
  
  dev.off()
  
  # See: http://www.inside-r.org/packages/cran/gplots/docs/heatmap.2 
  
  ##Next prepare the union matrix and the heatmap
  
  temp<-data.frame(unlist(UnionUD_gs_id)) 
  colnames(temp)<-c("UnionUD_gs_id")
  temp_UnionUD_AExp<-temp
  temp_UnionUD_p<-temp
  
  n<-nrow(mt)
  
  for (i in 1:n) { #Loop over the experiments specified in the match table
    
    MT<-read.csv(file=as.character(mt$ModT[i]),header=TRUE) #Read the ModT results table for the experiment
    ClV<-read.csv(file=as.character(mt$ClassV[i]),header=TRUE) #Read the Class Vector table for the experiment
    
              if(keratins == T) {
                KRT<-read.csv("human_keratins_as_downloaded_from_HGNC_06242016.csv",header =T)
                KRT<-KRT$Approved
                
                w<-which(MT$geneSymbol %in%  KRT)
                
                MT<-MT[-w,]
              }
    
    MT$geneSymbol<-as.character(MT$geneSymbol)
    
    mG<-which(MT$geneSymbol==c("")) #any rows that has missing gene symbol?
    
    if(length(mG)>0) {
      
      MT$geneSymbol[mG]<-as.character(MT$accession_number[mG]) #Replace any missing gene symbols with Uniprot IDs
      
    }
    
    
    #Exract individual treatments from this experiment
    
    a<-1
    n1<-nrow(ClV)/2
    
    for(i in 1:n1){ #Loop over the individual treatments within each experiment
      
      
      #Extract relevant column indexes for the specific treatment
      AExp<-which(colnames(MT)==paste("AveExpr",ClV[a,2],sep="."))
      p<-which(colnames(MT)==paste("P.Value",ClV[a,2],sep="."))
      gs<-which(colnames(MT)=="geneSymbol")
      id<-which(colnames(MT)=="accession_number") #finds the uniprot ID column
      
      gs_merge_id<-paste(as.character(MT[,"geneSymbol"]),as.character(MT[,id]),sep="____") #creates a new character vector with the same size as MT where gene symbols and uniprot ids are conctatenated to create a unique column
      
      MT<-data.frame(MT,gs_merge_id) #adds MT a new last column of "gs_merge_id"
      
      gs_AExp<-data.frame(MT[,"gs_merge_id"],MT[,AExp])
      colnames(gs_AExp)<-c("gs_merge_id",paste(ClV[a,2],sep="."))
      gs_p<-data.frame(MT[,"gs_merge_id"],MT[,p])
      colnames(gs_p)<-c("gs_merge_id",paste(ClV[a,2],sep="."))
      
      temp_UnionUD_AExp<-merge(temp_UnionUD_AExp,gs_AExp,all.x=T,by.x="UnionUD_gs_id",by.y="gs_merge_id") #Merge the Av.expression values of the UnionUD list genes into a data frame
      temp_UnionUD_p<-merge(temp_UnionUD_p,gs_p,all.x=T,by.x="UnionUD_gs_id",by.y="gs_merge_id") #Merge the nominal p values of the UnionUD list genes into a data frame
      
      a<-a+2 #Move to the next treatment within the experiment, until no treatment left
    }
    
    
    #Move to the next experiment in the match table, until all of them completed 
  }
  
  x<-as.matrix(temp_UnionUD_AExp[2:ncol(temp_UnionUD_AExp)],rownames.force = T)
  
  y<-sub('\\_.*', '', temp_UnionUD_AExp[,1])
  
  rownames(x)<-y
  
  
  m<-which(apply(is.na(x),1,sum)>0) #Finding the rows that contain at least one missing data (NA)
  
  x<-x[-m,] # remove the missing values to be able to cluster the genes
  
  colfuncUPDN <- colorRampPalette(c("blue","white","Red"))
  
  pdf("Commonly_identified_Top10_UP_DOWN_regulated_proteins_clustered.pdf",width =11, height = 11)
  
  
  
  heatmap.2(x, Rowv = T, Colv = T, trace='n', 
            dendrogram = 'both', notecol='black',
            margins=c(10,15), density.info = 'none',col=colfuncUPDN(1000), 
            symbreaks=T,symkey = T,keysize = 1,
            srtCol = 45,
            key.title=F, key.xlab = "Average Expression (log2)",na.color = "grey",symm = F,
            lwid=c(6,8),lhei = c(2,8)
  )
  
  
  dev.off()
  
  
  
  
  ##Next, preparing scatterplots where the significant hits are color marked
  
  
  pdf("Scatterplots.pdf",width=11,height=8)
  
  par(mfrow=c(1,3)) #pdf will be 3 pages long (rows)
  
  PVALUES<-c(0.05,0.1,0.25)    
  
  for(i in 1:3) { #Loop over the selected p-values
    
    Pv<-PVALUES[i]
    
    n<-nrow(mt)
    
    par(mfrow=c(n,4), mar=c(3,3,3,3)) # graphic parameters for each plot for a given page of pdf. 12 plots is a reasonable maximum
    
    for (i in 1:n) { #Loop over the experiments specified in the match table
      
      
      
      
      MT<-read.csv(file=as.character(mt$ModT[i]),header=TRUE) #Read the ModT results table for the experiment
      ClV<-read.csv(file=as.character(mt$ClassV[i]),header=TRUE) #Read the Class Vector table for the experiment
      
                  if(keratins == T) {
                    KRT<-read.csv("human_keratins_as_downloaded_from_HGNC_06242016.csv",header =T)
                    KRT<-KRT$Approved
                    
                    w<-which(MT$geneSymbol %in%  KRT)
                    
                    MT<-MT[-w,]
                  }
      
      
      a<-1
      n1<-nrow(ClV)/2
      
      for(i in 1:n1){ #plot the TMT10 ratio pairs of individual experiments in the ModT table
        
        
        
        adjp<-MT[,paste("adj.P.Val", ClV[a,2],sep=".")]
        
        x<-which(colnames(MT)==ClV[a,1])
        y<-which(colnames(MT)==ClV[a+1,1])
        
        plot(MT[,x],MT[,y],pch=20, col=ifelse(adjp<Pv,"red","black"), cex=1,xlim=c(-2,2),ylim=c(-2,2),main=ClV[a,2],cex.lab=0.8, xlab="Rep1",ylab="Rep2")
        abline(v=0,h=0,lwd=1,lty=2)
        
        a<-a+2
        
        text(paste("adj.p.val<",Pv,sep=""),x=1,y=-1.8,col="Red")
      } 
      
    }
    
  }
  
  dev.off()
}   