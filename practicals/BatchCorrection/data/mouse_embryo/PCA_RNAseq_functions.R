# Created by Asa Björklund 1510112
# functions for plotting and running PCA on RNAseq data


run.pca<-function(data,seln=0,log.vals=TRUE,samples.col=TRUE,center=TRUE){
	# data should be a matrix with expression values (counts/rpkms/fpkms)
	# with genes as rows, samples as columns

	# options:
	# seln - select top varied genes to use in analysis (default is 0 == use all genes)
	# log.vals=TRUE, if values are already in log-format
	# samples.col=TRUE, set to false if data has samples as columns
	# center=TRUE, if pca should be run with/without centering.

	if (!samples.col){ data<-t(data) }

	# remove zero read genes
	z<-which(rowSums(data)==0)
	if (length(z>0)){
	    data<-data[-z,]
	}

	if (log.vals) { data <- log2(data+1) }
	
	# select top varied genes
	if (seln>0){
	    sD<-apply(data,1,sd)
	    o<-order(sD,decreasing=T)
	    data<-data[o[1:seln],]
	}

	# run PCA
        myPca <- prcomp(t(data),center=center,scale.=FALSE)

	# add in contribution to variance in PC-names
        vars<- myPca$sdev^2
        vars<- vars/sum(vars)
        pcnames<-mat.or.vec(1,length(vars))
        for (i in 1:length(vars)) {
            pcnames[i]<-sprintf("PC%d\n%.5f",i,vars[i])
        }

	myPca$pc.names<-as.vector(pcnames)
        return(myPca)
}
		

 
pca.plot <- function (data,col="black",pch='*', main="PCA",selpc=c(1,2),cex=0.7, seln=0,log.vals=TRUE,center=TRUE,samples.col=TRUE){
	 # data is either a prcomp object or a data matrix
	 # all variables for run.pca function can be set here as well
	 # also defines color (col), symbol (pch), title (main), which PCs to use (selpc), and text scaling (cex)

	 if (class(data) != "prcomp"){
	    data<-run.pca(data,seln=seln,log.vals=log.vals,samples.col=samples.col,center=center)
	 }
	 
	 tmpPca <- as.data.frame(data$x[,selpc])
	 colnames(tmpPca)<-data$pc.names[selpc]
    	 plot(tmpPca, col=col,pch=pch,main=main,cex=cex)
         return(invisible(data))
}


pca.loadings.plot <- function(data,main="pca loadings",seln=0,log.vals=TRUE,samples.col=TRUE,nPlot=20,selPC=1:5,center=TRUE,horizontal=TRUE){
         # data is either a prcomp object or a data matrix

         if (class(data) != "prcomp"){
            data<-run.pca(data,seln=seln,log.vals=log.vals,samples.col=samples.col,center=center)
         }


         aload <- abs(data$rotation[,selPC])
         contr<- sweep(aload, 2, colSums(aload), "/")
         par(mfrow=c(length(selPC),2),mar=c(1,10,2,1),oma=c(1,5,1,1),cex=0.5)
         for (i in selPC) {
             top<-order(data$rotation[,i],decreasing=T)[1:nPlot]
             bottom<-order(data$rotation[,i],decreasing=F)[1:nPlot]
             barplot(contr[top,i],main=sprintf("genes on pos axis PC%d",i),ylab="% contr",las=2,horiz=horizontal)
             barplot(contr[bottom,i],main=sprintf("genes on neg axis PC%d",i),ylab="% contr",las=2,horiz=horizontal)
         }
         invisible(data)
}


pca.contribution.plot <- function(data,nPlot=20){
	 # data is either a prcomp object or a data matrix

	 if (class(data) != "prcomp"){
	    data<-run.pca(data,seln=seln,log.vals=log.vals,samples.col=samples.col,center=center)
	 }
         par(mfrow=c(1,1))
         vars<- data$sdev^2
         vars<- vars/sum(vars)
         barplot(vars[1:nPlot],names.arg=1:nPlot)	 
	 return(invisible(data))
}    


plot.pca.biplot <- function(PC,selpc=1:2,add=FALSE,mode="dist",nPlot=50,selection=NULL,scale.factor=1,cex=0.5){
    # select genes by defining mode:
    # dist -  maximum distance from origo
    # both - top genes on both PCs
    # all - plot all genes
    # selection - plot genes as defined in vector "selection"

    x<-selpc[1]
    y<-selpc[2]
    data<-PC$x[,c(x,y)]
    rot<-PC$rotation[,c(x,y)]
    #multiplication factor
    mult <- min(
        (max(data[,2]) - min(data[,2])/(max(rot[,2])-min(rot[,2]))),
        (max(data[,1]) - min(data[,1])/(max(rot[,1])-min(rot[,1])))
        )
    v1<-rot[,1]*mult*scale.factor
    v2<-rot[,2]*mult*scale.factor

    sel.genes<-c()
    nG<-length(v1)
    if (mode == "dist"){
    # take the genes with top dist? SELL not included, will be on 53rd place
      dist<-v1^2+v2^2
      o<-order(dist,decreasing=T)
      sel.genes<-o[1:nPlot]
    }else if (mode == "both") {
      #select top genes for each pc instead
      nSel<-round(nPlot/4)
      top1<-c(head(order(v1),nSel), tail(order(v1),nSel))
      top2<-c(head(order(v2),nSel), tail(order(v2),nSel))
      sel.genes<-union(top1,top2)
    }else if (mode == "all") {
      sel.genes<-1:nG
    }else if (mode == "selection") {
      sel.genes<-match(selection,names(v1))
    }else{
      stop(sprintf("Wrong definition of mode: %s",mode))
    }

    if (!add) {
       # make a white plot
        plot(data,pch=16,col="white")
    }
    text(v1[sel.genes],v2[sel.genes],names(v1[sel.genes]),cex=cex,col="black")
}