ApplyModel=F

Args <- commandArgs(trailingOnly=FALSE)

for (i in 1:length(Args)){
  if(Args[i]=="-f") ScriptPath<-Args[i+1]
  argSplit <- strsplit(Args[i], "=")
  if(argSplit[[1]][1]=="--file") ScriptPath <- argSplit[[1]][2]
}

for (arg in Args){
  argSplit <- strsplit(arg, "=")
  argSplit[[1]][1]
  argSplit[[1]][2]
  if(

ScriptPath<-dirname(dirname(ScriptPath))
source(file.path(ScriptPath,"LoadRequiredCode.r"))

aa=Sys.time()
library(parallel)
library(rgdal)
library(raster)

load(modelWorkspace, envir = environment())

if(ApplyModel!=F){
  csv = read.csv(ApplyModel, header=T, stringsAsFactors = F)
  out$dat$input$ParcTemplate = csv[2,1]
  csv = csv[2,4:length(csv)]
  out$dat$tif.ind = csv
  out$input$output.dir = file.path(dirname(ApplyModel))
  out$dat$bname = file.path(dirname(ApplyModel),out$input$script.name)
}

#pre-defined variables
out$input$make.binary.tif=FALSE
out$input$make.p.tif=FALSE
out$input$make.bin.tif=FALSE
out$input$MESS = TRUE
out$input$ResidMaps = FALSE
Model = out$input$script.name
make.p.tif=FALSE
make.bin.tif=FALSE
make.binary.tif=FALSE

#proc.tiff variables
model=out$mods$final.mod
vnames=names(out$dat$ma$train$dat)[-1]
tif.dir=out$dat$tif.dir$dname
filenames=out$dat$tif.ind
factor.levels=out$dat$factor.levels
make.binary.tif=make.binary.tif
thresh=out$mods$auc.output$thresh
make.p.tif=make.p.tif
outfile.p=paste(out$dat$bname,"_prob_map.tif",sep="")
outfile.bin=paste(out$dat$bname,"_bin_map.tif",sep="")
tsize=50.0
NAval=-3000
fnames=out$dat$tif.names
out=out
Model=out$input$script.name

# deconstruct proc.tiff
if(is.null(factor.levels)) factor.levels<-NA
MESS=MOD=out$input$MESS
if(is.null(thresh)) thresh<-.5
nvars<-length(vnames)
vnames.final.mod<-out$mods$vnames
nvars.final<-length(vnames.final.mod)
names(filenames)<-sub("_categorical","",names(filenames))
fullnames <- as.character(filenames[match(vnames,names(filenames))])
goodfiles <- file.access(fullnames)==0
if(!all(goodfiles)) stop(paste("ERROR: the following image files are missing:",paste(fullnames[!goodfiles],collapse=", ")))
if(nvars.final<1) MESS=FALSE
if(nvars.final==1) MOD=FALSE
gi <- GDALinfo(fullnames[1])
dims <- as.vector(gi)[1:2]
ps <- as.vector(gi)[6:7]
ll <- as.vector(gi)[4:5]
pref<-attr(gi,"projection")
RasterInfo=raster(fullnames[1])
RasterInfo@file@datanotation<-"FLT4S"
NAval<- -3.399999999999999961272e+38
if(!is.na(match("AREA_OR_POINT=Point",attr(gi,"mdata")))){                                                                          
  xx<-RasterInfo  #this shifts by a half pixel
  nrow(xx) <- nrow(xx) - 1
  ncol(xx) <- ncol(xx) - 1
  rs <- res(xx)
  xmin(RasterInfo) <- xmin(RasterInfo) - 0.5 * rs[1]
  xmax(RasterInfo) <- xmax(RasterInfo) - 0.5 * rs[1]
  ymin(RasterInfo) <- ymin(RasterInfo) + 0.5 * rs[2]
  ymax(RasterInfo) <- ymax(RasterInfo) + 0.5 * rs[2]
}
MB.per.row<-dims[2]*nvars*32/8/1000/1024
if(MESS) MB.per.row<-MB.per.row*8 #use more blocks for mess
nrows<-min(round(tsize/MB.per.row),dims[1])
bs<-c(nrows,dims[2])
chunksize<-bs[1]*bs[2]
tr<-blockSize(RasterInfo,chunksize=chunksize)
FactorInd<-which(!is.na(match(vnames,names(factor.levels))),arr.ind=TRUE)
if((nvars-length(FactorInd))==0) MESS<-MOD<-FALSE
if(tr$n<10 | getRversion()<2.14) multCore<-FALSE
tile.start<-seq(from=1,to=tr$n,by=ceiling(tr$n/(detectCores()-1))) 
outfile.p=file.path(out$input$output.dir,"ProbTiff","_prob_map.tif")


varlist = c('tr','dims','MESS','MOD','nvars','fullnames','nvars.final','vnames','NAval'
            ,'factor.levels','model','Model','pred.fct','make.binary.tif','make.p.tif','RasterInfo'
            ,'outfile.p','outfile.bin','thresh','vnames.final.mod', 'make.bin.tif','ScriptPath')
source(parRasterFile, local = F)

if(MESS) dir.create(file.path(out$input$output.dir,"MESSTiff"))
if(MOD) dir.create(file.path(out$input$output.dir,"ModTiff"))
no_cores = detectCores()-1
cl<-makeCluster(no_cores)
clusterExport(cl=cl, varlist = varlist, envir = environment())
clusterCall(cl=cl,fun =  function() { source(file.path(ScriptPath,"CalcMESS.r"), local = F) })
clusterCall(cl=cl,fun =  function() { source(file.path(ScriptPath,"pred.fct.r"), local = F) })
clusterEvalQ(cl=cl, expr = library(raster))
clusterEvalQ(cl=cl, expr = library(foreign))
clusterEvalQ(cl=cl, expr = library(gbm))
parLapply(cl,X = tile.start, fun=parRaster
          , dims=dims
          , factor.levels=factor.levels
          , fullnames=fullnames
          , maDir=out$input$ma.name
          , make.p.tif=make.p.tif
          , MESS=MESS
          , MOD=MOD
          , model=model
          , Model=Model
          , NAval=NAval
          , nToDo=ceiling(tr$n/(detectCores()-1))
          , nvars.final=nvars.final
          , outfile.bin=outfile.bin
          , outfile.p=outfile.p
          , pred.fct=pred.fct
          , RasterInfo=RasterInfo
          , thresh=thresh
          , tr=tr
          , train.dat=out$dat$ma$train$dat
          , vnames.final.mod=vnames.final.mod
          , vnames=vnames
          , template=out$dat$input$ParcTemplate
)
stopCluster(cl)

# Stitch rasters together in an efficient way.
library(gdalUtils)
library(foreign)

if(MESS){
  folder = file.path(out$input$output.dir,"MESSTiff")
  tifs = list.files(folder, full.names = T, pattern = '*.tif$')
  
  #Build VRT first to speed things up.
  vrt_out = paste0(folder,'/tmp.vrt',sep='')
  gdalbuildvrt(gdalfile = tifs, output.vrt = vrt_out, verbose = F, overwrite = T)
  
  #Then merge the files with gdal_translate
  translate_out = file.path(out$input$output.dir,paste0(out$input$script.name,"_mess_map.tif",sep=''))
  gdal_translate(src_dataset = vrt_out, dst_dataset = translate_out
                 , r = "near", co = c("COMPRESS=LZW","TILED=YES"))
  
  #Finally, clean up (delete) the folder since we dont' need those files anymore
  unlink(folder, recursive = T)
}

if(MOD){
  folder = file.path(out$input$output.dir,"ModTiff")
  tifs = list.files(folder, full.names = T, pattern = '*.tif$')
  
  #Build VRT first to speed things up.
  vrt_out = paste0(folder,'/tmp.vrt',sep='')
  gdalbuildvrt(gdalfile = tifs, output.vrt = vrt_out, verbose = F, overwrite = T)
  
  #Then merge the files with gdal_translate
  translate_out = file.path(out$input$output.dir,paste0(out$input$script.name,"_MoD_map.tif",sep=''))
  gdal_translate(src_dataset = vrt_out, dst_dataset = translate_out
                 , r = "near", co = c("COMPRESS=LZW","TILED=YES"))
  
  #Write out longest dbf file to capture MoD metadata
  dbfs = list.files(folder, full.names = T, pattern = '*.dbf$')
  for(i in dbfs){
    if(nrow(read.dbf(i))==max(nrow(read.dbf(dbfs)))){
      write.dbf(i, file.path(out$input$output.dir,paste0(out$input$script.name,'_MoD_map.tif.vat.dbf',sep='')))
      break
    }
  }
  #Finally, clean up (delete) the folder since we dont' need those files anymore
  unlink(folder, recursive = T)
}
bb=Sys.time()
bb-aa