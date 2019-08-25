################################################################
###############    Set the parameters here    ##################

CSV_Folder =   'G:/AnnalsCode/SpikeDetection_Rcode/DataFiles/CSV'           # The path of CSV (Electrophysiology_LFP) files
RoiXlsFolder = 'G:/AnnalsCode/SpikeDetection_Rcode/DataFiles/2P_Analysis'   # ROI .XLS 2-p calcium imaging data files (analysis file)
fs = 10000                   # original voltage sampling rate (Hz)
fsd = 50                     # voltage sampling rate (Hz) after downsampling

peakWindow = 0.1             # seconds
baselineWindow = 2           # seconds
minSpacingBefore = 0.5       # minimum spacing before (seconds)
minSpacingAfter = 0          # minimum spacing after (seconds)
maxFWHM=0.5                  # specify the maximal half width in seconds
trialsToAnalyze = 'all'      # specify trials to plot (e.g. 40:50), or set to "all"
FramePeriod= 0.032979438     # specify the frame acquisition period. 



################################################################
############# Load libraries, define basic functions ###########

loadLibrary = function(x) {
  if (!require(x, character.only=TRUE)) 
    install.packages(x)
  library(x, character.only=TRUE)
}

loadLibrary("caTools")          # caTools is needed for runquantile()
loadLibrary("readxl")

# define a simple string concatenation operator
"%+%" = function(x,y) paste(x,y,sep='')

# function to determine which points are local maxima
is.localmax = function(i,context, k, FUN=max, ...){
  context = c(rep(NA,k), context, rep(NA,k))           # pad with NAs
  i = i + k                                            # account for the NA padding
  contextMatrix = sapply(i, function(x) {context[(x-k):(x+k)]})
  z = context[i] == apply(contextMatrix, 2, FUN, ...)
  z[is.na(z)] = FALSE
  z
}

# function to find peaks (local maxima exceeding a threshold)
getPeaks = function(x,threshold,k) {
  i = which(x> threshold)
  if (length(i)>0)
    as.vector(na.omit(i[is.localmax(i,context=x, k=k,FUN=max)]))
  else
    NA
}

# for plotting purposes, to offset traces
offset = function(x, offset=TRUE, adj=1){
  if (is.numeric(offset)){
    if (length(offset) == 1)
      offset = rep(offset, dim(x)[2])
  } else {
    offset= cumsum(c(0,apply(x[, -ncol(x)],2, max, na.rm=T)))
  }
  t(t(x)+offset*adj)
}


################################################################
############# Load all LFP files, detect peaks #################

# get all the names of the lfp files (all the files in the folder that include ".csv")
lfpFiles = list.files(CSV_Folder, pattern='.csv', full.names=TRUE, recursive=TRUE)
lfpFiles = lfpFiles[grep('VoltageRecording', lfpFiles)]
if (length(grep('BACKUP', lfpFiles))>0)           # ignore backup files
  lfpFiles = lfpFiles[-grep('BACKUP', lfpFiles)]
if (is.numeric(trialsToAnalyze)){                 # option to select only specific trials to analyze
  lfpFiles = lfpFiles[trialsToAnalyze]
}

# read in the data, downsample it (takes up to a minute).
vData = unname(sapply(lfpFiles, function(f) {
  vData = read.csv(f)
  colMeans(matrix(vData[,'Input.7'],nrow=fs/fsd))
}))



### If you didnt change the folder but changed the parameters, start from here. ###


# estimate the baseline using a moving window filter that calculates the 90th percentile for each window. 
baseline = runquantile(vData, k=baselineWindow*fsd+1, probs=0.9)
vDataF = baseline-vData                  # subtract baseline, call this filtered signal (F)
threshold = 4 * mad(as.vector(vDataF))   # define threshold based on 4* median absolute deviation (robust stdev estimator)
peaks = apply(vDataF, 2, getPeaks, threshold=threshold, k=fsd*peakWindow) # find the peaks in each trial


FWHMs = sapply(1:ncol(vData), function(trial) {  # iterate over entire trial
  if (length(na.omit(peaks[[trial]]))==0)        # if there are no peaks, return NULL
    return(NULL)
  sapply(peaks[[trial]], function(p) {           # iterate over all the peaks in the trial
           min( which(vDataF[p:1,trial] < (0.5*vDataF[p,trial]))) + 
           min( which(vDataF[p:nrow(vData),trial] < (0.5*vDataF[p,trial]))) - 2  
         } 
  )
})


# select out peaks that meet the criterion of minimum spacing before and minimum spacing after
goodPeaks = sapply(1:ncol(vData), function(trial) {
  isGood = which(diff(c(peaks[[trial]],Inf)) > minSpacingAfter*fsd & 
                 diff(c(0,peaks[[trial]])) > minSpacingBefore*fsd & FWHMs[[trial]] < maxFWHM*fsd ) 
  peaks[[trial]][isGood]
})

twentyPercentPoints = sapply(1:ncol(vData), function(trial) {
  if (length(na.omit(goodPeaks[[trial]]))==0) return(NULL)
  sapply(goodPeaks[[trial]], function(p) max( which(vDataF[1:p,trial] < (0.2*vDataF[p,trial])) ))
})


#######################################################3
##################### Make LFP plots ##################

matplot(1:nrow(vDataF)/fsd, offset(vDataF), type='l', lty=1, col=1, xlab='time (s)', ylab='')
y =  offset(vDataF)[1,]
text(x=rep(0,length(y)), y=y, 1:length(y), cex=0.8)
for(j in 1:ncol(vDataF)){
  if(length(na.omit(peaks[[j]])) == 0) next
  points(peaks[[j]]/fsd, offset(vDataF)[peaks[[j]],j], col='red', pch=16,cex=0.9)
}
for(j in 1:ncol(vDataF)){
  if(length(na.omit(goodPeaks[[j]])) == 0) next
  points(goodPeaks[[j]]/fsd, offset(vDataF)[goodPeaks[[j]],j], col='blue', pch=16,cex=0.9)
}
for(j in 1:ncol(vDataF)){
  if(length(na.omit(twentyPercentPoints[[j]])) == 0) next
  points(twentyPercentPoints[[j]]/fsd, offset(vDataF)[twentyPercentPoints[[j]],j], col='green', pch=16,cex=0.9)
}
abline(v=1:80)

# Set up matrix containing 1) trial number 2) time of peak (seconds) and 3) frame of peak
peakTimes = do.call(rbind, lapply(1:length(twentyPercentPoints), function(i) { 
  if(length(na.omit(twentyPercentPoints[[i]]))>0) {
    cbind(i, twentyPercentPoints[[i]]/fsd) 
  } else {
    return(NULL)
  }
}))
colnames(peakTimes) = c('Trial Number', 'peakTimeSeconds')

Tstart = peakTimes[,'peakTimeSeconds'] * (1/FramePeriod) - 5
peakTimes = cbind(peakTimes, Tstart) 


##########################################################################
######## Now load in fluorescence data for each roi from each trial  #####

trialFluorFiles=list.files(RoiXlsFolder, pattern='ome.xls', full.names=TRUE, recursive=FALSE)

roi_df_list=list()

if (is.numeric(trialsToAnalyze)){
  trialFluorFiles = trialFluorFiles[trialsToAnalyze]
}
for(i in 1:length(trialFluorFiles)){
  mat=read_excel(trialFluorFiles[i],sheet = 2, col_types= ('numeric'))
  roi_df_list[[i]]=mat
}

# Calculate mean correlations between all ROIs for each peak
meanCorrelations = vector('numeric', nrow(peakTimes))
for (i in 1:nrow(peakTimes)){
  trialNumber = peakTimes[i,1]
  startFrame = peakTimes[i,3]
  if (startFrame+11 > nrow(roi_df_list[[trialNumber]]) | startFrame < 1){
    meanCorrelations[i] = NA
    next
  }
  mat = roi_df_list[[trialNumber]][startFrame+1:11, ] # want all the columns, but only 11 rows (timepoints)
  mat = as.matrix(mat)
  correlationMatrix = cor(mat)
  diag(correlationMatrix) = NA
  meanCorrelations[i] = mean(correlationMatrix, na.rm=TRUE)
  cat("event " %+% i %+%" in trial " %+% trialNumber %+%", mean correlation is " %+% meanCorrelations[i] %+% "\n")
  
}

mean(meanCorrelations)
sd(meanCorrelations)
