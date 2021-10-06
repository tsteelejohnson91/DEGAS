#***************************************************************
#***************************************************************
# Single Cell Transfer Learning Toolbox
# Travis S Johnson and Zhi Huang
# Compatible with OSX and Linux

#***************************************************************
#***************************************************************
# Initializing the package

initDEGAS <- function(){
  DEGAS.pyloc <<- "python3"
  DEGAS.toolsPath <<- paste0(.libPaths()[1],"/DEGAS/tools/")
  DEGAS.train_steps <<- 2000
  DEGAS.scbatch_sz <<- 200
  DEGAS.patbatch_sz <<- 50
  DEGAS.hidden_feats <<- 50
  DEGAS.do_prc <<- 0.5
  DEGAS.lambda1 <<- 3.0
  DEGAS.lambda2 <<- 3.0
  DEGAS.lambda3 <<- 3.0
  DEGAS.seed <<- "NULL"
}

#***************************************************************
#***************************************************************
# General utility functions

# Return filename extension
getExtension <- function(file){
  ex <- strsplit(basename(file), split="\\.")[[1]]
  return(ex[length(ex)])
}

# Return operating system
checkOS <- function(){
  return(Sys.info()['sysname'])
}

# Return the python version information
checkForPy <- function(){
  return(system(paste0(DEGAS.pyloc," -V")))
}

# Return if tensorflow loading information
checkForTF <- function(){
  return(system(paste0(DEGAS.pyloc," -c 'import tensorflow'")))
}

# Manually reset the python path
setPython <- function(path2python){
  DEGAS.pyloc <<- path2python
  #Sys.setenv(PATH = paste(c(path2python,Sys.getenv("PATH")),collapse = .Platform$path.sep))
}

# Manually reset the number of training steps
set_training_steps <- function(inp){
  DEGAS.train_steps <<- inp
}

# Manually reset the single cell batch size
set_single_cell_batch_size <- function(inp){
  DEGAS.scbatch_sz <<- inp
}

# Manually reset the single patient batch size
set_patient_batch_size <- function(inp){
  DEGAS.patbatch_sz <<- inp
}

# Manually reset the number of hidden features
set_hidden_feature_number <- function(inp){
  DEGAS.hidden_feats <<- inp
}

# Manually reset the dropout keep percentage (the percentage of nodes to keep)
set_dropout_keep_fraction <- function(inp){
  DEGAS.do_prc <<- inp
}

# Manually reset the L2 regularization term (lambda 3)
set_l2_regularization_term <- function(inp){
  DEGAS.lambda1 <<- inp
}

# Manually reset the patient loss term (lambda 1)
set_patient_loss_term <- function(inp){
  DEGAS.lambda2 <<- inp
}

# Manually reset the MMD loss term (lambda 2)
set_MMD_loss_term <- function(inp){
  DEGAS.lambda3 <<- inp
}

# Manually reset the seed
set_seed_term <- function(inp){
  if(is.null(inp)){
    DEGAS.seed <<- "NULL"
   }else if (is.numeric(inp) & length(inp)==1){
     DEGAS.seed <<- floor(inp)
   }else{
    message("ERROR: Input proper seed (NULL or integer)")
    message("Non-integer values will equal floor(value)")
   }
 }

#***************************************************************
#***************************************************************
# ccModel class and related functions

# ccModel class
setClass("ccModel",slots=list(Bias="list",Theta="list",Activation="list",
                              Depth="numeric",Model_type="character",Architecture="character"))

# zscore normalization
normFunc <- function(x){return((x-mean(x, na.rm = T))/(sd(x, na.rm = T)+1e-3))}

# scaling from 0-1
scaleFunc <- function(x){return((x- min(x, na.rm = T)) /(max(x, na.rm = T)-min(x, na.rm = T)+1e-3))}

# Preprocess count data
preprocessCounts <-function(X){
  return(t(apply(t(apply(as.matrix(t(log2(X+1))),1,normFunc)),1,scaleFunc)))
}

# center to 0
centerFunc <- function(x){return(x-mean(x,na.rm=T))}

# Activation functions and utilities

# Sigmoid activation function
sigmoid <- function(x) {
  1 / (1 + exp(-x))
}

# Log sum exp transformation (for softmax)
logsumexp <- function (x) {
  y = max(x)
  y + log(sum(exp(x - y)))
}

# Softmax activation function
softmax <- function (X) {
  return(t(apply(X,1,function(x) exp(x - logsumexp(x)))))
}

# Onehot labels to list of label names
fromOneHot <- function(labMat){
  return(apply(labMat,1,function(x) colnames(labMat)[which(x==1)]))
}

# List of label names to a onehot matrix with labels as column names
toOneHot <- function(labels){
  labs = unique(labels)
  out = matrix(0,length(labels),length(labs))
  colnames(out) = labs
  row.names(out) = row.names(labels)
  for(i in 1:length(labels)){
    out[i,labels[i]] = 1
  }
  return(out)
}

# Convert matrix of output weights to max value for each row (row max = 1 and not row max = 0)
probtoOneHot <- function(probMat){
  idx = apply(probMat,1,function(x) which(x==max(x)))
  probMat = probMat*0
  for (i in 1:length(idx)){
    probMat[i,idx[i]] = 1
  }
  return(probMat)
}

# IO utilities for training
# Write input to tmp files for tensorflow python algorithm from R
writeInputFiles <- function(scExp,scLab,patExp,patLab,tmpDir){
  write.table(scExp,file=paste0(tmpDir, '/scExp.csv'), row.names=FALSE, sep=',')
  if(!is.null(scLab)){write.table(scLab,file=paste0(tmpDir, '/scLab.csv'), row.names=FALSE, sep=',')}
  write.table(patExp,file=paste0(tmpDir, '/patExp.csv'), row.names=FALSE, sep=',')
  if(!is.null(patLab)){write.table(patLab,file=paste0(tmpDir, '/patLab.csv'), row.names=FALSE, sep=',')}
}

# Read output files from tensorflow python algorithm into R
readOutputFiles <- function(tmpDir,Model_type,architecture){
  activations = read.table(file=paste0(tmpDir, 'Activations.csv'), row.names=NULL, header=FALSE, sep=',',stringsAsFactors=FALSE)
  activations = activations[[1]]
  depth = length(activations)
  cnt=0
  for (activation in activations){
    cnt=cnt+1
    eval(parse(text=paste0("Bias",as.character(cnt)," = as.matrix(read.table(file=paste0(tmpDir, 'Bias",as.character(cnt),".csv'), row.names=NULL, header=FALSE, sep=','))")))
    eval(parse(text=paste0("Theta",as.character(cnt)," = as.matrix(read.table(file=paste0(tmpDir, 'Theta",as.character(cnt),".csv'), row.names=NULL, header=FALSE, sep=','))")))
  }
  Biases = list()
  Thetas = list()
  Activations = list()
  cnt=0
  for (activation in activations){
    cnt=cnt+1
    eval(parse(text=paste0("Biases[[",as.character(cnt),"]]=","Bias",as.character(cnt))))
    eval(parse(text=paste0("Thetas[[",as.character(cnt),"]]=","Theta",as.character(cnt))))
    eval(parse(text=paste0("Activations[[",as.character(cnt),"]]=activation")))
  }
  return(new('ccModel',Bias=Biases,Theta=Thetas,Activation=Activations,Depth=depth,Model_type=Model_type,Architecture=architecture))
}

# Make python executable for standard (feedforward) implementation
makeExec <- function(tmpDir,FFdepth,model_type){
  if (model_type != 'ClassClass' && model_type != 'ClassCox' && model_type != 'ClassBlank' && model_type != 'BlankClass' && model_type!='BlankCox'){
    stop("Please specify either 'BlankClass', 'ClassBlank', 'BlankCox', ClassClass' or 'ClassCox' for the model_type")
  }
  system(paste0('cp ',DEGAS.toolsPath,model_type,'MTL_p1.py ',tmpDir))
  system(paste0('cp ',DEGAS.toolsPath,model_type,'MTL_p3.py ',tmpDir))
  outlines = c()
  if (FFdepth == 1){
    outlines[length(outlines)+1] = "layerF=add_layer(xs,Fsc,hidden_feats,activation_function=tf.sigmoid,dropout_function=True,lambda1=lambda1, keep_prob1=kprob)"
  }else{
    for (i in 1:FFdepth){
      if (i == 1){
        outlines[length(outlines)+1] = paste0("layer",as.character(i),"=add_layer(xs,Fsc,hidden_feats,activation_function=tf.sigmoid,dropout_function=True,lambda1=lambda1, keep_prob1=kprob)")
      }else if (i < FFdepth){
        outlines[length(outlines)+1] = paste0("layer",as.character(i),"=add_layer(layer",as.character(i-1),",hidden_feats,hidden_feats,activation_function=tf.sigmoid,dropout_function=True,lambda1=lambda1, keep_prob1=kprob)")
      }else{
        outlines[length(outlines)+1] = paste0("layerF=add_layer(layer",as.character(i-1),",hidden_feats,hidden_feats,activation_function=tf.sigmoid,dropout_function=True,lambda1=lambda1, keep_prob1=kprob)")
      }
    }
  }
  fout = file(paste0(tmpDir,model_type,'MTL_p2.py'))
  writeLines(outlines,fout)
  close(fout)
  outlines = c()
  outlines[length(outlines)+1] = "#***********************************************************************"
  outlines[length(outlines)+1] = "# extracting coefficients from TF graph"
  if(model_type=='ClassClass' || model_type=='ClassCox'){additional_layers = 3}
  else{additional_layers=0}
  for (i in 1:(FFdepth+1+additional_layers)){
    if (i==1){
      outlines[length(outlines)+1] = paste0("Theta1 = sess.run(tf.get_collection(tf.GraphKeys.TRAINABLE_VARIABLES,'Variable:0'))[0]")
      outlines[length(outlines)+1] = paste0("Bias1 = sess.run(tf.get_collection(tf.GraphKeys.TRAINABLE_VARIABLES,'Variable_1:0'))[0]")
    }else{
      outlines[length(outlines)+1] = paste0("Theta",as.character(i)," = sess.run(tf.get_collection(tf.GraphKeys.TRAINABLE_VARIABLES,'Variable_",as.character((2*(i-1))),":0'))[0]")
      outlines[length(outlines)+1] = paste0("Bias",as.character(i)," = sess.run(tf.get_collection(tf.GraphKeys.TRAINABLE_VARIABLES,'Variable_",as.character((2*(i-1))+1),":0'))[0]")
    }
  }
  outlines[length(outlines)+1] = "#***********************************************************************"
  outlines[length(outlines)+1] = "# Saving model coefficients to files"
  for (i in 1:(FFdepth+1+additional_layers)){
    outlines[length(outlines)+1] = paste0("np.savetxt(data_folder+'Theta",as.character(i),".csv',Theta",as.character(i),",delimiter=',')")
    outlines[length(outlines)+1] = paste0("np.savetxt(data_folder+'Bias",as.character(i),".csv',Bias",as.character(i),",delimiter=',')")
  }
  outlines[length(outlines)+1] = "with open(data_folder+'Activations.csv','w') as f:"
  for (i in 1:FFdepth){
    outlines[length(outlines)+1] = "    f.write('sigmoid\\n')"
  }
  if(model_type == 'BlankCox'){
    outlines[length(outlines)+1] = "    f.write('sigmoid\\n')"
  }else{
    outlines[length(outlines)+1] = "    f.write('softmax\\n')"
  }
  if (model_type == 'ClassBlank' || model_type == 'BlankClass' || model_type == 'BlankCox'){
    #outlines[length(outlines)+1] = "    f.write('softmax\\n')"
  }else{
    if (model_type == 'ClassClass'){
      outlines[length(outlines)+1] = "    f.write('softmax\\n')"
    }else{
      outlines[length(outlines)+1] = "    f.write('sigmoid\\n')"
    }
    outlines[length(outlines)+1] = "    f.write('sigmoid\\n')"
    if (model_type == 'ClassClass'){
      outlines[length(outlines)+1] = "    f.write('softmax\\n')"
    }else{
      outlines[length(outlines)+1] = "    f.write('sigmoid\\n')"
    }
  }
  fout = file(paste0(tmpDir,model_type,'MTL_p4.py'))
  writeLines(outlines,fout)
  close(fout)
  system(paste0("cat ",tmpDir,model_type,"MTL_p1.py ",tmpDir,model_type,"MTL_p2.py ",tmpDir,model_type,"MTL_p3.py ",tmpDir,model_type,"MTL_p4.py > ",tmpDir,model_type,"MTL.py"))
}

# Make python executable for densenet implementation
makeExec2 <- function(tmpDir,FFdepth,model_type){
  if (model_type != 'ClassClass' && model_type != 'ClassCox' && model_type != 'ClassBlank' && model_type != 'BlankClass' && model_type!='BlankCox'){
    stop("Please specify either 'BlankClass', 'ClassBlank', 'BlankCox', ClassClass' or 'ClassCox' for the model_type")
  }
  system(paste0('cp ',DEGAS.toolsPath,model_type,'MTL_p1.py ',tmpDir))
  system(paste0('cp ',DEGAS.toolsPath,model_type,'MTL_p3.py ',tmpDir))
  outlines = c()
  if (FFdepth == 1){
    outlines[length(outlines)+1] = "layerF=add_layer(xs,Fsc,hidden_feats,activation_function=tf.sigmoid,dropout_function=True,lambda1=lambda1, keep_prob1=kprob)"
  }else{
    inpsz_str = "Fsc"
    for (i in 1:FFdepth){
      if (i == 1){
        outlines[length(outlines)+1] = paste0("layer",as.character(i),"=add_layer(xs,",inpsz_str,",hidden_feats,activation_function=tf.sigmoid,dropout_function=True,lambda1=lambda1, keep_prob1=kprob)")
      }else if (i < FFdepth){
        inpsz_str = paste0(inpsz_str,"+hidden_feats")
        layer_str = "tf.concat([xs"
        for(j in 1:i){
          if (j==i){
            layer_str = paste0(layer_str,'],1)')
          }else{
            layer_str = paste0(layer_str,',layer',as.character(j))
          }
        }
        outlines[length(outlines)+1] = paste0("layer",as.character(i),"=add_layer(",layer_str,",",inpsz_str,",hidden_feats,activation_function=tf.sigmoid,dropout_function=True,lambda1=lambda1, keep_prob1=kprob)")
      }else{
        inpsz_str = paste0(inpsz_str,"+hidden_feats")
        layer_str = "tf.concat([xs"
        for(j in 1:FFdepth){
          if (j==i){
            layer_str = paste0(layer_str,'],1)')
          }else{
            layer_str = paste0(layer_str,',layer',as.character(j))
          }
        }
        outlines[length(outlines)+1] = paste0("layerF=add_layer(",layer_str,",",inpsz_str,",hidden_feats,activation_function=tf.sigmoid,dropout_function=True,lambda1=lambda1, keep_prob1=kprob)")
      }
    }
  }
  fout = file(paste0(tmpDir,model_type,'MTL_p2.py'))
  writeLines(outlines,fout)
  close(fout)
  outlines = c()
  outlines[length(outlines)+1] = "#***********************************************************************"
  outlines[length(outlines)+1] = "# extracting coefficients from TF graph"
  if(model_type=='ClassClass' || model_type=='ClassCox'){additional_layers = 3}
  else{additional_layers=0}
  for (i in 1:(FFdepth+1+additional_layers)){
    if (i==1){
      outlines[length(outlines)+1] = paste0("Theta1 = sess.run(tf.get_collection(tf.GraphKeys.TRAINABLE_VARIABLES,'Variable:0'))[0]")
      outlines[length(outlines)+1] = paste0("Bias1 = sess.run(tf.get_collection(tf.GraphKeys.TRAINABLE_VARIABLES,'Variable_1:0'))[0]")
    }else{
      outlines[length(outlines)+1] = paste0("Theta",as.character(i)," = sess.run(tf.get_collection(tf.GraphKeys.TRAINABLE_VARIABLES,'Variable_",as.character((2*(i-1))),":0'))[0]")
      outlines[length(outlines)+1] = paste0("Bias",as.character(i)," = sess.run(tf.get_collection(tf.GraphKeys.TRAINABLE_VARIABLES,'Variable_",as.character((2*(i-1))+1),":0'))[0]")
    }
  }
  outlines[length(outlines)+1] = "#***********************************************************************"
  outlines[length(outlines)+1] = "# Saving model coefficients to files"
  for (i in 1:(FFdepth+1+additional_layers)){
    outlines[length(outlines)+1] = paste0("np.savetxt(data_folder+'Theta",as.character(i),".csv',Theta",as.character(i),",delimiter=',')")
    outlines[length(outlines)+1] = paste0("np.savetxt(data_folder+'Bias",as.character(i),".csv',Bias",as.character(i),",delimiter=',')")
  }
  outlines[length(outlines)+1] = "with open(data_folder+'Activations.csv','w') as f:"
  for (i in 1:FFdepth){
    outlines[length(outlines)+1] = "    f.write('sigmoid\\n')"
  }
  if(model_type == 'BlankCox'){
    outlines[length(outlines)+1] = "    f.write('sigmoid\\n')"
  }else{
    outlines[length(outlines)+1] = "    f.write('softmax\\n')"
  }
  if (model_type == 'ClassBlank' || model_type == 'BlankClass' || model_type == 'BlankCox'){
    #outlines[length(outlines)+1] = "    f.write('softmax\\n')"
  }else{
    if (model_type == 'ClassClass'){
      outlines[length(outlines)+1] = "    f.write('softmax\\n')"
    }else{
      outlines[length(outlines)+1] = "    f.write('sigmoid\\n')"
    }
    outlines[length(outlines)+1] = "    f.write('sigmoid\\n')"
    if (model_type == 'ClassClass'){
      outlines[length(outlines)+1] = "    f.write('softmax\\n')"
    }else{
      outlines[length(outlines)+1] = "    f.write('sigmoid\\n')"
    }
  }
  fout = file(paste0(tmpDir,model_type,'MTL_p4.py'))
  writeLines(outlines,fout)
  close(fout)
  system(paste0("cat ",tmpDir,model_type,"MTL_p1.py ",tmpDir,model_type,"MTL_p2.py ",tmpDir,model_type,"MTL_p3.py ",tmpDir,model_type,"MTL_p4.py > ",tmpDir,model_type,"MTL.py"))
}

# train model wrapper function
runCCMTL <- function(scExp,scLab,patExp,patLab,tmpDir,model_type,architecture,FFdepth){
  system('pwd')
  system(paste0('rm -rf ',tmpDir))
  system(paste0('mkdir ',tmpDir))
  message(paste0(as.character(FFdepth), "-layer ", architecture, " ",model_type, " DEGAS model"))
  if(architecture=="DenseNet"){
    makeExec2(tmpDir, FFdepth, model_type)
  }else if(architecture=="Standard"){
    makeExec(tmpDir, FFdepth, model_type)
  }else{
    stop('Incorrect architecture argument')
  }
  writeInputFiles(scExp,scLab,patExp,patLab,tmpDir)
  message(checkForPy())
  system(paste0(DEGAS.pyloc," ",tmpDir,model_type,"MTL.py ", tmpDir, " ",DEGAS.train_steps," ",DEGAS.scbatch_sz," ",DEGAS.patbatch_sz," ",DEGAS.hidden_feats," ",DEGAS.do_prc," ",DEGAS.lambda1," ",DEGAS.lambda2," ",DEGAS.lambda3," ",DEGAS.seed))
  ccModel1 = readOutputFiles(tmpDir,model_type,architecture)
  return(ccModel1)
}

# Bootstrap aggregation wrapper for model training
runCCMTLBag <- function(scExp,scLab,patExp,patLab,tmpDir,model_type,architecture,FFdepth,Bagdepth){
  orig_degas_seed = DEGAS.seed
  out <- list()
  for(i in 1:Bagdepth){
    DEGAS.seed <<- orig_degas_seed + (i-1)
    out[[i]] <- runCCMTL(scExp,scLab,patExp,patLab,tmpDir,model_type,architecture,FFdepth)
  }
  DEGAS.seed <<- orig_degas_seed
  return(out)
}

# Make predictions based on a trained DEGAS model
predClass <- function(ccModel1,Exp,scORpat){
  if(ccModel1@Architecture=="DenseNet"){
    return(predClass2(ccModel1,Exp,scORpat))
  }else if (ccModel1@Architecture=="Standard"){
    return(predClass1(ccModel1,Exp,scORpat))
  }else{
    stop("Incorrect architecture argument")
  }
}

# Prediction from trained standard architecture model
predClass1 <- function(ccModel1,Exp,scORpat){
  Z = Exp
  rm(Exp)
  if (ccModel1@Model_type=='BlankClass' || ccModel1@Model_type=='ClassBlank' || ccModel1@Model_type=='ClassBlank'){
    for (i in 1:(ccModel1@Depth)){
      calcZ = paste0(ccModel1@Activation[[i]],"(sweep((as.matrix(Z) %*% ccModel1@Theta[[",as.character(i),"]]),2,ccModel1@Bias[[",as.character(i),"]],'+'))")
      Z = eval(parse(text=calcZ))
    }
    return(Z)
  }else{
    for (i in 1:(ccModel1@Depth-4)){
      calcZ = paste0(ccModel1@Activation[[i]],"(sweep((as.matrix(Z) %*% ccModel1@Theta[[",as.character(i),"]]),2,ccModel1@Bias[[",as.character(i),"]],'+'))")
      Z = eval(parse(text=calcZ))
    }
  }
  if (toupper(scORpat)=='SC'){
    calcPred = paste0(ccModel1@Activation[[ccModel1@Depth-3]],"(sweep((Z %*% ccModel1@Theta[[",as.character(ccModel1@Depth-3),"]]),2,ccModel1@Bias[[",as.character(ccModel1@Depth-3),"]],'+'))")
  }else if (toupper(scORpat)=='PAT'){
    calcPred = paste0(ccModel1@Activation[[ccModel1@Depth-2]],"(sweep((Z %*% ccModel1@Theta[[",as.character(ccModel1@Depth-2),"]]),2,ccModel1@Bias[[",as.character(ccModel1@Depth-2),"]],'+'))")
  }else{
    stop("Incorrect prediction argument. Please use 'sc' or 'pat'")
  }
  return(eval(parse(text=calcPred)))
}

# Prediction from trained densenet architecture model
predClass2 <- function(ccModel1,Exp,scORpat){
  Z = Exp
  rm(Exp)
  if (ccModel1@Model_type=='BlankClass' || ccModel1@Model_type=='ClassBlank' || ccModel1@Model_type=='BlankCox'){
    for (i in 1:(ccModel1@Depth)){
      calcZ = paste0(ccModel1@Activation[[i]],"(sweep((as.matrix(Z) %*% ccModel1@Theta[[",as.character(i),"]]),2,ccModel1@Bias[[",as.character(i),"]],'+'))")
      if(i<ccModel1@Depth-1){
        Z = cbind(Z,eval(parse(text=calcZ)))
      }else{
        Z = eval(parse(text=calcZ))
      }
    }
    return(Z)
  }else{
    for (i in 1:(ccModel1@Depth-4)){
      calcZ = paste0(ccModel1@Activation[[i]],"(sweep((as.matrix(Z) %*% ccModel1@Theta[[",as.character(i),"]]),2,ccModel1@Bias[[",as.character(i),"]],'+'))")
      if(i<ccModel1@Depth-4){
        Z = cbind(Z,eval(parse(text=calcZ)))
      }else{
        Z = eval(parse(text=calcZ))
      }
    }
  }
  if (toupper(scORpat)=='SC'){
    calcPred = paste0(ccModel1@Activation[[ccModel1@Depth-3]],"(sweep((Z %*% ccModel1@Theta[[",as.character(ccModel1@Depth-3),"]]),2,ccModel1@Bias[[",as.character(ccModel1@Depth-3),"]],'+'))")
  }else if (toupper(scORpat)=='PAT'){
    calcPred = paste0(ccModel1@Activation[[ccModel1@Depth-2]],"(sweep((Z %*% ccModel1@Theta[[",as.character(ccModel1@Depth-2),"]]),2,ccModel1@Bias[[",as.character(ccModel1@Depth-2),"]],'+'))")
  }else{
    stop("Incorrect prediction argument. Please use 'sc' or 'pat'")
  }
  return(eval(parse(text=calcPred)))
}

# Predict from bootstrap aggregated models
predClassBag <- function(ccModel,Exp,scORpat){
  out = list()
  for(i in 1:length(ccModel)){
    out[[i]] <- predClass(ccModel[[i]],Exp,scORpat)
  }
  out = Reduce("+", out) / length(out)
  return(out)
}
  
# Predict patient class from proportions of single cell classes
predPatClassFromSCClass <- function(ccModel1,Exp){
  Z1 = sigmoid(sweep((Exp %*% ccModel4@Theta4),2,ccModel1@Bias4,'+'))
  return(softmax(sweep((Z1 %*% ccModel1@Theta5),2,ccModel1@Bias5,'+')))
}

#***************************************************************
# Other functions

# Generates sets for k fold cross validation
splitKfoldCV <- function(N,k){
  if(k<3){
    stop("Please use 3 or more folds")
  }
  Idx = as.numeric(sample(1:N,N,replace=FALSE))
  sz = rep(floor(N/k),k)
  rem = N-sum(sz)
  if(rem>0){
    cntr=0
    for(i in 1:rem){
      if(cntr==k){
        cntr = 1
      }else{
        cntr = cntr+1
      }
      sz[cntr] = sz[cntr]+1
    }
  }
  cntr = 0
  grpIdx = list()
  for (i in 1:k){
    grpIdx[[i]] = Idx[(cntr+1):(cntr+sz[i])]
    cntr = cntr + sz[i]
  }
  return(grpIdx)
}

# Get a feature vector from a dataframe
getFeat = function(vec,df,colm,colo){
  tmp = vec;
  for (i in 1:length(vec)){
    tmp[i] = df[which(df[,colm]==vec[i])[1],colo]
  }
  return(tmp)
}

# returns duplicate row names to remove
remDupIdx <- function(X,dup_rnames,rnames){
  rem = c()
  for (dup_rname in dup_rnames){
    tmp = which(rnames==dup_rname)
    Xtmp = X[tmp,]
    Xmean = rowMeans(as.matrix(Xtmp[,2:dim(Xtmp)[2]]))
    Xmean[is.na(Xmean)] = 0
    Xmean = abs(Xmean)
    rem = c(rem,tmp[which(Xmean!=max(Xmean,na.rm=TRUE))])
  }
  return(rem)
}

#***************************************************************
# Post-processing functions

# Quantile normalization
quantNorm <- function(df,center='median',rescale=TRUE,rescale_mult=1e4){
  df_rank <- apply(df,2,rank,ties.method="min")
  df_sorted <- data.frame(apply(df, 2, sort))
  df_mean <- apply(df_sorted, 1, mean)
  index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }
  df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
  rownames(df_final) <- rownames(df)
  meds = eval(parse(text=paste0("apply(df_final, 2, ",center,")")))
  df_final = sweep(df_final, 2, meds, '-')
  if(rescale){
    df_final[df_final>0] = eval(parse(text=paste0("log2(",rescale_mult,"*df_final[df_final>0])")))
    df_final[df_final<0] = eval(parse(text=paste0("-log2(-",rescale_mult,"*df_final[df_final<0])")))
    meds = eval(parse(text=paste0("apply(df_final, 2, ",center,")")))
    df_final = sweep(df_final, 2, meds, '-')
  }
  return(df_final)
}

# Return euclidean distance between two points
euclDist <- function(loc1,loc2){
  return(sqrt(sum((loc1 - loc2)^2)))
}

# Return a matrix of all pairwise distances
pairDist <- function(locs){
  N = dim(locs)[1]
  out = matrix(NA,N,N)
  for(i in 1:N){
    for(j in 1:i){
      out[i,j] = out[j,i] = euclDist(locs[i,],locs[j,])
    }
  }
  return(out)
}

# Return k-nearest-neighbor smoothed probabilites
knnSmooth <- function(probs,locs,k=5){
  out = probs
  dists = pairDist(locs)
  if(class(probs)=="numeric"){
    N = length(probs)
    for(i in 1:N){
      idx = order(dists[i,])
      out[i] = mean(probs[idx[1:k]])
    }
  }else{
    N = dim(locs)[1]
    for(i in 1:N){
      idx = order(dists[i,])
      out[i,] = colMeans(probs[idx[1:k],])
    }
  }
  return(out)
}

# return random sample (of s) which is evenly distributed across sample groups (g)
# where each group has n samples.
# Note: If a group has less than n samples, then all samples in that group are used.
evenSamp <- function(s,g,n){
  groups = unique(g)
  out = list()
  for (group in groups){
    if(sum(g==group)>=n){
      out[[group]] = sample(s[g==group],n,replace=FALSE)
    }else if(sum(g==group)>0){
      out[[group]] = sample(s[g==group],sum(g==group),replace=FALSE)
    }else{
      #Adding nothing
    }
  }
  out = unlist(out)
  return(out)
}

# Convert DEGAS output [0,1] to an association [-1,1]
toCorrCoeff <- function(probs){
  k=dim(probs)[2]
  if(k<2 || is.null(k)){k=2}
  l=2
  return(2*((probs-1/k)/(l-l/k) + 1/l)-1)
}
