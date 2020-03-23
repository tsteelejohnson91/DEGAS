#***************************************************************
# Single Cell Transfer Learning Toolbox
# Travis S Johnson and Zhi Huang

#***************************************************************
# Initializing the package

initDEGAS <- function(){
  DEGAS.pyloc_ <<- "python"
  DEGAS.toolsPath <<- paste0(.libPaths()[1],"/DEGAS/tools/")
  DEGAS.train_steps <<- 2000
  DEGAS.scbatch_sz <<- 200
  DEGAS.patbatch_sz <<- 50
  DEGAS.hidden_feats <<- 50
  DEGAS.do_prc <<- 0.5
  DEGAS.lambda1 <<- 3.0
  DEGAS.lambda2 <<- 3.0
  DEGAS.lambda3 <<- 3.0
}

#***************************************************************
# General utility functions

getExtension <- function(file){
  ex <- strsplit(basename(file), split="\\.")[[1]]
  return(ex[length(ex)])
}

checkOS <- function(){
  return(Sys.info()['sysname'])
}

checkForPy <- function(){
  return(system(paste0(DEGAS.pyloc," -V")))
}

checkForTF <- function(){
  return(system("python -c 'import tensorflow'"))
}

setPython <- function(path2python){
  DEGAS.pyloc <<- path2python
  #Sys.setenv(PATH = paste(c(path2python,Sys.getenv("PATH")),collapse = .Platform$path.sep))
}

set_training_steps <- function(inp){
  DEGAS.train_steps <<- inp
}

set_single_cell_batch_size <- function(inp){
  DEGAS.scbatch_sz <<- inp
}

set_patient_batch_size <- function(inp){
  DEGAS.patbatch_sz <<- inp
}

set_hidden_feature_number <- function(inp){
  DEGAS.hidden_feats <<- inp
}

set_dropout_keep_fraction <- function(inp){
  DEGAS.do_prc <<- inp
}

set_l2_regularization_term <- function(inp){
  DEGAS.lambda1 <<- inp
}

set_patient_loss_term <- function(inp){
  DEGAS.lambda2 <<- inp
}

set_MMD_loss_term <- function(inp){
  DEGAS.lambda3 <<- inp
}

TFsetup <- function(){
  hasPy <- checkForPy()
  hasTF <- checkForTF()
  if (!hasPy && !hasTF){
    message("Python and TensorFlow implementations found")
  }
  if (hasPy || hasTF){
    resp <- readline(prompt="Are you currently running TensorFlow from python3? [Y|N]")
    if (toupper(resp) == 'Y'){
      resp <- readline(prompt="Please input the path to the python version fro which you are currently running TensorFlow.")
      setPython(resp)
      hasPy <- checkForPy()
      hasTF <- checkForTF()
    }else if (toupper(resp) == 'N'){

    }else{
      stop("Incorrect answer to prompt. Please input either 'Y' or 'N'.")
    }
  }
}

#***************************************************************
# ccModel class and related functions

# ccModel class
setClass("ccModel",slots=list(Bias="list",Theta="list",Activation="list",
                              Depth="numeric",Model_type="character",Architecture="character"))

normFunc <- function(x){return((x-mean(x, na.rm = T))/(sd(x, na.rm = T)+1e-3))}

scaleFunc <- function(x){return((x- min(x)) /(max(x)-min(x)+1e-3))}

centerFunc <- function(x){return(x-mean(x,na.rm=T))}

# Activation functions and utilities
sigmoid <- function(x) {
  1 / (1 + exp(-x))
}

logsumexp <- function (x) {
  y = max(x)
  y + log(sum(exp(x - y)))
}

softmax <- function (X) {
  return(t(apply(X,1,function(x) exp(x - logsumexp(x)))))
}

fromOneHot <- function(labMat){
  return(apply(labMat,1,function(x) colnames(labMat)[which(x==1)]))
}

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

probtoOneHot <- function(probMat){
  idx = apply(probMat,1,function(x) which(x==max(x)))
  probMat = probMat*0
  for (i in 1:length(idx)){
    probMat[i,idx[i]] = 1
  }
  return(probMat)
}

# IO utilities for training
writeInputFiles <- function(scExp,scLab,patExp,patLab,tmpDir){
  write.table(scExp,file=paste0(tmpDir, '/scExp.csv'), row.names=FALSE, sep=',')
  if(!is.null(scLab)){write.table(scLab,file=paste0(tmpDir, '/scLab.csv'), row.names=FALSE, sep=',')}
  write.table(patExp,file=paste0(tmpDir, '/patExp.csv'), row.names=FALSE, sep=',')
  if(!is.null(patLab)){write.table(patLab,file=paste0(tmpDir, '/patLab.csv'), row.names=FALSE, sep=',')}
}

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
  if(model_type=='ClassClass' || model_type=='CoxClass'){additional_layers = 3}
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
  outlines[length(outlines)+1] = "    f.write('softmax\\n')"
  if (model_type == 'ClassBlank' || model_type == 'BlankClass'){
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
  outlines[length(outlines)+1] = "    f.write('softmax\\n')"
  if (model_type == 'ClassBlank' || model_type == 'BlankClass'){
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

# runClassClassMTL training model function
runCCMTL <- function(scExp,scLab,patExp,patLab,tmpDir,model_type,architecture,FFdepth){
  system('pwd')
  system(paste0('rm -rf ',tmpDir))
  system(paste0('mkdir ',tmpDir))
  if(architecture=="DenseNet"){
    makeExec2(tmpDir, FFdepth, model_type)
  }else if(architecture=="Standard"){
    makeExec(tmpDir, FFdepth, model_type)
  }else{
    stop('Incorrect architecture argument')
  }
  writeInputFiles(scExp,scLab,patExp,patLab,tmpDir)
  message(checkForPy())
  system(paste0(DEGAS.pyloc," ",tmpDir,model_type,"MTL.py ", tmpDir, " ",DEGAS.train_steps," ",DEGAS.scbatch_sz," ",DEGAS.patbatch_sz," ",DEGAS.hidden_feats," ",DEGAS.do_prc," ",DEGAS.lambda1," ",DEGAS.lambda2," ",DEGAS.lambda3))
  ccModel1 = readOutputFiles(tmpDir,model_type,architecture)
  return(ccModel1)
}

predClass <- function(ccModel1,Exp,scORpat){
  if(ccModel1@Architecture=="DenseNet"){
    return(predClass2(ccModel1,Exp,scORpat))
  }else if (ccModel1@Architecture=="Standard"){
    return(predClass1(ccModel1,Exp,scORpat))
  }else{
    stop("Incorrect architecture argument")
  }
}

# Prediction for standard architecture
predClass1 <- function(ccModel1,Exp,scORpat){
  Z = Exp
  rm(Exp)
  if (ccModel1@Model_type=='BlankClass' || ccModel1@Model_type=='ClassBlank'){
    #message('A')
    for (i in 1:(ccModel1@Depth)){
      calcZ = paste0(ccModel1@Activation[[i]],"(sweep((as.matrix(Z) %*% ccModel1@Theta[[",as.character(i),"]]),2,ccModel1@Bias[[",as.character(i),"]],'+'))")
      Z = eval(parse(text=calcZ))
    }
    return(Z)
  }else{
    #message('B')
    for (i in 1:(ccModel1@Depth-4)){
      #message(i)
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

# Prediction for DenseNet architecture
predClass2 <- function(ccModel1,Exp,scORpat){
  Z = Exp
  rm(Exp)
  if (ccModel1@Model_type=='BlankClass' || ccModel1@Model_type=='ClassBlank'){
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
