import tensorflow as tf				#NEED
if int(tf.__version__[0])>1:
	import tensorflow.compat.v1 as tf
	tf.disable_v2_behavior()  
from functools import partial		#NEED
#slim = tf.contrib.slim
import numpy as np					#NEED
#import pandas as pd
import math							#NEED
#from scipy import stats				
#from sklearn import preprocessing
#import scipy.io as sio
#from os import listdir
import os							#NEED
import sys							#NEED
from os.path import isfile, join	#NEED

#***********************************************************************
# Gaussian kernel matrix

def compute_pairwise_distances(x, y):
    if not len(x.get_shape()) == len(y.get_shape()) == 2:
        raise ValueError('Both inputs should be matrices.')
    if x.get_shape().as_list()[1] != y.get_shape().as_list()[1]:
        raise ValueError('The number of features should be the same.')
    norm = lambda x: tf.reduce_sum(tf.square(x), 1)
    return tf.transpose(norm(tf.expand_dims(x, 2) - tf.transpose(y)))

def gaussian_kernel_matrix(x, y, sigmas):
    beta = 1. / (2. * (tf.expand_dims(sigmas, 1)))
    dist = compute_pairwise_distances(x, y)
    s = tf.matmul(beta, tf.reshape(dist, (1, -1)))
    return tf.reshape(tf.reduce_sum(tf.exp(-s), 0), tf.shape(dist))

#***********************************************************************
# Add weights to hidden layer
def get_weight(shape, lambda1): 
	var = tf.Variable(tf.random_normal(shape), dtype=tf.float32)
	if tf.__version__[0]==1:
		tf.add_to_collection('losses', tf.contrib.layers.l2_regularizer(lambda1)(var))
	else:
		tf.add_to_collection('losses', tf.keras.regularizers.L2(lambda1)(var))  
	return var

#***********************************************************************
# Add hidden layer to network
def add_layer(input,in_size,out_size,activation_function=None,dropout_function=False,lambda1=0, keep_prob1=1):
	Weights= get_weight([in_size,out_size], lambda1) 
	biases=tf.Variable(tf.zeros([1,out_size])+0.1)
	Wx_plus_b=tf.matmul(input,Weights)+biases
	if dropout_function==True:
		Wx_plus_b=tf.nn.dropout(Wx_plus_b,keep_prob=keep_prob1)
	else:
		pass
	if activation_function is None:
		outputs=Wx_plus_b
	else:
		outputs=activation_function(Wx_plus_b)
	return outputs

#***********************************************************************
# Generate R matrix
def Rmatrix(surv):
	out = np.zeros([len(surv), len(surv)], dtype=int)
	for i in range(len(surv)):
		for j in range(len(surv)):
			out[i,j] = surv[j] >= surv[i]
	return(out)

#***********************************************************************
# CORAL loss for the hidden layer output from jindongwang github
def CORAL_loss(source,target):
	d = tf.cast(source.shape[1],tf.float32)
	# source covariance
	xm = tf.reduce_mean(source, axis=0,keep_dims=True) - source
	xc = tf.matmul(tf.transpose(xm), xm)
	# target covariance
	xmt = tf.reduce_mean(target, axis=0,keep_dims=True) - target
	xct = tf.matmul(tf.transpose(xmt), xmt)
	# frobenius norm between source and target
	loss = tf.sqrt(tf.reduce_sum(tf.pow(xc - xct,2)))
	loss = loss/(4*d*d)
	return loss

#***********************************************************************
# MMD loss for the hidden layer output from jindongwang github
def maximum_mean_discrepancy(x, y, kernel=gaussian_kernel_matrix):
    """Computes the Maximum Mean Discrepancy (MMD) of two samples: x and y.
    Maximum Mean Discrepancy (MMD) is a distance-measure between the samples of
    the distributions of x and y. Here we use the kernel two sample estimate
    using the empirical mean of the two distributions.
    MMD^2(P, Q) = || \E{\phi(x)} - \E{\phi(y)} ||^2
                = \E{ K(x, x) } + \E{ K(y, y) } - 2 \E{ K(x, y) },
    where K = <\phi(x), \phi(y)>,
      is the desired kernel function, in this case a radial basis kernel.
    Args:
      x: a tensor of shape [num_samples, num_features]
      y: a tensor of shape [num_samples, num_features]
      kernel: a function which computes the kernel in MMD. Defaults to the
              GaussianKernelMatrix.
    Returns:
      a scalar denoting the squared maximum mean discrepancy loss.
    """
    with tf.name_scope('MaximumMeanDiscrepancy'):
        # \E{ K(x, x) } + \E{ K(y, y) } - 2 \E{ K(x, y) }
        cost = tf.reduce_mean(kernel(x, x))
        cost += tf.reduce_mean(kernel(y, y))
        cost -= 2 * tf.reduce_mean(kernel(x, y))
    # We do not allow the loss to become negative.
    cost = tf.where(cost > 0, cost, 0, name='value')
    return cost


def mmd_loss(source_samples, target_samples, scope=None):
    """Adds a similarity loss term, the MMD between two representations.
    This Maximum Mean Discrepancy (MMD) loss is calculated with a number of
    different Gaussian kernels.
    Args:
      source_samples: a tensor of shape [num_samples, num_features].
      target_samples: a tensor of shape [num_samples, num_features].
      scope: optional name scope for summary tags.
      Returns:
        a scalar tensor representing the MMD loss value.
    """
    sigmas = [
      1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1, 5, 10, 15, 20, 25, 30, 35, 100,
      1e3, 1e4, 1e5, 1e6
      ]
    gaussian_kernel = partial(
      gaussian_kernel_matrix, sigmas=tf.constant(sigmas))
    loss_value = maximum_mean_discrepancy(
      source_samples, target_samples, kernel=gaussian_kernel)
    loss_value = tf.maximum(1e-4, loss_value)
    assert_op = tf.Assert(tf.is_finite(loss_value), [loss_value])
    with tf.control_dependencies([assert_op]):
        tag = 'MMD Loss'
        if scope:
            tag = scope + tag
        tf.summary.scalar(tag, loss_value)
        tf.losses.add_loss(loss_value)
    return loss_value

#***********************************************************************
# Multigroup MMD
def multi_mmd(d1, d2, d3, d4):
    return((mmd_loss(d1,d2)+mmd_loss(d1,d3)+mmd_loss(d1,d4)+mmd_loss(d2,d3)+mmd_loss(d2,d4)+mmd_loss(d3,d4))/6)

#***********************************************************************
# pairwise dist from mbsariyildiz github
def squared_dist(A): 
    expanded_a = tf.expand_dims(A, 1)
    expanded_b = tf.expand_dims(A, 0)
    distances = tf.reduce_sum(tf.squared_difference(expanded_a, expanded_b), 2)
    return distances
#***********************************************************************
# pairwise dist from mbsariyildiz github
def pairwise_dist_loss (A,labels):
    D = tf.sqrt(squared_dist(A))
    L = tf.matmul(labels,tf.transpose(labels))
    return(tf.reduce_mean(D*L))

#***********************************************************************
# CIndex for hazards etc.
def CIndex(hazards,labels,survtime_all):
    concord = 0.0
    total = 0.0
    N_test = labels.shape[0]
    labels = np.asarray(labels, dtype=bool)
    for i in range(N_test):
        if labels[i] == 1:
            for j in range(N_test):
                if survtime_all[j] > survtime_all[i]:
                    total = total + 1
                    if hazards[j] < hazards[i]: concord = concord + 1
                    elif hazards[j] < hazards[i]: concord = concord + 0.5
    return(concord/total)

def genSCdat(fpath,flist):
	for f in flist:
		print(f)
		Xtmp = np.genfromtxt(join(fpath,f),delimiter='\t',skip_header = 1)
		Xtmp = Xtmp[:,1:]
		Ftmp = []
		with open(join(fpath,f),'r') as fin:
			SCtmp = fin.readline().rstrip('\n').split('\t')
			line = fin.readline()
			Ptmp = [os.path.splitext(f)[0]]*Xtmp.shape[1]
			while line != '':
				line = line.rstrip('\n').split('\t')
				Ftmp.append(line[0])
				line = fin.readline()
		if f == flist[0]:
			outdata = Xtmp
			outcell = SCtmp
			outfeat = Ftmp
			outpati = Ptmp
		else:
			outdata = np.concatenate((outdata,Xtmp),axis=1)
			outcell = list(outcell) + list(SCtmp)
			outpati = list(outpati) + list(Ptmp)
	outfeat = dict(zip(outfeat,range(len(outfeat))))
	return(outdata,outcell,outfeat,outpati)

def genSCpat(dat,patlist):
	pats = list(set(patlist))
	out = np.zeros((dat.shape[0],len(pats)))
	idx = 0;
	for pat in pats:
		out[:,idx] = np.mean(dat[:,[x for x in range(len(patlist)) if patlist[x] == pat]],axis=1)*1000
		idx = idx+1
	return(out,pats)

def resample(prc_cut,Y,train):
    add = list()
    rem = list()
    train = np.squeeze(train)
    colsums = np.sum(Y[train,:],axis=0);
    cutoff = math.ceil(np.percentile(colsums,prc_cut));
    for i in range(len(colsums)):
        if colsums[i] == 0:
            pass
        elif colsums[i] < cutoff:
            idx = np.squeeze(np.where(Y[train,i]>=1));
            choice = np.random.choice(train[idx],int(cutoff-colsums[i]))
            add = add + choice.tolist();
        elif colsums[i] == cutoff:
            pass
        else:
            idx = np.squeeze(np.where(Y[train,i]>=1));
            choice = np.random.choice(train[idx],int(colsums[i]-cutoff),replace=False)
            rem = rem + choice.tolist()
    return list(set(train)-set(rem))+add
    #return np.concatenate((list([val for val in train if val not in rem]),add));		# slower for smaller datasets

#*******TESTING BELOW!!!!********************
def resample_mixGamma(X,Y,train,nsamp,depth):
    add = list()
    train = np.squeeze(train)
    colsums = np.sum(Y[train,:],axis=0)
    samp_per_class = round(nsamp/len(colsums))
    idx = list()
    for i in range(len(colsums)):
        idx = idx + [np.squeeze(np.where(Y[train,i]>=1)).tolist()];
        if samp_per_class > colsums[i]:
            choice = np.random.choice(train[idx[i]],int(samp_per_class),replace=True)
        else:
            choice = np.random.choice(train[idx[i]],int(samp_per_class),replace=False)
        add = add + choice.tolist()
    tmpX = np.zeros([nsamp,X.shape[1]])
    tmpY = np.zeros([nsamp,Y.shape[1]])
    for i in range(nsamp):
        percBinom = np.random.gamma(shape=1,size=len(colsums))
        percBinom = percBinom/sum(percBinom)
        intBinom = np.round(percBinom*depth)
        tmpIdx = list()
        for j in range(len(colsums)):
            if int(intBinom[j]) > colsums[j]:
                tmpIdx = tmpIdx + np.random.choice(train[idx[j]], int(intBinom[j]),replace=True).tolist()
            else:
                tmpIdx = tmpIdx + np.random.choice(train[idx[j]], int(intBinom[j]),replace=False).tolist()
        tmpX[i,:] = np.mean(X[tmpIdx,:],axis=0)+1e-3
        tmpY[i,:] = intBinom/sum(intBinom)
    #scaler = preprocessing.MinMaxScaler()        # CHANGED 20201213
    #tmpX = np.transpose(scaler.fit_transform(np.transpose(zscore(tmpX,axis=0))))     # CHANGED 20201212
    #return(tmpX,tmpY)		#CHANGED 20201211
    return(np.concatenate((X[add,:],tmpX), axis=0),np.concatenate((Y[add,:],tmpY)))		#CHANGED 20201211

def intersect(lst1,lst2):
	return(list(set(lst1) & set(lst2)))

def rank(inputdat, axis=-1):
	return np.argsort(np.argsort(inputdat, axis=axis), axis=axis)


#***********************************************************************
# IO to tmp folders
data_folder = sys.argv[1]
Xsc = np.loadtxt(data_folder+'scExp.csv',delimiter=',',skiprows=1)
Ysc = np.loadtxt(data_folder+'scLab.csv',delimiter=',',skiprows=1)
Nsc = Ysc.shape[0]
Fsc = Xsc.shape[1]
Lsc = Ysc.shape[1]
idx_sc = np.arange(Nsc)
np.random.shuffle(idx_sc)

Xpat = np.loadtxt(data_folder+'patExp.csv',delimiter=',',skiprows=1)
#Ypat = np.loadtxt(data_folder+'patLab.csv',delimiter=',',skiprows=1)
Npat = Xpat.shape[0]
Fpat = Xpat.shape[1]
#Lpat = Ypat.shape[1]
idx_pat = np.arange(Npat)
np.random.shuffle(idx_pat)

#scaler = preprocessing.MinMaxScaler()
#Xsc = np.transpose(scaler.fit_transform(np.transpose(stats.zscore(Xsc,0))))
#Xpat = np.transpose(scaler.fit_transform(np.transpose(stats.zscore(Xpat,0))))

#***********************************************************************
# Hyperparameters
train_steps = int(sys.argv[2])
scbatch_sz = int(sys.argv[3])
patbatch_sz = int(sys.argv[4])
hidden_feats = int(sys.argv[5])
do_prc = float(sys.argv[6])
lambda1 = float(sys.argv[7])
lambda2 = float(sys.argv[8])
lambda3 = float(sys.argv[9])

#***********************************************************************
# Building network
kprob = tf.placeholder(tf.float32)
xs = tf.placeholder(tf.float32, [None,Fsc])
ys_sc = tf.placeholder(tf.float32, [None,Lsc])
#ys_pat = tf.placeholder(tf.float32, [None,Lpat])
es = tf.placeholder(tf.float32, [None,Lsc])
#ps = tf.placeholder(tf.float32, [None,Lpat])
lsc = tf.placeholder(tf.int32, shape=())
lpat = tf.placeholder(tf.int32, shape=())

