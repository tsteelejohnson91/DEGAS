layer1=add_layer(xs,Fsc,hidden_feats,activation_function=tf.sigmoid,dropout_function=True,lambda1=lambda1, keep_prob1=kprob)
layer2=add_layer(tf.concat([xs,layer1],1),Fsc+hidden_feats,hidden_feats,activation_function=tf.sigmoid,dropout_function=True,lambda1=lambda1, keep_prob1=kprob)
layerF=add_layer(tf.concat([xs,layer1,layer2],1),Fsc+hidden_feats+hidden_feats,hidden_feats,activation_function=tf.sigmoid,dropout_function=True,lambda1=lambda1, keep_prob1=kprob)
