#***********************************************************************
# extracting coefficients from TF graph
Theta1 = sess.run(tf.get_collection(tf.GraphKeys.TRAINABLE_VARIABLES,'Variable:0'))[0]
Bias1 = sess.run(tf.get_collection(tf.GraphKeys.TRAINABLE_VARIABLES,'Variable_1:0'))[0]
Theta2 = sess.run(tf.get_collection(tf.GraphKeys.TRAINABLE_VARIABLES,'Variable_2:0'))[0]
Bias2 = sess.run(tf.get_collection(tf.GraphKeys.TRAINABLE_VARIABLES,'Variable_3:0'))[0]
#***********************************************************************
# Saving model coefficients to files
np.savetxt(data_folder+'Theta1.csv',Theta1,delimiter=',')
np.savetxt(data_folder+'Bias1.csv',Bias1,delimiter=',')
np.savetxt(data_folder+'Theta2.csv',Theta2,delimiter=',')
np.savetxt(data_folder+'Bias2.csv',Bias2,delimiter=',')
with open(data_folder+'Activations.csv','w') as f:
    f.write('sigmoid\n')
    f.write('softmax\n')
