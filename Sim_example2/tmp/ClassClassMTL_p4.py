#***********************************************************************
# extracting coefficients from TF graph
Theta1 = sess.run(tf.get_collection(tf.GraphKeys.TRAINABLE_VARIABLES,'Variable:0'))[0]
Bias1 = sess.run(tf.get_collection(tf.GraphKeys.TRAINABLE_VARIABLES,'Variable_1:0'))[0]
Theta2 = sess.run(tf.get_collection(tf.GraphKeys.TRAINABLE_VARIABLES,'Variable_2:0'))[0]
Bias2 = sess.run(tf.get_collection(tf.GraphKeys.TRAINABLE_VARIABLES,'Variable_3:0'))[0]
Theta3 = sess.run(tf.get_collection(tf.GraphKeys.TRAINABLE_VARIABLES,'Variable_4:0'))[0]
Bias3 = sess.run(tf.get_collection(tf.GraphKeys.TRAINABLE_VARIABLES,'Variable_5:0'))[0]
Theta4 = sess.run(tf.get_collection(tf.GraphKeys.TRAINABLE_VARIABLES,'Variable_6:0'))[0]
Bias4 = sess.run(tf.get_collection(tf.GraphKeys.TRAINABLE_VARIABLES,'Variable_7:0'))[0]
Theta5 = sess.run(tf.get_collection(tf.GraphKeys.TRAINABLE_VARIABLES,'Variable_8:0'))[0]
Bias5 = sess.run(tf.get_collection(tf.GraphKeys.TRAINABLE_VARIABLES,'Variable_9:0'))[0]
Theta6 = sess.run(tf.get_collection(tf.GraphKeys.TRAINABLE_VARIABLES,'Variable_10:0'))[0]
Bias6 = sess.run(tf.get_collection(tf.GraphKeys.TRAINABLE_VARIABLES,'Variable_11:0'))[0]
Theta7 = sess.run(tf.get_collection(tf.GraphKeys.TRAINABLE_VARIABLES,'Variable_12:0'))[0]
Bias7 = sess.run(tf.get_collection(tf.GraphKeys.TRAINABLE_VARIABLES,'Variable_13:0'))[0]
#***********************************************************************
# Saving model coefficients to files
np.savetxt(data_folder+'Theta1.csv',Theta1,delimiter=',')
np.savetxt(data_folder+'Bias1.csv',Bias1,delimiter=',')
np.savetxt(data_folder+'Theta2.csv',Theta2,delimiter=',')
np.savetxt(data_folder+'Bias2.csv',Bias2,delimiter=',')
np.savetxt(data_folder+'Theta3.csv',Theta3,delimiter=',')
np.savetxt(data_folder+'Bias3.csv',Bias3,delimiter=',')
np.savetxt(data_folder+'Theta4.csv',Theta4,delimiter=',')
np.savetxt(data_folder+'Bias4.csv',Bias4,delimiter=',')
np.savetxt(data_folder+'Theta5.csv',Theta5,delimiter=',')
np.savetxt(data_folder+'Bias5.csv',Bias5,delimiter=',')
np.savetxt(data_folder+'Theta6.csv',Theta6,delimiter=',')
np.savetxt(data_folder+'Bias6.csv',Bias6,delimiter=',')
np.savetxt(data_folder+'Theta7.csv',Theta7,delimiter=',')
np.savetxt(data_folder+'Bias7.csv',Bias7,delimiter=',')
with open(data_folder+'Activations.csv','w') as f:
    f.write('sigmoid\n')
    f.write('sigmoid\n')
    f.write('sigmoid\n')
    f.write('softmax\n')
    f.write('softmax\n')
    f.write('sigmoid\n')
    f.write('softmax\n')
