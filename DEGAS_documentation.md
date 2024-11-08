# DEGAS Function Documentation

## 1) Preprocessing

**normFunc( x )**<br>
**Inputs**<br>
**x:** A vector of numbers to convert to z-scores. Generally each cell or patient sample would use this conversion in the DEGAS workflow.<br>
**Output**<br>
A vector of z-scores based on the original numeric vector<br>
<br>
**scaleFunc( x )**<br>
**Inputs**<br>
**x:** A vector of number to scale to [0,1]. Generally each cell or patient sample would use this conversion after conversion to z-scores using normFunc.<br>
**Output**<br>
A vector of [0,1] scaled values<br>
<br>
**preprocessCounts( X )**<br>
**Inputs**<br>
**X:** A matrix or dataframe of counts where the rows represent genes and the columns represent cells/samples<br>
**Output**<br>
The output is matrix of normalized (log2 and zscore, i.e. log2  and normFunc) and scaled expression ([0-1] scaling, i.e. ScaleFunc) values to use as input for DEGAS model training. The resulting matrix is transposed such that the rows are cells/samples and the columns are genes.<br>
<br>
**initDEGAS()**<br>
This function must be run to initialize the DEGAS parameters before training the model. The following functions can be used to change the initialized DEGAS parameters to the user’s preference.<br>
**Inputs**<br>
None<br>
**Output**<br>
None<br>
<br>
**setPython( path2python )**<br>
**Inputs**<br>
**path2python:** A string containing the path to the location of the preferred python executable. This is important if the tensorflow package is only available for a specific version of python on the user’s computer.<br>
**Output**<br>
None<br>
<br>
**set_training_steps( inp )**<br>
**Inputs**<br>
**inp:** The number of training steps during DEGAS model training.<br>
**Output**<br>
None<br>
<br>
**set_single_cell_batch_size( inp )**<br>
**Inputs**<br>
**inp:** The number of cells to use in the minibatch<br>
**Output**<br>
None<br>
<br>
**set_patient_batch_size( inp )**<br>
**Inputs**<br>
**inp:** The number of patients to use in the minibatch<br>
**Output**<br>
None<br>
<br>
**set_hidden_feature_number( inp )**<br>
**Inputs**<br>
**inp:** The number of features in each hidden layer of the DEGAS model.<br>
**Output**<br>
None<br>
<br>
**set_dropout_keep_fraction( inp )**<br>
**Inputs**<br>
**inp:** The percentage of nodes to keep during dropout.<br>
**Output**<br>
None<br>
<br>
**set_l2_regularization_term( inp )**<br>
**Inputs**<br>
**inp:** The regularization term weight lambda for L2 regularization.<br>
**Output**<br>
None<br>
<br>
**set_patient_loss_term( inp )**<br>
**Inputs**<br>
**inp:** The term weight lambda for the patient label loss.<br>
**Output**<br>
None<br>

**set_MMD_loss_term( inp )**<br>
**Inputs**<br>
**inp:** The term weight lambda for the MMD loss.<br>
**Output**<br>
None<br>

**set_seed_term( inp )**<br>
**Inputs**<br>
**inp:** The seed to set for the DEGAS model.<br>
**Output**<br>
None<br>

## 2) Model training and prediction

**runCCMTLBag( scExp, scLab, patExp, patLab, tmpDir, model_type, architecture, FFdepth, Bagdepth )**<br>
**Inputs**<br>
**scExp:** A matrix (rows=cells, columns=genes) of expression values from scRNA-seq data. The columns (i.e. genes) should be in the same order as the patExp matrix.<br>
**scLab:** A matrix (rows=cells, columns=labels) of binary cell labels corresponding to the scExp cells. This is used in the case where cell type or cell cluster labels exist for the scRNA-seq data. Each row is onehot encoded such that each row contains a one in a single column and zeros in the rest of the columns. In the case that no such labels exist NULL should be passed to this argument.<br>
**patExp:** A matrix (rows=patient samples, columns=genes) of expression values from bulk expression data. The columns (i.e. genes) should be in the same order as the scExp matrix.<br>
**patLab:** A matrix (rows=patient samples, columns=labels) of patient sample labels corresponding to the patExp samples. In the case of classification tasks there should exist two or more columns where each row contains a single one with all other columns in that row containing a zero. In the case of survival tasks there should be two columns. The first column will contain times and the second column should contain event status (1=event, 0=censored). In the case that no such labels exist NULL should be passed to this argument.<br>
**tmpDir:** This string argument indicates the path to the tmp directory used for processing (e.g. “~/Desktop/tmp/”)<br>
**model_type:** This string argument indicates the type of model used. Options include “BlankClass”, “ClassBlank”, “ClassClass”, “BlankCox”, and “ClassCox”. This indication should match the input matrices for scLab and patLab. As an example, for a BlankClass model, scLab should be NULL and patLab should be a onehot matrix of patient sample labels.<br>
**architecture:** This string argument can either be “Standard” or “DenseNet”. This specification indicates which type of model to train where “Standard” is a feed forward network and DenseNet is a dense net network.<br>
**FFdepth:** This numeric argument indicates the number of layers in the model. This number must be greater than or equal to 1. Please note that very large numbers may take long times to run or use all of the available RAM. Also, note that Standard and DenseNet models are the same if this argument is set to 1.<br>
**Bagdepth:** This numeric argument indicates the number of times to bootstrap aggregate the models. This argument must be greater than or equal to 1. Please not that if equal to 1, there is no bootstrap aggregation.<br>
**Output**<br>
A trained bootstrap aggregated DEGAS model.<br>
<br>
**predClassBag( ccModel, Exp, scORpat )**<br>
**Inputs**<br>
**ccModel:** This is a bootstrap aggregated DEGAS model that was previously trained with the runCCMTLBag function.<br>
**Exp:** This is a matrix (rows=cells/samples, columns=genes) of expression values from which labels will be predicted. The order of the column genes should be identical for the prediction dataset and the training dataset.<br>
**scORpat:** This is a string (“sc” or “pat”) to tell the algorithms whether to output the predicted patient sample label (“pat”) or the predicted cellular label (“sc”). It is worth noting that for BlankClass and BlankCox models the cellular label can not be predicted. For ClassBlank models the patient sample label cannot be predicted.<br>
**Output**<br>
A matrix of label probabilities based on the new input data, the trained bootstrap aggregated DEGAS model, and the output label of choice.<br>
<br>
**runCCMTL( scExp, scLab, patExp, patLab, tmpDir, model_type, architecture, FFdepth)**<br>
See runCCMTLBag. The only difference is that this function does not train a bootstrap aggregated model. The inputs and outputs are the same as runCCMTLBag except that the model is not bootstrap aggregated.<br>
<br>
**predClass( ccModel1, Exp, scORpat)**<br>
See predClassBag. The only difference is that it predicts outcomes from a singular DEGAS model, i.e. not bootstrap aggregated. The inputs and outputs are the same as predClassBag except that the model is not bootstrap aggregated.<br>

## 3) Post processing

**knnSmooth( probs, locs, k )**<br>
**Inputs**<br>
**probs:** A numeric matrix (rows=cells/samples, columns=labels) which will be smoothed. For the purposes of the general DEGAS workflow, this will be the association matrix which can be smoothed to remove some of the random noise in the labels.<br>
**locs:** A numeric matrix of locations for the row in the probs matrix. For the purposes of the general DEGAS workflow, this will generally be PCA/tSNE/UMAP coordinates of the sample association matrix.<br>
**k:** A number indicating the number of neighbors to smooth with. Default is 5.<br>
**Output**<br>
A numeric matrix with the same dimensions as probs that has been smoothed based on it’s neighbors.<br>
<br>
**toCorrCoeff( probs )**<br>
**Inputs**<br>
**probs:** A numeric matrix (rows=cells/samples, columns=labels) which will be smoothed. For the purposes of the general DEGAS workflow, this will be the association matrix which can be smoothed to remove some of the random noise in the labels.<br>
**Output**<br>
A numeric matrix where the probability matrix values [0,1] have been converted to associations [-1,1].<br>
