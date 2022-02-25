### Frequent trouble shooting steps
1) If you cannot follow the examples provided on github  
please check the following.  
  
Open the terminal on your computer.
  
Check for python 3 on your computer:  
``$python3``   

If python 3 is available on your computer you will  
see the prompt when it starts:  
``>>>``  
  
If you do not see the prompt, python 3 must be  
installed.  
  
If python starts, check for the tensorflow package:  
``>>>import tensorflow as tf``  
``>>>tf.__version__``  

If a version number appears then tensorflow is  
installed on your computer. If not, then tensorflow  
needs to be installed.

The same can be performed in python for the other
required packages numpy, math, and functools.  
``>>>import numpy as np``  
``>>>np.__version__``  
``>>>import math as mt``  
``>>>mt.__version__``  
``>>>import functools as ft``  
``>>>ft.__version__``  
  
If version numbers appear after the following  
python code, then the packages are installed.  
Otherwise install them with pip3 from terminal.  
  
``$pip3 nameofpackage``

2) When using the toOneHot() function please do not use factors as input. Please convert to character before use.
3) Errors when training the model
Other common errors will occur when training the model (i.e. runCCMTLBag step). Frequently an error will read:
"cannot open file '/path/to/tmp/Activations.csv': No such file or directory"
This means that the python code did not execute properly and did not finish training the model. The most frequent cause
will be the python version you are running does not have a working version of tensorflow. After the runCCMTLBag() function
is called, it will print the version of python on your computer that it is running, e.g. Python 3.6.13 :: Anaconda, Inc.
Please ensure that when you open this version of python that you can import tensorflow:
``>>>import tensorflow as tf``
If you followed the setup steps and installed tensorflow but still get this error. Your tensorflow installation is not
associated with the python version you are running in the DEGAS package. You either need to install tensorflow for the version
of python displayed by the runCCMTLBAg output or you need to specify the path the the python version with tensorflow using the
``>set_python_path('/path/to/pythonWithTensorflow')`` function.
 
### Debugging log
