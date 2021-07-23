This code implements algorithm from 
Structure learning for continuous time Bayesian networks via penalized likelihood (2020 Miasojedow, Rejchel, Cakala)

and tests its performance on simulated data.

To compile use provided Makefile. C++ 14 or higher version is required.

Use variable NUM_THREADS in utils/parameters.cpp to set number of threads which will be used for inference. 
For fast computations value of at least 50 is advised.

If macro MEMORY_HUNGRY in utils/parameters.h is defined (it is defined by default), 
the program will use more memory unefficient subroutines. This will result in a significant increase in computation time.

Use variable TRIES in CTBN.cpp to change number of tries per model.

Average scores of inference will be send to standard output. 


