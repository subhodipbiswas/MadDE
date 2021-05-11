## LSHADE_cnEpSin
*******************************************
This package is a MATLAB/Octave source code of [LSHADE_cnEpSin](https://doi.org/10.1109/CEC.2017.7969336) which is a new version of LSHADE-EpSin.

#### Author(s): Noor H. Awad & Mostafa Z. Ali
#### Last Modified: Debanjan Saha
#### Date: 9/Jun/2021
*******************************************
#### Change Log: 
    12/Jan/2021: linked LSHADE-cnEpSin with CEC 2021 benchmark suite

### Papers
*******************************************
Please see the following papers:

* **LSHADE_cnEpSin**:
  ```Noor H. Awad, Mostafa Z. Ali, Ponnuthurai N. Suganthan, Ensemble Sinusoidal Differential Covariance Matrix Adaptation with Euclidean Neighborhood  for Solving CEC2017 Benchmark Problems, in Proc. IEEE Congr. Evol. Comput. CEC 2017, June, Donostia - San Sebastián, Spain.```

* **LSHADE-EpSin**:
  Noor H. Awad, Mostafa Z. Ali, Ponnuthurai N. Suganthan and Robert G. Reynolds: An Ensemble Sinusoidal Parameter Adaptation incorporated with L-SHADE for Solving CEC2014 Benchmark Problems, in Proc. IEEE Congr. Evol. Comput. CEC 2016, Canada, July, 2016.

* **LSHADE**:
  Ryoji Tanabe and Alex Fukunaga: Improving the Search Performance of SHADE Using Linear Population Size Reduction,  Proc. IEEE Congress on Evolutionary Computation (CEC-2014), Beijing, July, 2014.

* **JADE**:
  J. Zhang, A.C. Sanderson: JADE: Adaptive differential evolution with optional external archive,” IEEE Trans Evol Comput, vol. 13, no. 5, pp. 945–958, 2009.

### Source Code
*******************************************
For this package, Dr. Tanabe mentioned downloading JADE's source code from Dr. Q. Zhang's website (http://dces.essex.ac.uk/staff/qzhang) and modified it.
The code for the algorithm is publicly available at Suganthan's repository


### Execution Instructions
*******************************************

**Step 1**. Run the file ``make_cec2021.m`` to generate the .mex files corresponding to the .cpp files. 


**Step 2**. Run the file ``main.m`` to simulate the trial runs of **LSHADE-cnEpSin** on CEC 2021 benchmark. The results are stored in the ``Results`` folder.


### Queries
*******************************************
If you have any questions, please contact the original authors in the paper.
