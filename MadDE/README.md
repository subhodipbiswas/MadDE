## MadDE code release
*******************************************
#### Author(s): S. Biswas and D. Saha 
#### Date: 17/Jun/2021
*******************************************
##### Change Log: 
    17/Jun/2021:    MadDE v1.0.0 - Initial Code-Release
*******************************************
### Description
*******************************************
This is an MATLAB implementation of the **M**ultiple **ad**aptation
**D**ifferential **E**volution (**MadDE**)
algorithm that was accepted at the CEC 2021 Special Session and 
Competition on single-objective bound-constrained numerical optimization problem

### Papers
*******************************************
For more details on MadDE, see the following paper(s):

* **MadDE**:
  ```S. Biswas, D. Saha, S. De, A. D. Cobb, S. Das and B. A. Jalaian, "A multiple adaptation-based Differential Evolutionfor global numerical optimization," 2021 IEEE Congress on Evolutionary Computation (CEC), Kraków, Poland, 2021.```

* **LSHADE**:
  Ryoji Tanabe and Alex Fukunaga: Improving the Search Performance of SHADE Using Linear Population Size Reduction,  Proc. IEEE Congress on Evolutionary Computation (CEC-2014), Beijing, July, 2014.

* **JADE**:
  J. Zhang, A.C. Sanderson: JADE: Adaptive differential evolution with optional external archive,” IEEE Trans Evol Comput, vol. 13, no. 5, pp. 945–958, 2009


### Source Code
*******************************************
For this package, we downloaded source codes from multiple repositories publicly available. Our intention has been to contribute towards open-sourcing the codes for replication and experiment. 
The code implementation is based on the following websites:

* JADE's source code from Dr. Q. Zhang's [website](http://dces.essex.ac.uk/staff/qzhang).
* LSHADE's source code from Dr. R. Tanabe's [website](https://ryojitanabe.github.io/publication).

The source code uses a random seeds present inside folder `input_data/Rand_Seeds.txt` along with `Mersenne Twister` random number generator. 
A simple example would be

`% set the random number generator settings` \
`seed = 100                 % Set a seed` \
`s = rng(seed, 'twister')   % store the generator settings`

     Type: 'twister'
     Seed: 100
     State: [625×1 uint32]

Using this seed helps in replication of the experimental results.

### Environment
*******************************************
System configurations in our experimental environment:

    - OS: Windows 10 Pro
    - CPU: AMD Ryzen 7 (3.9GHz)
    - RAM: 16GB
    - Language: Matlab


### Execution Instructions
*******************************************

Run the file ``main.m`` to simulate the sample runs of **MadDE** on CEC 2021 benchmark and store the output in the ``Results`` folder.

### Code Modification
*******************************************
Dimension size and all the parameter settings of the algorithm are easily changeable by rewriting source code in *MadDE.m*.

### Queries
*******************************************
If you have any questions, please feel free to contact us (sub17was [at] gmail.com, debanjansh [at] gmail.com).
