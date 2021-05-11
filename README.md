# MadDE
This GitHub repository is for the **M**ultiple **ad**aptation **D**ifferential **E**volution (**MadDE**) algorithm that was accepted at the CEC 2021 Special Session and Competition on single-objective bound-constrained numerical optimization problem.

We have also provided the codes of baseline methods that we compared **MadDE** against the following:
```
AGSK
IMODE
j2020
LSHADE
LSHADE_cnEpSin
```

## SUrrogate-assisted Bayesian Hyperparameter Optimizer (SUBHO)
`SUBHO` is a hyperparameter optimizer used for finding optimal configuration of hyperparameters for a general optimizers by using the idea of Bayesian Optimization. Here, we use **SUBHO** to tune the hyperparameters of the proposed **MadDE** algorithm using MATLAB's `bayesopt` package. The framework is general and can be applied to similarly tune the hyperparameters of any optimization algorithm. For more details kindly refer to our published paper.


## Folder Structure
  ```
  ./
  │
  ├── README.md - Overview of the code repository
  │
  ├── AGSK/
  │   ├── main.m - code to execute AGSK
  │   └── Results - results of simulation run
  │
  ├── IMODE/
  │   ├── main.m - code to execute IMODE
  │   └── Results - results of simulating IMODE
  │
  ├── j2020/
  │   ├── Makefile - run this file to simulate j2020
  │   └── Results - results of simulating j2020
  │
  ├── LSHADE/
  │   ├── main.m - code to execute LSHADE
  │   └── Results - results of simulating LSHADE
  │
  ├── LSHADE_cnEpSin/
  │   ├── main.m - code to execute LSHADE_cnEpSin algorithm
  │   └── Results - results of simulating LSHADE_cnEpSin
  │
  ├── MadDE/ - code for our algorithm
  │   ├── main.m - code to execute MadDE
  │   └── Results - results of simulating MadDE
  │  
  └── calc_score.m - Calculate the score of the baseline algorithms
  │
  ├── SUBHO/ -
  │   ├── SUBHO.m - code to simulate hyperprameter tuning
  │   └── Hyperparameters/ - hyperparameter configurations (if any) found by SUBHO earlier 
  │   └── readHPTs - reads the scores corresponding the results stored in 'Hyperparameters'
  │   └── getScoreopt.m - scores a hyperparameter configuration found by SUBHO w.r.t. to manual version
  │   └── tune_hyperparameters.m - the scoring function supplied to MATLAB's "bayesopt" library for tuning parameters
  ```

### Comparing the optimizers

We have provided the results of the baseline algorithms. You can find their relative performance by running the file `calc_score.m`. **Note:** The results might slightly vary from machine to machine due to the random number generators. We suggest you simulate all the algorithms on your own machine before computing their *relative performance*.





## Citation
If you use this data/code for your work, please consider citing our paper:

```
S. Biswas, D. Saha, S. De, A. D. Cobb, S. Das and B. A. Jalaian, "Improving Differential Evolution through Bayesian
Hyperparameter Optimization," 2021 IEEE Congress on Evolutionary Computation (CEC), Kraków, Poland, 2021.
```
## Help
Should you have queries, please reach out to (sub17was, debanjansh) [at] gmail.com.
