# SCN

Code, datasets and results for the paper "Resilient Supply Chain Network Design using Multi-item Production" <br/> 
We are using Matlab R2022a, Yalmip Toolbox and Gurobi solver. <br/> 

Codes:<br/> 
main.m: the main script that provides some key parameters and calls the function prepare_run(). <br/> 
prepare_run.m: the function that reads in the datasets, prepares the cost parameters, and calls the function concave_minimization(). <br/> 
concave_minimization.m: the function that initializes and updates the set $S_c$ in Algorithm 1. <br/> 
model_solving_w: the function that computes the optimal solution for Pro- $S_c$. <br/> 
model_solving_cont: the function that computes the optimal solution for Pro- $S_c$ with given choices of $\alpha$ and $\beta$. <br/>
