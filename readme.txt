# read me

1_maddfitz_true_simulate.R: a script used to simulate data according to darwin's scenario, unreplicated bursts scenario, modified darwin's scenario (inside), modified Darwin's scenario (outside), modified Darwin's scenario (both). This script also fits all models to our set and outputs a file: maddfitz_res.Rsave. This output provides the bulk of the results presented in the main paper. 

2_maddfitz_true_analyze.R: a script which analyzes maddfitz_res.Rsave, producing both numerical results of the text, Figure 5, Table 1, and several supplemental figures.

3_likelihood_contour.R: a script that analyzes maddfitz_res.Rsave to produce a likelihood surface for particular parameter estimates and demonstrates that the rates to and from intermediate states are unidentifiable. 

4_maddfitz_simmed_simulate.R: a script used to simulate data according to three simplified models (simplified correlated, simplified independent, and simplified hidden Markov independent). This script also fits these models to each of their datasets generating the result file: maddfitz_simmed_res.Rsave

5_maddfitz_simmed_analyze.R: a script which analyzes maddfitz_simmed_res.Rsave and provides evidence that the generating model is chosen when data is simulated under the simplified models.