#### PhenoCam_data_download_and_modeling
This R code uses the phenocamr and phenor packages to download and combine data from eight different PhenoCam sites in the Northeastern US that can be used to train regional phenology models.

The spring and fall phenology models were calibrated using the phenoR package (Hufkens et al. 2018). We first parameterized models for the beginning and end of the season based on dates extracted from PhenoCam: 10% of the seasonal amplitude in the spring 10% of the seasonal amplitude in the fall. We chose the threshold of 10% because it is as close to the “true” onset period as we can detect (Richardson et al. 2018). Then, to model the duration of the transition period, we trained models to predict the end of the spring transition, and the beginning of the fall transition by holding all other parameter values constant except for the critical value of accumulating temperature marking the end of the transition period. By doing so, we are effectively modeling 4 dates for each year – the initiation and completion of both spring onset and fall senescence for each model. It is important to model the beginning and end of the transition separately because the rate of leaf development or senescence can vary from year-to-year (Klosterman et al. 2018). We optimized models using a generalized simulated annealing algorithm (Xiang et al. 2013) in the phenoR package used to identify the global minima of a combination of multiple different parameter values. For each model, we ran 25 chains of 100,000 iterations and identified the parameter sets with the lowest root mean square error (RMSE) as the parameters that best fit the data. Then the best phenology models were selected based on the lowest Akaike information criterion (AIC) values. 



##### Incorportating new phenology models into PnET

The the C++ code included is an altered phenology subroutine from the PnET-CN biogeochemical model. It uses the optimized parameter values from the PhenoCam data to update the model's ability to predict spring and fall transitions. For the most part, we were able to write the new phenology model to run with minimal changes in the rest of the code. However, the accumulating warming and chilling temperature variables must first be declared in the pnet_model.h subroutine (code can be found below).

To incorporate the new phenology in PnET, copy and paste the included phenology subroutine (pnet_Phenology.cpp) into the PnET-CN project, To select the spring and fall models, change the numbers referring to spring_mod_select or fall_mod_select

Before the models will run, you first need to declare the shared variables that accumulate overtime in pnet_model.h., and set the initial values in pnet_init_vars

1.) Add these two lines to pnet_model.h:

double CDDTot;		// total chilling degree days // AFT

double Tsum;		// secondary GDD calculation // AFT

2.) Add these lines to set the initial values of the new shared variables to 0 in pnet_init_vars

share->GDDTot=0; //AFT

share->CDDTot=0; //AFT

share->Tsum=0; //AFT

