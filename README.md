# School closures reduced social mixing of children during COVID-19 with implications for transmission risk and school reopening policies
Analysis of school closure and reopening policies in COVID-19 outcomes in the Bay Area (California)
# Documentation

## Introduction

This document explains the process for reproducing the data and analysis described in the paper entitled "Model-based assessment of SARS-CoV-2 Delta variant transmission dynamics within partially vaccinated K-12 school populations". The objectives of this study were to estimate increase in the total number of symptomatic infections among students and teachers/staff between in-school and remote instruction over a 128-day semester when 1) various NPIs (mask use, cohorts, and weekly testing of students/teachers) were implemented in schools; 2) across various community-wide vaccination coverages (50%, 60%, 70%), and 3) student ($\geq$ 12 years) and teacher/staff vaccination coverages (50% - 95%); and 4) for universal masking vs. masking only among the unvaccinated.

## Downloading code and data

Code and data files are available from from this repository at: http://www.github.com/jrhead/COVIDandVaccinatedSchools. The code is modified from a previously published model (doi: https://doi.org/10.1098/rsif.2020.0970), with code available from: http://www.github.com/jrhead/COVIDandSchools.

A description of each folder and file are below:

### data

The **data** folder contains social contact data obtained from both the Bay Area social contact survey, as well as the POLYMOD social contact surveys from the UK. The folder is sub-divided into two directories, with the following contents:

1. **POLYMOD_contact_data**: directory which contacts code and age-structured contact matrices during pre-pandemic times

    * **polymod_community_matrix_generation.R** contains the source code to download the POLYMOD age-structured contact matrices by location, using the socialmixr package.

    * **polymod_dat.RData** contains an age-structured contact matrix representing non-household contacts summed across all locations.

    * **polymod_dat.CF2_4.RData** contains an age-structured contact matrix representing non-household, non-work, non-school, and non-transportation contacts. It is combined with the work, school, and transportation matrices from the Bay Area synthetic population and the Bay Area survey to generate the community contact matrix for counterfactual scenarios 2 and 4 (workplaces open).

    * **polymod_dat.CF3_5.RData** contains an age-structured contact matrix representing transportation contacts. It is combined with the work and school, matrices from the Bay Area synthetic population and the non-transportation, non-household contacts from the Bay Area survey to generate the community contact matrix for counterfactual scenarios 3 and 5 (socializing permitted).

    * **polymod_dat.CF6.RData** contains an age-structured contact matrix representing non-household, non-work, and non-school contacts. It is combined with the work and school matrices from the Bay Area synthetic population to generate the community contact matrix for counterfactual scenario 6 (workplaces & schools open; socializing permitted)

2. **survey_contact_data**: directory which contains the age-structred contact matrices by location during the pandemic in California Bay Area, May 4 - June 1, 2020.

    * **community-matrix_060820.RData**: contains an age-structured contact matrix representing the mean number of non-household contacts that occurred in the past 24 hours in the Bay Area at all locations.

    * **childcare_060820.RData**: contains an age-structured contact matrix representing the mean number of non-household contacts that occurred in the past 24 hours in the Bay Area at a child care setting.

    * **transit_060820.RData**: contains an age-structured contact matrix representing the mean number of non-household contacts that occurred in the past 24 hours in the Bay Area on public transport.

    * **work_060820.RData**: contains an age-structured contact matrix representing the mean number of non-household contacts that occurred in the past 24 hours in the Bay Area at the respondent's place of work.
    
3. **survey_contact_data_Feb**: directory which contains the age-structred contact matrices by location during the pandemic in California Bay Area, February 8 - April 1, 2021.

    * **FebComMatrix_med.RData**: contains an age-structured contact matrix representing the mean number of non-household contacts that occurred in the past 24 hours in the Bay Area at all locations.

    * **childcareFeb.RData**: contains an age-structured contact matrix representing the mean number of non-household contacts that occurred in the past 24 hours in the Bay Area at a child care setting.

    * **transitFeb.RData**: contains an age-structured contact matrix representing the mean number of non-household contacts that occurred in the past 24 hours in the Bay Area on public transport.

    * **workFeb.RData**: contains an age-structured contact matrix representing the mean number of non-household contacts that occurred in the past 24 hours in the Bay Area at the respondent's place of work.

### code

The **code** folder contains all the R scripts necessary to run all instances of the SEIR model. It is divided into three folders; a) **functions** contains the functions called in the main script to implement the SEIR;  b) **starting_states** contains .RData files containing the 1,000 various observations of the synthetic population (N = 16,000), the fates of each individual in the synthetic population (e.g., whether the individual will be asymptomatic, symptomatic, hospitalized, etc.), and the initial states of the 16,000 individuals in the synthetic population; c) **model_simulations** contains the main scripts that implement the SEIR These scripts are modifiable based on model parameters.

1. **functions** directory: There are two files in this main directory, which contain functions that are called in the main scripts (see information under #3, model_simulations). 

    * **contact_matrixVM.R** contains a function to return an NxN sparse matrix indicating contacts between pairs in the synthetic population.

    * **model_functionVM.R** contains a series of functions needed in the agent based model for simulating reopening strategies. These functions include:
    
        1. *Fate_func*:  updates fates of individuals, and assigns waiting times
        2. *gen_COMmat*: generates and NxN community matrix for the full population
        3. *SEIRQ_func*: runs the SEIR model, including an intervention where symptomatic individuals isolate
        4. *SEIRtest_func*: adds to the SEIRQ_func to permit testing on a particular interval, following by quarantine/isolation
        5. *modelled_outcomes*: counts/summarizes important outcomes (e.g. cumulative infections, daily hospitalizations)


2. **starting_states** directory: There are three .RData files in this main directory. Each contains a list of 1,000 realizations of a particular outcome for the synthetic population. The lists are: 

    * **synth_pops_1000reps.RData** This .RData file is a list of 1,000 realizations of a dataframe holding the a synthetic population of size N=16,000, where each agent is assigned an age and membership in a household, community, and -- depending on age -- a workplace or school/grade/class.

    * **fate_1000reps.RData** This .RData file is a list of 1,000 realizations of a vector holding one of four fates for each of the 16,000 agents in the synthetic population. These fates determine the pathway that the agent would follow should they be infected with SARS-CoV-2. The fate is drawn from independent realizations of a Bernoulli distribution based on age-conditional probabilities (see Table S1; Figure 1B in manuscript). The fates are: A: asymptomatic infection; C: clinical infection; H: infection resulting in survival and hospitalization; D: infection resulting in death.

    * **starting_sates_100reps.RData** This .RData file is a list of 1,000 realizations of a vector holding the starting state for each of the 16,000 agents in the synthetic population. These states include: V: successfully vaccinated; S: susceptible; R: recovered; A: infected, asymptomatic; C: infected, symptomatic; H1/H2: infected, hospitalized; D1/D2: infected, died.


3. **model_simulations** directory: There are three files in this main directory. Each can be updated and executed according to the desired parameterization of the model. Each script can be directly modified to change the proportion of the wider community and the school community who is vaccinated, the relative susceptibility of children to SARS-CoV-2 compared to adults, the vaccine effectiveness, and other natural history parameters..

    * **interventions.R** The main R script which implements the SEIR under a school closure scenario, a school open scenario, and a school open with NPIs scenario: masks; masks + cohorts; masks + weekly testing. 

    * **who_masks.R** The main R script which implements the SEIR under a school open scenario where everyone wears masks and a school open scenario where only the unvaccinated wear masks. 

    * **within_school_coverages.R** The main R script which implements the SEIR under a school open scenarios without NPIs but increasing levels of vaccination coverage among the eligible student and teacher population. 


## Running the scripts

### Examine the effect of various non-pharmaceutical interventions for various vaccination coverages.

To examine the added infections when schools are open compared to when schools are remote, across various NPIs (masks, masks + cohorts, masks + testing) and vaccine coverages, options at the top of the R script, **interventions.R**, can be modified. Some are shown here, including the ratio of asymptomatic to symptomatic transmission ($\alpha$), the susceptibility ratio of children less than a specified age to children above that age, the vaccine effectivness, the vaccine coverage, and $R_0$:

```
# MODEL PARAMETERS VARIED IN ANALYSES 
alpha = 0.5 #ratio of asymptomatic to symptomatic transmission; set to 0.5 in article

susceptRatio = 1 #ratio of the susceptibility of children < (susceptAgeSplit) to SARS-CoV-2 vs. adults; varied between 0.5 and 1
susceptAgeSplit = 10 #age where susceptibility differs; here, agents < 10 assumed half as susceptible 

vacc_eff = 0.85 #vaccination effectiveness -- set to 0.85 in article
vacc_prop = 0.50 #proportion of the population vaccinated -- Make 0 for no vacc situation

R0_init = 5*0.844 + 2.5*(1-0.844) #basic reproduction number, weighted by which is Delta
```

Other parameters may be varied throughout, including other population level parameters.

To load needed functions into the main R script, follow the below code, making sure the correct working directory is specified.

```
library(here)
here() #check for correct working directory

source("..//functions//contact_matrix_VM.R") #make contact matrices
source("..//functions//model_functions_VM.R")

#Load starting values -- list of fates, synthetic population, and starting states
load("..//starting_states//fate_1000reps.RData") #fates -- fate_list
load("..//starting_states//starting_states_1000reps.RData") #starting conditions -- start_state_list
load("..//starting_states//synth_pops_1000reps.RData") #synthetic pop -- synth_pop_list

#Load the community matrices for the third wave of the survey
load("..//..//data//survey_contact_data_Feb//FebComMatrix_med.RData")
load("..//..//data//survey_contact_data_Feb//workFeb.RData")
```

This script was implemented in parallel on a computer with multiple cores. To change the code to fit with different computer specifications, see these lines of code:

```
##set up parallel processing
nCores <- 5
registerDoParallel(nCores)
```

To update the name of the simulations that will be saved in an .RData file, update the name in the final line in the script:

```
save(outcomes, file = "Int_50cov_85VE_HS.RData")
```

### Examine the effect of universal masking compared to only masking unvaccinated individuals

To examine the added infections when everyone masks compared to when only the unvaccinated mask, options at the top of the R script, **who_masks.R**, can be modified, as shown previously.

The key parameter that tells the SEIRQ_func to only apply masks among the unvaccinated is the who_mask argument in the `SEIRQ_func` function. The default value for who_mask is 1 for everyone, which tells the model to apply the mask effectiveness shown (e.g., `mask_ef_teacher`) to everyone in the population group. By updating the who_mask argument in the `SEIRQ_func` function to `as.numeric(got_v == 0)`, it tells the function to only apply masks to those who did not recieve the vaccine. The vector `got_v` here includes those who got the vaccine, but were not succesfully immunized. 

Mask effectiveness can also be modified in these arguments as well. `mask_ef_teacher` represents the effectiveness of masks among teachers, `mask_ef_elem` represents the effectiveness among elementary schools students, and so on.

```
StatesPhase_MASK_UV <- SEIRQ_func(synth_pop_df, 
                                  start_day = 1, #Aug 15
                                  end_day = 128, #Dec 20
                                  contacts_orig_list = OPEN_contact_list, 
                                  alpha,
                                  beta = wait_times[["beta"]],
                                  time_state1 = wait_times[["time_state1"]], 
                                  time_next_state1 = wait_times[["time_next_state1"]],
                                  time_state2 = wait_times[["time_state2"]], 
                                  time_next_state2 = wait_times[["time_next_state2"]],
                                  time_state3 = wait_times[["time_state3"]], 
                                  time_next_state3 = wait_times[["time_next_state3"]], 
                                  fate = fate_i,
                                  state_full = state0,
                                  state = state0,
                                  propComplyCase = 0.80*0.6, 
                                  propComplyHH = 0.25,
                                  redContacts = 0.75,
                                  mask_ef_teacher = 0.5, #wearing masks
                                  mask_ef_elem = 0.15,
                                  mask_ef_middle = 0.25,
                                  mask_ef_HS = 0.35,
                                  who_mask = as.numeric(got_v == 0))
```
Other parameters may be varied throughout, including other population level parameters.

### Examine the effect of increasing vaccination coverage within school population, in absence of others NPIs

To examine the added infections when everyone masks compared to when only the unvaccinated mask, options at the top of the R script, **within_school_coverages.R**, can be modified. These are shown here:

```
# MODEL PARAMETERS VARIED IN ANALYSES
teach_vacc_prop = c(0.8, 0.9, 0.95) #proportion of the teachers vaccinated 
stud_vacc_prop = c(0.8, 0.9, 0.95) #proportion of the students >=12  vaccinated 

```
Other parameters may be varied throughout, including other population level parameters.
