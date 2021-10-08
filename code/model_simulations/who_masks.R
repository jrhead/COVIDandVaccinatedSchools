#####################################################################################  
#### Script to evaluate infections, illnesses, hopsitalization, and deaths for:######
#       1) Everyone masks                                                           # 
#       2) Only people who did not receive the vaccine masks                        # 
#        -- across levels of vaccine effectiveness                                  # 
#####################################################################################  

########################## 
#Last modified: 9/30/2021#
########################## 

#############################################################
######  1. DEFINE MODEL PARAMETERS            ###############
#############################################################

# MODEL PARAMETERS VARIED IN ANALYSES 
alpha = 0.5 #ratio of asymptomatic to symptomatic transmission; set to 0.5 in article

susceptRatio = 1 #ratio of the susceptibility of children < (susceptAgeSplit) to SARS-CoV-2 vs. adults; varied between 0.5 and 1
susceptAgeSplit = 10 #age where susceptibility differs; here, agents < 10 assumed half as susceptible 

#Vaccine effectiveness for various disease endpoints
VE_any = 0.77 #Higdon
VE_symp = 0.85 #Lopez Bernal, NEJM
VE_sev = 0.93 #CDC MMWR

vacc_prop = 0.50 #proportion of the population vaccinated -- Make 0 for no vacc situation

R0_init = 5*0.844 + 2.5*(1-0.844) #basic reproduction number, weighted by which is Delta

##Define waiting times -- parameters in a Weibull distribution
#E --> I; I --> R or D; I --> H; H --> R or D
shapeEI = 4; scaleEI = 6
shapeIR = 7; scaleIR = 8
shapeIH = 7; scaleIH = 7.8
shapeHRD = 13; scaleHRD = 15

# COMMUNITY TRANSMISSION PARAMETERS VARIED IN ANALYSES
propMidEssentialTot = 0.85 #Proportion of adults still working; 
#set to 85% of in-person work for high transmission scenario; 

CommContactIncr = 2 #indicates increases in socializing over summer; 
#set to 2 (twice as many community contacts as observed during survey) 

N = 16000 # number of agents
NWork <- round((N/2)/20,0) #number of work places (avg size = 20)

#################################################################################
######   2. Load functions, and data for community contact matrices #############
#################################################################################
### LOAD R SCRIPTS AND FUNCTIONS ###
library(here)
here() #check for correct working directory

source(".//code//functions//contact_matrix_VM.R") #make contact matrices
source(".//code//functions//model_functions_VM.R")

#Load starting values -- list of fates, synthetic population, and starting states
load(".//code//starting_states//fate_1000reps.RData") #fates -- fate_list
load(".//code//starting_states//starting_states_1000reps.RData") #starting conditions -- start_state_list
load(".//code//starting_states//synth_pops_1000reps.RData") #synthetic pop -- synth_pop_list

#Load the community matrices for the third wave of the survey
load(".//data//survey_contact_data_Feb//FebComMatrix_med.RData")
load(".//data//survey_contact_data_Feb//workFeb.RData")

#Calculate the observed community contact matrix in February
Obs.comm_feb <- FebComMatrix_med - array_work_feb #subtract work because work gets added back in in the work matrix

#Load the polymod data for the early epidemic period
polymod_dat <- readRDS(".//data//POLYMOD_contact_data//polymod_dat.RData")
matrix <- as.matrix(polymod_dat$matrix)

######################################################################
######   3. Run simulations for each reopening scenario ##############
######################################################################

##set up parallel processing
library(doParallel)
library(foreach)

nCores <- 5
#nCores <- as.numeric(Sys.getenv('SLURM_CPUS_ON_NODE'))
registerDoParallel(nCores)
set.seed(111)

#Run all situations 1000 times, in parallel
outcomes <- foreach(reps = 1:1000, .packages = "dplyr") %dopar% {
  
  ################## 0. PREP BY CREATING SYNTHETIC POPULATION AND THEIR FATES ################################
  #Generate synetic population
  
  synth_pop_df <- synth_pop_list[[reps]]
  synth_pop_df$Age <- as.numeric(as.character(synth_pop_df$Age))
  
  #Define number in each age group
  Age_Vect <- c(sum(synth_pop_df$AgeCat == "Under5"),
                sum(synth_pop_df$AgeCat == "5-10"),
                sum(synth_pop_df$AgeCat == "11-13"),
                sum(synth_pop_df$AgeCat == "14-17"),
                sum(synth_pop_df$AgeCat == "18-64"),
                sum(synth_pop_df$AgeCat == "65+"))
  
  # For constant contact rate between age groups
  AgeContRates <- matrix(c(
    matrix[1,]/Age_Vect[1],
    matrix[2,]/Age_Vect[2],
    matrix[3,]/Age_Vect[3],
    matrix[4,]/Age_Vect[4],
    matrix[5,]/Age_Vect[5],
    matrix[6,]/Age_Vect[6]), ncol = 6)
  
  rownames(AgeContRates) <- c("Under5", "5-10", "11-13", "14-17", "18-64", "65+")
  
  #define new age categories

  Age_Vect2 <- c(sum(synth_pop_df$AgeCat2 == 1),
                 sum(synth_pop_df$AgeCat2 == 2),
                 sum(synth_pop_df$AgeCat2 == 3),
                 sum(synth_pop_df$AgeCat2 == 4),
                 sum(synth_pop_df$AgeCat2 == 5),
                 sum(synth_pop_df$AgeCat2 == 6))

  
  #Define original contact matrix under a no intervention situation
  contacts_orig_list <- Rcontact_matrix(synth_pop_df, 
                                        HHcontact = 1, 
                                        SCHcontact = .2/7,
                                        WORKcontact = 5/7 ,
                                        GRADEcontact = 1/7-.2/7, 
                                        CLASScontact = 5/7- 1/7-.2/7, 
                                        AgeContRates,
                                        propEssential=0.28,
                                        propMidEssential = propMidEssentialTot-0.28,
                                        NWork)
  
  #calculate eigenvalue of contact matrix (used in calculation of beta below)
  #expressed as the weighted mean of total contacts (mean # of contacts from someone who has just been contacted)
  numCont <- rowSums(contacts_orig_list[[1]]) + 
    rowSums(contacts_orig_list[[2]]) + 
    rowSums(contacts_orig_list[[3]]) + 
    rowSums(contacts_orig_list[[4]]) + 
    rowSums(contacts_orig_list[[5]]) + 
    rowSums(contacts_orig_list[[6]])
  MAXEIG <- weighted.mean(numCont, w = numCont)
 
  state0 <- start_state_list[[reps]]
  
  #Re-allocate the vaccination within the start states to fit with initial parameters for vaccination...
  # Create the vaccination vector -- who gets vaccinated and for whom is vaccination effective?
  got_v <- rep(0,N)
  got_v[synth_pop_df$Age < 12] <- 0
  got_v[synth_pop_df$Age >= 12] <- rbinom(length(got_v[synth_pop_df$Age >= 12]), size = 1, prob = vacc_prop)
  
  state0[state0 == "V"] <- "S"
  state0[state0 == "S" & got_v == 1] <- "V"
  
  VE = got_v*VE_any
  
  #Define conditional probabilities for clinical outcomes 
  #rho: prob infection is clinical
  rho = rep(0.69, N) #prob case is clinical for under 18 is 21%, 69% for those older...
  rho[synth_pop_df$AgeCat == "Under5" | synth_pop_df$AgeCat == "5-10" |
        synth_pop_df$AgeCat == "11-13"| synth_pop_df$AgeCat == "14-17"] <- 0.21
  
  #update P(clin|age,vac) for thos who are vaccinated
  rho[got_v == 1] <- rho[got_v == 1] * ((1-VE_symp)/(1-VE_any))
  
  #Define conditional probability of hospitalization
  h_rate <- rep(NA, N) #hospitalization rates
  #SOURCE: https://www.thelancet.com/journals/laninf/article/PIIS1473-3099(20)30243-7/fulltext
  h_rate[synth_pop_df$Age <10] <- 0.00001
  h_rate[synth_pop_df$Age >= 10 & synth_pop_df$Age < 20] <- 0.000408
  h_rate[synth_pop_df$Age >= 20 & synth_pop_df$Age < 30] <- 0.0104
  h_rate[synth_pop_df$Age >= 30 & synth_pop_df$Age < 40] <- 0.0343
  h_rate[synth_pop_df$Age >= 40 & synth_pop_df$Age < 50] <- 0.0425
  h_rate[synth_pop_df$Age >= 50 & synth_pop_df$Age < 60] <- 0.0816
  h_rate[synth_pop_df$Age >= 60 & synth_pop_df$Age < 70] <- 0.118
  h_rate[synth_pop_df$Age >= 70 & synth_pop_df$Age < 80] <- 0.166
  h_rate[synth_pop_df$Age >= 80] <- 0.184
  
  #update P(clin|age,vac) for thos who are vaccinated
  h_rate[got_v == 1] <- h_rate[got_v == 1] * ((1-VE_sev)/(1-VE_symp))
  
  #define conditional probability of death
  mu = rep(NA, N) #death rates, among those hospitalized
  #SOURCE: https://www.bmj.com/content/369/bmj.m1923
  mu[synth_pop_df$Age <10] <- 0.017#0.5*(0.024 + 0.017)
  mu[synth_pop_df$Age >= 10 & synth_pop_df$Age < 20] <- 0.017 #0.5*(0.024 + 0.017)
  mu[synth_pop_df$Age >= 20 & synth_pop_df$Age < 30] <- 0.026 #0.5*(0.036 + 0.026)
  mu[synth_pop_df$Age >= 30 & synth_pop_df$Age < 40] <- 0.036 #0.5*(0.059 + 0.036)
  mu[synth_pop_df$Age >= 40 & synth_pop_df$Age < 50] <- 0.062 #0.5*(0.095 + 0.062)
  mu[synth_pop_df$Age >= 50 & synth_pop_df$Age < 60] <- 0.099 #0.5*(0.143 + 0.099)
  mu[synth_pop_df$Age >= 60 & synth_pop_df$Age < 70] <- 0.151 #0.5*(0.221 + 0.151)
  mu[synth_pop_df$Age >= 70 & synth_pop_df$Age < 80] <- 0.238 #0.5*(0.363 + 0.238)
  mu[synth_pop_df$Age >= 80] <- 0.365#0.5*(0.365+0.538)
  
  ##update fates for each agent
  fate_pre_vacc <- fate_list[[reps]]
  
  ##determine fates for each agent
  wait_times <- update_fates(N, state0,
                             fate_pre_vacc,
                             synth_pop_df,
                             R0_init, MAXEIG, 
                             alpha, susceptRatio, susceptAgeSplit,
                             shapeEI, scaleEI, 
                             shapeIR, scaleIR, 
                             shapeIH, scaleIH, 
                             shapeHRD, scaleHRD, 
                             rho, 
                             h_rate, 
                             mu)
  

  ################## Reopening (Aug 9 - Dec 17) ##################
  ################## Sim1: Schools Fully Closed ############################
  
  #Update contact matrices
  CLOSED_contact_list <- contacts_orig_list
  
  EWrk.cont.mid <- contacts_orig_list[[11]]
  MidEssentialWorkMat <- matrix(0, ncol = N, nrow = N)
  for (i in 1:dim(EWrk.cont.mid)[1]){
    MidEssentialWorkMat[EWrk.cont.mid[i,1],EWrk.cont.mid[i,2]] <- 5/7
    MidEssentialWorkMat[EWrk.cont.mid[i,2],EWrk.cont.mid[i,1]] <- 5/7
  }
  MidEssentialWorkMatsp <- as(MidEssentialWorkMat, "sparseMatrix")
  rm(MidEssentialWorkMat)
  
  CLOSED_contact_list[[2]] <- MidEssentialWorkMatsp #Work
  CLOSED_contact_list[[3]] <- 0*CLOSED_contact_list[[3]] #School
  CLOSED_contact_list[[4]] <- 0*CLOSED_contact_list[[4]] #Grade
  CLOSED_contact_list[[5]] <- 0*CLOSED_contact_list[[5]] #Class
  CLOSED_contact_list[[6]] <- gen_COMmat(Obs.comm_feb, Age_Vect2, synth_pop_df) #community, using survey data
  
  #Call SEIR function
  StatesPhase_CLOSED <- SEIRQ_func(synth_pop_df, 
                                    start_day = 1, #August 9
                                    end_day = 128, #Dec 17
                                    contacts_orig_list = CLOSED_contact_list, 
                                    alpha,
                                    beta = wait_times[["beta"]],
                                    time_state1 = wait_times[["time_state1"]], time_next_state1 = wait_times[["time_next_state1"]],
                                    time_state2 = wait_times[["time_state2"]], time_next_state2 = wait_times[["time_next_state2"]],
                                    time_state3 = wait_times[["time_state3"]], time_next_state3 = wait_times[["time_next_state3"]], 
                                    fate = wait_times[["fate_updated"]],
                                    state_full = state0,
                                    state = state0,
                                    VE = VE,
                                    propComplyCase = 0.80*0.6, 
                                    propComplyHH = 0.25,
                                    redContacts = 0.75,
                                    mask_ef_teacher = 0,
                                    mask_ef_elem = 0,
                                    mask_ef_middle = 0,
                                    mask_ef_HS = 0)
  
  state_full.CLOSED <- StatesPhase_CLOSED[["state_full"]]
  
  #extract summary indicators
  outcomeCLOSED <- modelled_outcomes(state_full.CLOSED, synth_pop_df)
  
  rm(state_full.CLOSED);
  rm(CLOSED_contact_list)
  gc()
  #####################################

  ################## Sim1: Mask for all persons ############################
  OPEN_contact_list <- contacts_orig_list
  
  OPEN_contact_list[[2]] <- MidEssentialWorkMatsp #Work
  OPEN_contact_list[[6]] <- gen_COMmat(Obs.comm_feb, Age_Vect2, synth_pop_df) #Community, using survey data
  
  #call SEIR funtion -- no testing
  StatesPhase_MASK <- SEIRQ_func(synth_pop_df, 
                                  start_day = 1, #Aug 15
                                  end_day = 128, #Dec 20
                                  contacts_orig_list = OPEN_contact_list, 
                                  alpha,
                                  beta = wait_times[["beta"]],
                                  time_state1 = wait_times[["time_state1"]], time_next_state1 = wait_times[["time_next_state1"]],
                                  time_state2 = wait_times[["time_state2"]], time_next_state2 = wait_times[["time_next_state2"]],
                                  time_state3 = wait_times[["time_state3"]], time_next_state3 = wait_times[["time_next_state3"]], 
                                  fate = wait_times[["fate_updated"]],
                                  state_full = state0,
                                  state = state0,
                                  VE = VE,
                                  propComplyCase = 0.80*0.6, 
                                  propComplyHH = 0.25,
                                  redContacts = 0.75,
                                  mask_ef_teacher = 0.5, #wearing masks
                                  mask_ef_elem = 0.15,
                                  mask_ef_middle = 0.25,
                                  mask_ef_HS = 0.35)
  
  
  state_full.MASK <- StatesPhase_MASK[["state_full"]]
  
  #extract summary indicators
  outcomeMASK <- modelled_outcomes(state_full.MASK, synth_pop_df)
  
  rm(state_full.MASK);
  gc()
  

  ################## Sim1: Mask for unvaccinated persons ############################
  
  #call SEIR funtion -- no testing
  StatesPhase_MASK_UV <- SEIRQ_func(synth_pop_df, 
                                  start_day = 1, #Aug 15
                                  end_day = 128, #Dec 20
                                  contacts_orig_list = OPEN_contact_list, 
                                  alpha,
                                  beta = wait_times[["beta"]],
                                  time_state1 = wait_times[["time_state1"]], time_next_state1 = wait_times[["time_next_state1"]],
                                  time_state2 = wait_times[["time_state2"]], time_next_state2 = wait_times[["time_next_state2"]],
                                  time_state3 = wait_times[["time_state3"]], time_next_state3 = wait_times[["time_next_state3"]], 
                                  fate = wait_times[["fate_updated"]],
                                  state_full = state0,
                                  state = state0,
                                  VE = VE,
                                  propComplyCase = 0.80*0.6, 
                                  propComplyHH = 0.25,
                                  redContacts = 0.75,
                                  mask_ef_teacher = 0.5, #wearing masks
                                  mask_ef_elem = 0.15,
                                  mask_ef_middle = 0.25,
                                  mask_ef_HS = 0.35,
                                  who_mask = as.numeric(got_v == 0)) #only people who did not get the vaccine
  
  
  state_full.MASKUV <- StatesPhase_MASK_UV[["state_full"]]
  
  #extract summary indicators
  outcomeMASKUV <- modelled_outcomes(state_full.MASKUV, synth_pop_df)
  
  rm(state_full.MASKUV);
  gc()
  
  #########################################################################
  # Combine output
  
  outcomes <- list("CLOSED" = outcomeCLOSED, 
                  "MASK" = outcomeMASK, 
                  "MASKUV" = outcomeMASKUV,
                  "fate" = wait_times[["fate_updated"]],
                  "got_v" = got_v)
  outcomes
}

#Save model simulations
save(outcomes, file = "VaccReopening_WhoMasks_45VE.RData")