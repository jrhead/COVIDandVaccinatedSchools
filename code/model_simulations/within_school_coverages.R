#####################################################################################  
#### Script to evaluate infections, illnesses, hopsitalization, and deaths for:######
#       1) No NPIs                                                                  # 
#       2) 70% vaccine coverage of the eligible non-school community                # 
#       3) 70 - 95% coverage of the eligible school community                       # 
#####################################################################################  

########################## 
#Last modified: 8/22/2022#
########################## 

#############################################################
######  1. DEFINE MODEL PARAMETERS            ###############
#############################################################

# MODEL PARAMETERS VARIED IN ANALYSES 
alpha = 0.5 #ratio of asymptomatic to symptomatic transmission; set to 0.5 in article

susceptRatio = 1 #ratio of the susceptibility of children < (susceptAgeSplit) to SARS-CoV-2 vs. adults; varied between 0.5 and 1
susceptAgeSplit = 10 #age where susceptibility differs; here, agents < 10 assumed half as susceptible 

vacc_eff = 0.85 #vaccination effectiveness 
vacc_prop = 0.50 #proportion of the community vaccinated -- Make 0 for no vacc situation
teach_vacc_prop = c(0.8, 0.9, 0.95) #proportion of teachers vaccinated 
stud_vacc_prop = c(0.8, 0.9, 0.95) #proportion of students >12 vaccinated 

R0_init = 5*0.844 + 2.5*(1-0.844) #basic reproduction number, weighted by which is Delta

##Define waiting times -- parameters in a Weibull distribution
#E --> I; I --> R or D; I --> H; H --> R or D
shapeEI = 4; scaleEI = 6
shapeIR = 7; scaleIR = 14
shapeIH = 7; scaleIH = 11
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
  synth_pop_df$Age <- as.numeric(as.character(synth_pop_df$Age))
  
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
 
  #Define conditional probabilities for clinical outcomes 
  #rho: prob infection is clinical
  rho = rep(0.69, N) #prob case is clinical for under 18 is 21%, 69% for those older...
  rho[synth_pop_df$AgeCat == "Under5" | synth_pop_df$AgeCat == "5-10" |
        synth_pop_df$AgeCat == "11-13"| synth_pop_df$AgeCat == "14-17"] <- 0.21
  
  synth_pop_df$Age <- as.numeric(as.character(synth_pop_df$Age))
  
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
  fate_i <- fate_list[[reps]]
  
  state0 <- start_state_list[[reps]]
  
  #Re-allocate the vaccination within the start states to fit with initial parameters for vaccination...
  # Create the vaccination vector -- who gets vaccinated and for whom is vaccination effective?

  CommV <- TeachV <- SchoolElg <- rep(0,N)
  CommV[synth_pop_df$Age > 18 & is.na(synth_pop_df$School)] <- 1
  TeachV[synth_pop_df$Age > 18 & !is.na(synth_pop_df$School)] <- 1
  SchoolElg[synth_pop_df$Age >= 12 & synth_pop_df$Age <= 18 & !is.na(synth_pop_df$School)] <- 1
  
  ## First get some consistent things -- 1) community vaccination
  V <- rep(0, N)
  V[CommV == 1] <- rbinom(length(V[CommV == 1]), size = 1, prob = vacc_eff*vacc_prop)
  
  state0[state0 == "V"] <- "S"
  state0[state0 == "S" & V == 1] <- "V"
  
  #initialize storage lists for open and closed
  outcomeCLOSEDLists <- outcomeOPENLists <- list()
  
  for (vi in 1:length(stud_vacc_prop)){
    
    #reset the vaccination and initial state vectors
    Vs <- V
    state0v<- state0
    
    #Update the vaccinated teachers and students
    Vs[TeachV == 1] <- rbinom(length(Vs[TeachV == 1]), size = 1, prob = vacc_eff*teach_vacc_prop[vi])
    Vs[SchoolElg == 1] <- rbinom(length(Vs[SchoolElg == 1]), size = 1, prob = vacc_eff*stud_vacc_prop[vi])
    
    #Overwrite some states as being vaccinated
    state0v[Vs==1] <- "V"
    
    ##determine fates for each agent
    wait_times <- update_fates(N, fate_i, state0v, 
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
     CLOSED_contact_list[[6]] <- gen_COMmat(Obs.comm_feb, Age_Vect2, synth_pop_df) #Community, using survey data
    
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
                                      fate = fate_i,
                                      state_full = state0v,
                                      state = state0v,
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
  
    ################## Sim2: Opening schools fully!!! ############################
    OPEN_contact_list <- contacts_orig_list
    
    OPEN_contact_list[[2]] <- MidEssentialWorkMatsp
    #OPEN_contact_list[[6]] <- gen_COMmat(CommContactIncr*Obs.comm) #Community, EXACTLY as observed, with work, daycare, public trans included
    OPEN_contact_list[[6]] <- gen_COMmat(Obs.comm_feb, Age_Vect2, synth_pop_df)
    
    #Call SEIR function
    StatesPhase_OPEN <- SEIRQ_func(synth_pop_df, 
                                    start_day = 1, #Aug 9
                                    end_day = 128, #Dec 17 
                                    contacts_orig_list = OPEN_contact_list, 
                                    alpha,
                                    beta = wait_times[["beta"]],
                                    time_state1 = wait_times[["time_state1"]], time_next_state1 = wait_times[["time_next_state1"]],
                                    time_state2 = wait_times[["time_state2"]], time_next_state2 = wait_times[["time_next_state2"]],
                                    time_state3 = wait_times[["time_state3"]], time_next_state3 = wait_times[["time_next_state3"]], 
                                    fate = fate_i,
                                    state_full = state0v,
                                    state = state0v,
                                    propComplyCase = 0.80*0.6, 
                                    propComplyHH = 0.25,
                                    redContacts = 0.75,
                                    mask_ef_teacher = 0,
                                    mask_ef_elem = 0,
                                    mask_ef_middle = 0,
                                    mask_ef_HS = 0)
    
    state_full.OPEN <- StatesPhase_OPEN[["state_full"]]
    
    #extract summary indicators
    outcomeOPEN <- modelled_outcomes(state_full.OPEN, synth_pop_df)
    
    rm(state_full.OPEN);
    
    gc()
    
    ##########################################
    
    outcomeCLOSEDLists[[vi]] <- append(outcomeCLOSEDLists, outcomeCLOSED)
    outcomeOPENLists[[vi]] <- append(outcomeOPENLists, outcomeOPEN)
  }

  rm(OPEN_contact_list)
  rm(CLOSED_contact_list)

  
  #########################################################################
  # Combine output
  
  outcomes <- list("CLOSED" = outcomeCLOSEDLists, 
                   "OPEN" = outcomeOPENLists)
  outcomes
}

#Save model simulations
save(outcomes, file = "IncrSchoolVacc_ES.RData")