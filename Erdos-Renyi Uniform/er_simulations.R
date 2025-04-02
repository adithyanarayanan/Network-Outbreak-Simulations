## Date Last Modified: 4/2/25
## Author: Adithya Narayanan


## Task: Simulate outbreaks, and augmented novel tracing strategies, on Erdos Renyi Networks with 500 nodes.


# Note: Retaining seed values as is will replicate the exact simulations. Changing them might alter final numbers slightly, but will still produce similar results.
# While we alter the outbreak step function from the covidhm library in this script, we do still need the library to be installed and locally available to access 
# other functionality, and data within the library.

rm(list = ls())

### Import Libraries

library(covidhm) 
library(igraph)


### Set Values for the following simulation parameters in this portion of the script
### To modify the script to alter values that are not in this portion of the script - please reach out to author.


## Pick a strategy (options: p- primary tracing, s - secondary tracing, sp - selective primary tracing, 
## ss - selective secondary tracing, ps - pseudo-secondary tracing)
strategy = "p"

## Number of days until which to simulate an outbreak
out_days = 251

## Number of Networks, and Number of Outbreaks simulated per network
# numnetworks = 100 
# There is a nuance to changing the number of networks parameter which requires some trial and error. 
# If you would like to discuss, please reach out to the author. For the ER Simulations, this is hard-coded to be 100 networks. 
numoutbreaks = 100

## Each outbreak in each network produces a .csv file which counts the number of isolated, infectious, and quarantined individuals at a day level granularity. 
## For instance, network 1, outbreak 3, would produce - counts_1_3.csv in a folder. Aggregation across all such files produces cumulative data. 
## Set the name of the directory in which these files will be stored. Ensure that the directory exists before beginning simulation.
countsdirectory = "countsdirectory"

## Set up tracing strategies based on choice of strategy above
if(strategy == "p"){
quarantine = TRUE 
isolation = TRUE
tracing = TRUE
secondary = FALSE
}

if(strategy == "s"){
quarantine = TRUE
isolation = TRUE
tracing = TRUE
secondary = TRUE
}

if(strategy == "sp"){
quarantine = TRUE
isolation = TRUE
tracing = TRUE
secondary = FALSE
}

if(strategy == "ss"){
quarantine = TRUE
isolation = TRUE
tracing = TRUE
secondary = TRUE
}

if(strategy == "ps"){
quarantine = TRUE
isolation = TRUE
tracing = TRUE
secondary = TRUE
}





### Outbreak Step Function (source: covidhm library Firth et al.)
## One key modification made to this function is to account for intersections of the highly influential nodes (top_priorities) 
# with the primary or secondary strategies to simulate the novel augmented strategies

outbreak_step_new <- function(day, case_data, net = haslemere,
                              prop.asym, incfn, delayfn,
                              prop.ascertain, presymrate, R, quarantine,
                              isolation, tracing,
                              secondary, outside,
                              testing = FALSE,
                              cap_max_tests = NULL,
                              test_neg) {
  
  
  # Rename some variables ---------------------------------------------------
  
  newnet <- net
  colnames(newnet) <- c("caseid","contact","rate")
  
  
  
  # Update isolation, quarantine and infection status -----------------------------------
  if(isolation) {
    
    new_isolations <- which(day > case_data$isolated_time)
    if(quarantine){
      new_quarantines <- which(day > case_data$quarantine_time)
    }
    new_releases <- which(day > case_data$release_time)
    
    
    case_data$isolated[new_isolations] <- TRUE
    if(quarantine){case_data$quarantined[new_quarantines] <- TRUE}
    case_data$isolated[new_releases] <- FALSE
    case_data$quarantined[new_releases] <- FALSE
    
    
    
    #Reset isolation, release and test time for individuals who didn't undergo full isolation or quarantine
    early_releases <- c(new_releases[case_data$release_time[new_releases] <
                                       (case_data$isolated_time[new_releases] + 14)],
                        new_releases[case_data$release_time[new_releases] <
                                       (case_data$quarantine_time[new_releases] + 14)])
    case_data$isolated_time[early_releases] <- Inf
    case_data$test_time[early_releases] <- Inf
    case_data$release_time[early_releases] <- Inf
    case_data$quarantine_time[early_releases] <- Inf
  }
  
  #Assign recovered status to recovered individuals
  recovered <- which(day > case_data$recovery_time)
  case_data$status[recovered] <- "R"
  
  
  
  
  
  
  # Add infections from outside ---------------------------------------------
  
  if(outside > 0) {
    potential_new_inf<- which(case_data$status == "S" & !case_data$isolated & !case_data$quarantined)
    new_infections <- potential_new_inf[rbernoulli(length(potential_new_inf),
                                                   p = outside)]
    
    case_data$exposure[new_infections] <- day
    case_data$onset[new_infections] <- day + incfn(length(new_infections))
    case_data$recovery_time[new_infections] <- case_data$onset[new_infections] + 7
    case_data$status[new_infections] <- "I"
    
    #Isolation times for symtpomatic new infections
    if(isolation)
    {
      sym_cases <- new_infections[!case_data$asym[new_infections]]
      case_data$isolated_time[sym_cases] <- case_data$onset[sym_cases] +
        delayfn(length(sym_cases))
      case_data$release_time[sym_cases] <- case_data$isolated_time[sym_cases] + 14
    }
  }
  
  
  # Pull out infectious inds who are not isolated
  infectors <- dplyr::filter(case_data, status == "I", !isolated)
  
  
  
  
  
  
  # New cases  ----------------------------
  
  #Get contacts of infectious inds from network
  new_cases <- dplyr::filter(newnet,
                             caseid %in% infectors$caseid,
                             rate > 0)
  
  new_inf_rows <- match(new_cases$contact,
                        case_data$caseid)
  
  #Only keep susceptible contacts who are not isolated
  new_cases <- new_cases[case_data$status[new_inf_rows] == "S" &
                           !case_data$isolated[new_inf_rows] &
                           !case_data$quarantined[new_inf_rows],]
  
  
  
  
  
  
  # Generate new infections -------------------------------------------------
  
  if(nrow(new_cases) > 0) {
    
    #Filter based on probability that each contact is infected
    infector_rows <- match(new_cases$caseid,
                           case_data$caseid)
    
    asymrate <- ifelse(case_data$asym[infector_rows],0.5,1)
    
    infected <- rbernoulli(nrow(new_cases),
                           p = inf_prob(day = rep(day,length(infector_rows)),
                                        inc_samp = case_data$onset[infector_rows],
                                        contactrate = new_cases$rate,
                                        theta = presymrate,
                                        infasym = asymrate,
                                        R = R))
    
    #Each contact can only be infected once
    new_cases <- new_cases[infected,] %>%
      group_by(contact) %>%
      slice_sample(n = 1)
    
  }
  
  
  
  
  
  
  # Compile data for all new infections -------------------------------------
  
  
  # Compile a data frame for all new cases, new_cases is the amount of people that each infector has infected
  
  if(nrow(new_cases) > 0){
    
    prob_samples <- match(new_cases$contact,case_data$caseid)
    case_data$infector[prob_samples] <- case_data$caseid[match(new_cases$caseid,
                                                               case_data$caseid)]
    case_data$exposure[prob_samples] <- day
    case_data$status[prob_samples] <- "I"
    case_data$onset[prob_samples] <- day + incfn(length(prob_samples))
    case_data$recovery_time[prob_samples] <- case_data$onset[prob_samples] + 7
    
    #Isolation times for symtpomatic cases
    if(isolation) {
      sym_cases <- prob_samples[!case_data$asym[prob_samples]]
      case_data$isolated_time[sym_cases] <- case_data$onset[sym_cases] +
        delayfn(length(sym_cases))
      case_data$release_time[sym_cases] <- case_data$isolated_time[sym_cases] + 14
      
    }
    
  }
  
  
  
  
  
  
  # Contact tracing ---------------------------------------------------------
  
  #Empty vector of contacts
  traced_contacts <- c()
  
  
  if(tracing) {
    #get contacts of symptomatic infectors who onset yesterday
    new_ill <- filter(case_data,
                      status == "I",
                      !asym,
                      onset < day,
                      onset > (day - 1))
    
    if(nrow(new_ill) > 0){
      
      #Get contacts of infectious inds from network
      case_contacts <- filter(newnet,
                              caseid %in% new_ill$caseid,
                              rate > 0)
      
      new_contact_rows <- match(case_contacts$contact,
                                case_data$caseid)
      
      #Only keep contacts who are not isolated and not recovered
      traced_contacts <- case_contacts[case_data$status[new_contact_rows] != "R" &
                                         !case_data$isolated[new_contact_rows] &
                                         !case_data$quarantined[new_contact_rows],]
      traced_contacts <- traced_contacts$contact[rbernoulli(length(traced_contacts$contact),
                                                            p = prop.ascertain)]
    }
    
  }
  
  
  ## Secondary contact tracing -------------
  
  
  
  # Executed only when secondary contacts have to be traced as per choices made above
  if(secondary) {
    if(nrow(new_ill) > 0) {
      if(nrow(case_contacts) > 0) {
        
        sec_contacts <- filter(newnet,
                               caseid %in% case_contacts$contact,
                               rate > 0)
        
        new_contact_rows <- match(sec_contacts$contact,
                                  case_data$caseid)
        
        
        #Only keep secondary contacts who are not isolated and not recovered
        traced_sec_contacts <- sec_contacts[case_data$status[new_contact_rows] != "R" &
                                              !case_data$isolated[new_contact_rows] &
                                              !case_data$quarantined[new_contact_rows],]
        
        #Filter based on ascertainment rate
        traced_sec_contacts <- traced_sec_contacts$contact[rbernoulli(length(traced_sec_contacts$contact),
                                                                      p = prop.ascertain)]
        

        # If Pseudo-secondary tracing, trace only most influential from secondary contacts
        if(strategy == "ps"){
            traced_sec_contacts = intersect(traced_sec_contacts, top_priorities)
        }
        
        
        #Add to list of traced contacts
        traced_contacts <- unique(c(traced_contacts,
                                    traced_sec_contacts))
        
      }
    }
  }
  
  # If selective primary tracing, trace only most influential from primary contacts
  if(strategy == "sp"){
    traced_contacts = intersect(traced_contacts, top_priorities)
  }

  # If selective secondary tracing, trace only most influential from primary and secondary contacts
  if(strategy == "ss"){
    traced_contacts = intersect(traced_contacts, top_priorities)
  }


  # Set quarantine times based on tracing -----------------------------------
  
  if(isolation & tracing & length(traced_contacts) > 0) {
    
    if(quarantine) {
      #If you are recovered and asymptomatic and traced, you isolate
      recovered_traced <- which(case_data$caseid %in% traced_contacts &
                                  case_data$status == "R" &
                                  case_data$asym)
      case_data$quarantine_time[recovered_traced] <- day+delayfn(length(recovered_traced))
      
      #If you are susceptible you quarantine on being traced
      susceptible_traced <- which(case_data$caseid %in% traced_contacts &
                                    case_data$status == "S")
      
      case_data$quarantine_time[susceptible_traced] <- day+delayfn(length(susceptible_traced))
      
      #If you are infectious you isolate if the delay is shorter than your current iso time
      infectious_traced <- which(case_data$caseid %in% traced_contacts &
                                   case_data$status == "I")
      new_iso_time <- day + delayfn(length(infectious_traced))
      case_data$quarantine_time[infectious_traced] <-  day + delayfn(length(infectious_traced))
    } else {
      
      #If you're susceptible and there's no quarantine and you're traced
      #You isolate immediately at symptom onset
      susceptible_traced <- which(case_data$caseid %in% traced_contacts &
                                    case_data$status == "S" &
                                    !case_data$asym)
      
      case_data$isolated_time[susceptible_traced] <- case_data$onset[susceptible_traced]
      
      #If you're infectious and traced and no quarantine, you also isolate at onset
      infectious_traced <- which(case_data$caseid %in% traced_contacts &
                                   case_data$status == "I" &
                                   !case_data$asym)
      
      case_data$isolated_time[infectious_traced] <- case_data$onset[infectious_traced]
      
    }
    
    case_data$release_time <- ifelse(case_data$isolated_time < case_data$quarantine_time,
                                     case_data$isolated_time + 14,
                                     case_data$quarantine_time + 14)
    
  }
  
  
  
  # Test and release quarantined and isolated cases ------------------------------------
  
  if(testing)
  {
    #Identify untested inds and assign a test time based on their isolation and quarantine time
    untested <- which(case_data$test_time == Inf)
    case_data$test_time[untested] <- ifelse(case_data$isolated_time[untested] < case_data$quarantine_time[untested],
                                            case_data$isolated_time[untested],
                                            case_data$quarantine_time[untested]) +
      delayfn(length(untested))
    
    #Today's tests
    new_tests <- which(day < case_data$test_time &
                         (day+1) > case_data$test_time)
    
    #If there's a test cap randomly sample tests to be excluded and reset their test time
    if(length(new_tests) > cap_max_tests) {
      not_tested <- new_tests[sample(1:length(new_tests),length(new_tests)-cap_max_tests)]
      case_data$test_time[not_tested] <- Inf
      new_tests <- new_tests[!(new_tests %in% not_tested)]
    }
    
    
    #Run tests
    #90% chance of testing positive
    #2% false positive rate
    
    test_results <- (case_data$status[new_tests] == "I" &
                       rbernoulli(length(new_tests),
                                  1-test_neg)) |
      rbernoulli(length(new_tests), 0.02)
    
    
    negative_tests <- new_tests[!test_results]
    
    case_data$release_time[negative_tests] <- case_data$test_time[negative_tests]
    
  }
  
  
  
  
  return(case_data)
}





### Generate Networks, track outbreak metrics

isolations = c()
quarantines = c()
infections = c()
network = c()
simulation = c()


# Identify 100 usable ER Networks for our simulations
# The igraph function erdos.renyi.game generates networks, sometimes with isolated nodes.
# From trial and error, it was identified that when seeded between 1:117, we have 100 networks which don't have isolated nodes. 
# In this portion, we identify the 100 networks which are usable, and discard the rest.
whicher = c() 

n = seq(1:117)
for(i in n){
  set.seed(i)
  karate = erdos.renyi.game(500, 1999, type = "gnm")
  if(length(which(degree(karate) == 0)) == 0){
    whicher = c(whicher, i)
  } 
}

print("Identified 100 Erdos-Renyi Networks")

ms = seq(1:117) # Tracker for current network
is = seq(1:numoutbreaks) # Tracker for current outbreak (i'th outbreak on m'th network)


# Generate specified number of networks - Barabasi Albert Preferential Attachment Networks
# All networks are generated with 500 nodes by default

for(m in ms){

  if(!m %in% whicher){next} # Skip networks which have isolated nodes
  
  set.seed(m)
  kgraph = erdos.renyi.game(500, 1999, type = "gnm")
  nodes = 500
  
  
  print(paste("Network", as.character(m), sep = " "))
  
  
  V(kgraph)$color = rep(0, nodes)
  V(kgraph)$size = rep(2, nodes)
  V(kgraph)$label = NA
  E(kgraph)$size = 0.1
  
  
  kgraph_mat = as_adjacency_matrix(
    kgraph,
    type = c("both"),
    attr = NULL,
    edges = FALSE,
    names = TRUE,
    sparse = igraph_opt("sparsematrices")
  )
  kgraph_mat = matrix(kgraph_mat, nrow=nodes, ncol=nodes)
  
  
  ###### Getting PRINCE priorities
  
  g = kgraph
  W <- kgraph_mat
  deg.vec <- colSums(W)
  set.seed(29)
  Y <- runif(nrow(kgraph_mat))
  n <- nrow(W)
  
  D <- matrix(0, nrow = n, ncol = n)
  D_sqr <- D
  diag(D_sqr) <- 1/sqrt(deg.vec)
  W_prime <- D_sqr %*% W %*% D_sqr
  
  F_old = Y
  F_store <- c(F_old)
  alph = .95
  for (i in 1:100){
    F_new = alph*W_prime %*% F_old + (1-alph)*Y
    F_store  = cbind(F_store, F_new)
    F_old <- F_new
  }
  
  hub_weights <- F_store[, i]
  num = 125
  top_priorities <- match(head(sort(hub_weights,decreasing = T),num),hub_weights) # stores the 125 most influential nodes in our generated network
  


  ## Set Up Outbreak Parameters  
  kgraph_formatted = format_network(kgraph_mat)
  
  incfn <- dist_setup(dist_shape = 2.322737,dist_scale = 6.492272)
  delayfn <- dist_setup(1,1.4)
  
  
  isolated = c()
  quarantined = c()
  infected = c()

  ## Begin simulating i outbreaks on the current network out of m networks  
  for(i in is){
    print(paste("Simulation", as.character(i), sep = " "))
    
    
    set.seed(i)
    case_data <- outbreak_setup(net = kgraph_formatted, num.initial.cases = 1,incfn,delayfn, prop.asym= 0.2, isolation = TRUE)
    case_data = arrange(case_data, caseid)
    
    #### Simulating COVID
    
    day = c()
    exposed = c()
    
    days = out_days
    cday = 0
    
    
    set.seed(29)
    
    ise_set = c()
    quar_set = c()
    inf_set = c()
    
    
    iso_today = c()
    quar_today = c()
    inf_today = c()
    
    while(cday <= days){
      
      case_data = outbreak_step_new(day = cday, case_data = case_data, 
      net = kgraph_formatted, prop.asym = 0.2, incfn = incfn, delayfn = delayfn, 
      prop.ascertain = 1, presymrate = 0.2, R = 0.8, 
      quarantine = TRUE, isolation = TRUE, tracing = TRUE, secondary = secondary, 
      outside = 0, testing = FALSE,test_neg = 0.1)
      case_data = arrange(case_data, caseid)
      
      
      ##############################################
      ## Notify progress of simulation every 10 days. Can comment out if verbose output it not required.
      
      if(cday%%10 == 0){
        print(paste("Day", as.character(cday), sep = " "))}
      
      
      ###########################################
      
      
      
      color = ifelse(case_data$isolated == TRUE, "red",
                     ifelse(case_data$quarantined == TRUE, "yellow", "white"))
      
      
      k = data.frame(case_data$caseid, case_data$isolated, case_data$quarantined, case_data$status, color)
      # print(which(k$case_data.quarantined == TRUE))
      # print(which(k$case_data.isolated == TRUE))
      
      ise_set = c(ise_set, which(k$case_data.isolated))
      quar_set = c(quar_set, which(k$case_data.quarantined))
      inf_set = c(inf_set, which(k$case_data.status == "I"))
      
      
      iso_current = which(k$case_data.isolated)
      quar_current = which(k$case_data.quarantined)
      quar_current = setdiff(quar_current, iso_current)
      
      iso_today = c(iso_today, length(iso_current))
      quar_today = c(quar_today, length(quar_current))
      inf_today = c(inf_today, length((which(k$case_data.status == "I"))))
      
      V(kgraph)$color = k$color
      
      cday = cday + 1
    }
    
    days_list = seq(0,days)
    simulations_id = rep(i, days+1)
    counts_of_stats = data.frame(days_list, simulations_id, iso_today, quar_today, inf_today)
    
    
    filename = paste0(countsdirectory, "/counts_", m, "_", i, ".csv")
    write.csv(counts_of_stats, filename)
    
    
    isolations = c(isolations, length(unique(ise_set)))
    quarantines = c(quarantines, length(unique(quar_set)))
    infections = c(infections, length(unique(inf_set)))
    network = c(network, m)
    simulation = c(simulation, i)
    
  }
  
  
  
  
}

## Print out summary of outbreak metrics from all simulations
all_data = data.frame(isolations, quarantines, infections, network, simulation)
summary(all_data)

sd(isolations)
sd(quarantines)
sd(infections)

## Optionally write all_data dataframe to a csv

# filename2 = "ps_primary_250_extra.csv"
# write.csv(all_data, filename2)

