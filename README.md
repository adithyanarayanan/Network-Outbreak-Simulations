# Network-Outbreak-Simulations
Code to recreate the simulations from 

_Narayanan A, Muldoon SF, Jehrio M, Hageman Blair R. Informing pandemic intervention strategies
through coupled contact tracing and network node prioritization. [In Revisions: March 2025]. PLOS
Complex Systems.
_

This repository contains code that produces the simulations described in the study cited above. Each directory contains the scripts to generate (where necessary) a unique type of network(s) and simulates outbreaks on them, which are then curtailed with the use of several existing and proposed novel augmented contact tracing strategies. 

## Abstract

SARS-CoV-2 has highlighted the challenges of social intervention measures for disease control, which are difficult to implement and highly disruptive to modern society. Simulation models have demonstrated the efficacy of primary and secondary tracing at SARS-CoV-2 disease control but at the cost of quarantining large proportions of the population. This paper develops novel tracing strategies that harness node (individual) influence in a social association network for contact-tracing approaches to disease
control. The overarching assumption is that an individual’s potential to spread diseasecan be modeled by their ability to propagate influence through a network. Models ofidea and influence propagation have been widely studied in the context of social networks but have limited application to disease models. The PRIoritization and Complex Elucidation (PRINCE) algorithm is leveraged to estimate an individual node’s influence score that reflects their ability to propagate disease based on network
connectivity. In this study, we propose novel augmented tracing strategies that leverage a node’s influence to assist with targeted tracing in its 1-hop and 2-hop neighborhoods: i) pseudo-secondary tracing (tracing and quarantining the immediate contacts and the influential contacts of contacts of an infectious symptomatic individual) and ii) selective secondary tracing (tracing and quarantining the influential immediate contacts, and influential contacts of contacts of an infectious symptomatic individual). Contagion dynamics on simulated and real-world networks, benchmarked with existing strategies, demonstrate that our novel strategies mitigate societal disruption by lowering the maximum number of people quarantined concurrently while also assisting the ease of on-ground deployment by reducing the number of individuals to be traced for every infectious individual detected when compared with the most effective existing tracing strategy. Novel approaches of this type that embed network influence into pandemic control provide an opportunity for disease control that ultimately lessens the disruption to society.

## Usage

### Installation
This package leverages the work of 
_Firth, Josh A., et al. "Using a real-world network to model localized COVID-19 control strategies." Nature medicine 26.10 (2020): 1616-1622.
_
and the associated code for package/library covidhm published at 
https://github.com/biouea/covidhm

Before commencing any work, please download the code from this repository locally to be able to install it. Please do not follow the instructions from the _covidhm_ repo as the dependencies have changed some of their functionality and may not allow covidhm to function accurately. 

Instead, after downloading the covidhm files, and in a clean R environment with no packages/libraries except the base installed, please run the script install_packages.R. If at any step R prompts you to update any of the packages you try to install, please skip updates. 


### Simulations
Each directory houses a simulation script which runs a unique kind of network simulations. Within each such script, please find in the first section, parameters that you can tune to simulate outbreaks as per your specifications. To simulate the results from the study, please use the study parameters (specified within scripts as the default values). 

Each script contains specific instructions which are unique to the network and simulation type and the comments can guide their execution. Each execution of the scripts takes hours and produces several csv files to store recorded data which is stored in a subdirectory by default, but can be modified to be stored elsewhere. 

Each directory also contains a second script, which helps generate a unique visualization for each network-outbreak configuration. This is generated from the data stored from the simulation scripts.

If any of the scripts do not work as expected, or raise any issues, please reach out to me and I will help with setting the simulations up. 

