library("markovchain")

nStates <- 3
P <- matrix(rep(1/nStates, nStates * nStates),
             nrow=nStates, byrow=TRUE)
stateNames <- c("Rain", "Nice", "Snow")
rownames(P) <- stateNames
colnames(P) <- stateNames

pi <-  c(2/5, 1/5, 2/5)
s <- matrix(rep(1, nStates * nStates),
            nrow=nStates, byrow=TRUE)

alpha <- matrix(0, nStates,nStates)
for (i in 1:nStates) {
    for (j in 1:nStates) {
        alpha[i, j] <- s[i, j] * (1 /( 1 + (pi[i]/pi[j])*(P[i,j]/P[j,i])))
    }
}

Q <- matrix(0, nStates,nStates)
for (i in 1:nStates) {
    for (j in 1:nStates) {
        if (i != j) {
            Q[i, j] <- P[i, j] * alpha[i, j]
        }
    }
    Q[i, i] <- 1 - sum( Q[i, ] ) 
}

hastMetWeather <- new("markovchain", states = stateNames,
                   transitionMatrix = Q,
                 name = "hastMetWeather")
cat( "Steady State Distribution: ", steadyStates(hastMetWeather), "\n" )

start <- 900
end <- 1000
hastMetWeatherSim <- rmarkovchain(n=end, object=hastMetWeather,
                                  t0="Nice")
cat("Simulation Results\n")
cat( "Rain:", sum( hastMetWeatherSim[(start+1):end] == "Rain")/(end-start), "\n")
cat( "Nice:", sum( hastMetWeatherSim[(start+1):end] == "Nice")/(end-start), "\n")
cat( "Snow:", sum( hastMetWeatherSim[(start+1):end] == "Snow")/(end-start), "\n")

## NAME: hastMet.R
##
## USAGE: within R, at interactive prompt: source("hastMet.R")
##        or at command line: Rscript hastMet.R
## REQUIRED ARGUMENTS: none
## OPTIONS: none
## DESCRIPTION: Illustration of the Hastings-Metropolis algorithm on
##         the simple three-state weather Markov chain.  The
##         illustration is to closely mimic the mathematical
##         statement of the algorithm without regard to efficiency.
## DIAGNOSTICS: none
## CONFIGURATION AND ENVIRONMENT: requires R markovchain library
## DEPENDENCIES:  requires R markovchain library
## INCOMPATIBILITIES: none known
## PROVENANCE: origianl from Steven Dunbar
## BUGS AND LIMITATIONS: none kown
## FEATURES AND POTENTIAL IMPROVEMENTS: make more efficient

## AUTHOR:  Steve Dunbar
## VERSION: Version 1.0 as of March 5, 2020
## KEYWORDS: Hastings-Metropolis algorithm, Markov chain




