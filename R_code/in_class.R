# 
# BYS 602 Mouse weight study
#

# Clear R's brain...
rm( list = ls() )

library(dplyr)
library(ggplot2)


weights = read.csv( "femaleMiceWeights.csv" )

pop = read.csv( "ControlsPopulation.csv" )


col = 2
control = weights[1:12,2]

treatment = weights[13:24,2]
print(control)
control_avg = mean(control)
print( control_avg )

print(treatment)
treatment_avg = mean(treatment)
print( treatment_avg )

diff = treatment_avg - control_avg



#sample the population for a set of controls
nsamples = 10000
null = vector( "numeric", nsamples )
for( i in 1:nsamples)
{
   controlnull = sample( pop[,1], 12)
   treatmentnull = sample( pop[,1], 12)
   null[i] = mean(treatmentnull) - mean(controlnull)
   
}#end for


hist( null )
abline(v=diff)

# What percent bigger than the experiment diff?
m = mean( null >= diff )


print( m )

clsh = read.csv( "class_heights.csv")

data = clsh[,1]

u = mean( data )
ssqr = var( data )
sig = sqrt( ssqr )

print( u )
print( ssqr )
print( sig )




##aggregate(weights,)

# EOF
