# 
# Learning R
#

# Clear R's brain...
rm( list = ls() )

library(dplyr)
library(ggplot2)


x <- seq( from=-10, to=10, by=0.1 )

y <- x*x

z =  y * y

qplot( x, y, geom="line")



die = 1:6


generate = function()
{
   s = sample( die, size=2, replace=TRUE )
   
   return(sum(s))
}

Toss = replicate( 10000, generate() )


qplot( Toss, binwidth=1 )

