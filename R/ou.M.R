# Function from paleoTS (License:	GPL-3)

ou.M <-
function(anc, theta, aa, tt) theta*(1 - exp(-aa*tt)) + anc*exp(-aa*tt)
