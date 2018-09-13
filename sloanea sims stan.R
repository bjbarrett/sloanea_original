library(rethinking)
library(rstan)
library(truncnorm)

setwd("~/Dropbox/Lomas Barbudal/Sloanea Data")
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

Softmax <- function(x){
exp(x)/sum(exp(x))
} #softmax function to simplify code

#initial values
n <- 20 # number of individuals
nbouts <- 100 #rounds foraged
##parameters
k.lambda <- 5
k.gam <- -1.3	 # logit weight of social info for group-- gamma in model
k.phi <- -1        # logit stickiness parameter (phi)--alpha in model
k.gam_i <- rnorm( n , mean=0 , sd=1 ) #varying effect weight of social information for individuals--g_indiv in model
k.phi_i <- rnorm( n , mean=0 , sd=1.5 ) #varying effects component of stickiness of attraction scores (phi)
k.fconf <- 1 #conformity inflection point 0 is 1; 1 is 2.7 , negative is anticonformity
k.fconf_i <- rnorm( n , mean=0 , sd=.5 ) #varying effects component of stickiness of attraction scores ()


group <- c(1:n)
numforg <- n
numobs <- n


techpr_i <- c(.25,.25,.25,.25)



#####SIMULATE  INDIVIDUAL(change k.gam to ~ -15 to make individual learning)
dsim_s <- data.frame( i=0 , bout=0 , tech=0 , y1=0 , y2=0, y3=0 , y4=0 , s1=0 , s2=0 , s3=0 , s4=0, A1=0 , A2=0 , A3=0 , A4=0)
therow <- 1
A4 <- A3 <- A2 <- A1 <- rep(.25,n) #attraction scores for each tech
#A4 <- rep(10,n) #attraction scores for each tech

S1 <- S2 <- S3 <- S4 <- rep(0,n+1) # num of individuals choosing each tech in previous bout
for ( r in 1:nbouts ) {
    for ( i in 1:n ) {        
        prtech_i <-  Softmax(k.lambda*as.vector(c(A1[i],A2[i],A3[i],A4[i]) ) )
        prtech <- prtech_i 
# choose tech
        tech <- sample( 1:4 , size=1 , prob=prtech)
# update attractions
        yields <- rep(0,4)
        yields[tech] <- 1
        effg <- logistic(k.phi + k.phi_i[i] )
        A1[i] <- (1-effg)*A1[i] + effg*yields[1]
        A2[i] <- (1-effg)*A2[i] + effg*yields[2]
        A3[i] <- (1-effg)*A3[i] + effg*yields[3]
        A4[i] <- (1-effg)*A4[i] + effg*yields[4]
       
# record data
        dsim_s[therow,] <- c( i , r , tech , yields[1] , yields[2] , yields[3] , yields[4] , S1[r] , S2[r] , S3[r] , S4[r] , A1[i] , A2[i] , A3[i] , A4[i] )
        therow <- therow + 1
    } #i
    S1[r+1] <- length( dsim_s$tech[dsim_s$tech==1 & dsim_s$bout==r] )
    S2[r+1] <- length( dsim_s$tech[dsim_s$tech==2 & dsim_s$bout==r] )
    S3[r+1] <- length( dsim_s$tech[dsim_s$tech==3 & dsim_s$bout==r] )
    S4[r+1] <- length( dsim_s$tech[dsim_s$tech==4 & dsim_s$bout==r] )
} #r

# sort by individual then round (bout)
o <- order( dsim_s$i )
dsim2 <- dsim_s[o,]

plot(s1/n ~ bout, data=dsim2, col="red" , ylim=c(0,1))
points(s2/n ~ bout, data=dsim2 , col="gold")
points(s3/n ~ bout, data=dsim2 , col="green")
points(s4/n ~ bout, data=dsim2 , col="blue")


#####SIMULATE FREQUENCY DEPEDENT LEARNING AND INDIVIDUAL(change k.gam to ~ -15 to make individual learning)
dsim_s <- data.frame( i=0 , bout=0 , tech=0 , y1=0 , y2=0, y3=0 , y4=0 , s1=0 , s2=0 , s3=0 , s4=0, A1=0 , A2=0 , A3=0 , A4=0)
therow <- 1
A4 <- A3 <- A2 <- A1 <- rep(0,n) #attraction scores for each tech
S1 <- S2 <- S3 <- S4 <- rep(0,n+1) # num of individuals choosing each tech in previous bout
for ( r in 1:nbouts ) {
    for ( i in 1:n ) {        
        prtech_i <-  Softmax(k.lambda*as.vector(c(A1[i],A2[i],A3[i],A4[i]) ) )
        my.s <- logistic( k.gam + k.gam_i[i] ) #my.s is how social information is weighted
        prtech_su <- c(S1[r],S2[r],S3[r],S4[r])
        fconf <- exp( k.fconf + k.fconf_i[i]) 
        if (sum(prtech_su)>0){
            # social learning
            prtech_s <- prtech_su^fconf/sum(prtech_su^fconf)
            prtech <- (1-my.s) * prtech_i + my.s * prtech_s
        }else{
            prtech <- prtech_i 
        }
# choose tech
        tech <- sample( 1:4 , size=1 , prob=prtech)
# update attractions
        yields <- rep(0,4)
        yields[tech] <- 1
        effg <- logistic(k.phi + k.phi_i[i] )
        A1[i] <- (1-effg)*A1[i] + effg*yields[1]
        A2[i] <- (1-effg)*A2[i] + effg*yields[2]
        A3[i] <- (1-effg)*A3[i] + effg*yields[3]
        A4[i] <- (1-effg)*A4[i] + effg*yields[4]
       
# record data
        dsim_s[therow,] <- c( i , r , tech , yields[1] , yields[2] , yields[3] , yields[4] , S1[r] , S2[r] , S3[r] , S4[r] , A1[i] , A2[i] , A3[i] , A4[i] )
        therow <- therow + 1
    } #i
    S1[r+1] <- length( dsim_s$tech[dsim_s$tech==1 & dsim_s$bout==r] )
    S2[r+1] <- length( dsim_s$tech[dsim_s$tech==2 & dsim_s$bout==r] )
    S3[r+1] <- length( dsim_s$tech[dsim_s$tech==3 & dsim_s$bout==r] )
    S4[r+1] <- length( dsim_s$tech[dsim_s$tech==4 & dsim_s$bout==r] )
} #r
# sort by individual then round (bout)
o <- order( dsim_s$i )
dsim2 <- dsim_s[o,]

plot(s1/n ~ bout, data=dsim2, col="red" , ylim=c(0,1))
points(s2/n ~ bout, data=dsim2 , col="gold")
points(s3/n ~ bout, data=dsim2 , col="green")
points(s4/n ~ bout, data=dsim2 , col="blue")

datalist_s <- list(
                    N = nrow(dsim2),
                    J = length( unique(dsim2$i) ),
                    K = length(unique(dsim2$tech) ),
                    tech = dsim2$tech,
                    y = cbind( dsim2$y1 , dsim2$y2 , dsim2$y3 , dsim2$y4 ),
                    s = cbind( dsim2$s1 , dsim2$s2 ,dsim2$s3 , dsim2$s4 ),
                    bout = dsim2$bout,
                    id = dsim2$i,
                    N_effects=3
                    )
parlist=c("a_id" , "mu", "lambda","phi_i" , "gamma_i" , "fconf_i" , "dev" , "log_lik")
fit_s <- stan( file = 'PN_freq_social.stan', data = datalist_s , iter = 1400, warmup=700 , chains=2 , pars=parlist , control=list(adapt_delta=0.99) )


post <- extract(fit_s) #(0,1)
dens(logistic(post$mu[,1]) , main="phi", show.HPDI=0.90) #phi
abline(v=logistic(k.phi) , lw=3 )
dens((post$mu[,2]) , main="log(fconf)" , show.HPDI=0.90) #f-conf
abline(v=k.fconf , lw=3 )
dens(logistic(post$mu[,3]) , main="gamma" , show.HPDI=0.90) #gamma
abline(v=logistic(k.gam) , lw=3 )
dens(post$lambda) #gamma
abline(v=k.lambda , lw=3 )

par(mfrow = c(3, 1))

for (i in 1:n){
    dens(logistic(post$mu[,1] + post$a_id[,i,1]) , xlim=c(0,1) , main="phi", show.HPDI=0.90) #delta
    abline(v=logistic(k.phi + k.phi_i[i]) , lw=3 )
    dens(exp((post$mu[,2]+ post$a_id[,i,2])) , main="log(fconf)" , show.HPDI=0.90 , xlim=c(0,5)) #f-conf
    abline(v=exp(k.fconf + k.fconf_i[i]) , lw=3)
    dens(logistic(post$mu[,3] + post$a_id[,i,3]) , main="gamma" , show.HPDI=0.90, xlim=c(0,1)) #gamma
    abline(v=logistic(k.gam + k.gam_i[i]) , lw=3 ) 
}
  plot.new()



