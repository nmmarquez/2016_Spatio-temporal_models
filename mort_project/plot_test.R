inf_term <- function(x, N0=1, lambda=-1, c=0){
    N0 * exp(lambda * x) + c
}

logit_scaled <- function(x, eta, scale){
    1 / (1 + exp(scale * (eta - x)))
}

ya_term <- function(x, eta=mean(x), scale=1, ceiling=1){
    ceiling * logit_scaled(x, eta, scale)
}

sns_term <- function(x, rho, m, b){
    (m * x + b) * logit_scaled(x, eta=rho, scale=3)
}

plot(0:100, inf_term(0:100, lambda=-1, N0=10, c=1.2))
plot(0:100, ya_term(0:100, ceiling=3.8, scale=.2))
plot(0:100, sns_term(0:100, 30, .03, b=-1.2))

plot(0:100, inf_term(0:100, lambda=-1, N0=10, c=1.2) + 
           ya_term(0:100, ceiling=3.8, scale=.3, eta=17) + 
           sns_term(0:100, 30, .05, b=-1.65), type="l")
