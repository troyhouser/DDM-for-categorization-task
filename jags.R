df = read.csv("~/Dropbox (University of Oregon)/TPC_DDM/TPC_ddm.csv")
df2 = cbind(df$subj_idx,df$response,df$rt,df$dist)
colnames(df2) = c("id","response","rt","condition")
df2 = as.data.frame(df2)
df2$condition = as.factor(df2$condition)
df2 = na.omit(df2)
df2 = df2[!df2$rt<0.2,]
library(brms)
formula = bf(rt | dec(response) ~ 0 + condition + (1|p|id),
             bs ~ 0 + condition + (1|p|id),
             ndt ~ 0 + condition + (1|p|id),
             bias = 0.5)

prior = get_prior(formula,
          data = df2, 
          family = wiener(link_bs = "identity", 
                          link_ndt = "identity", 
                          link_bias = "identity"))
tmp_dat <- make_standata(formula, 
                         family = wiener(link_bs = "identity", 
                                         link_ndt = "identity",
                                         link_bias = "identity"),
                         data = df2, prior = prior)
str(tmp_dat, 1, give.attr = FALSE)

initfun <- function() {
  list(
    b = rnorm(tmp_dat$K),
    b_bs = runif(tmp_dat$K_bs, 1, 2),
    b_ndt = runif(tmp_dat$K_ndt, 0.1, 0.15),
    #b_bias = rnorm(tmp_dat$bias, 0.5, 0.1),
    sd_1 = runif(tmp_dat$M_1, 0.5, 1),
    z_1 = matrix(rnorm(tmp_dat$M_1*tmp_dat$N_1, 0, 0.01),
                 tmp_dat$M_1, tmp_dat$N_1),
    L_1 = diag(tmp_dat$M_1)
  )
}
fit_wiener <- brm(formula, 
                  data = df2,
                  family = wiener(link_bs = "identity", 
                                  link_ndt = "identity",
                                  link_bias = "identity"),
                  prior = prior, init = initfun,
                  iter = 10000, warmup = 4000, 
                  chains = 4, cores = 4,
                  control=list(adapt_delta=0.99))
fit_wiener
NPRED <- 500
pred_wiener <- predict(fit_wiener, 
                       summary = FALSE, 
                       negative_rt = TRUE, 
                       nsamples = NPRED)