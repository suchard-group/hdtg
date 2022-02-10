library(devtools)

d = 20
n_iter = 1000
burnin = 0.2
mu = rep(0, d)
set.seed(666)
cov_mat = rWishart(1, 2 * d, diag(d))[ , , 1]
cor_mat = cov2cor(cov_mat)

cor_mat = matrix(data = 0.9, nrow = d, ncol = d)
diag(cor_mat) = rep(1, d)
prec = solve(cor_mat)
p0 <- runif(d, 0, 0.2)

load_all(export_all=FALSE)
res = rcmg(n = n_iter, mean = mu, prec = prec, constraits = rep(1, d), t = 50, burnin = 0, p0 = p0, cpp_flg = T, nuts_flg = T, debug_flg = F)

plot(ts(res[3,(burnin*n_iter):n_iter]))
mean(res[1,(burnin*n_iter):n_iter])
