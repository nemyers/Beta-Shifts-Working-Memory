function ll = mle_mixpdf(x, mu, param)

K = param(1);
Pm = param(2);
ll = -sum(log(mixpdf(x, mu, K, Pm)));