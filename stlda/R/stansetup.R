stansetup <-
function () 
{
    rstan_options(auto_write = TRUE)
    options(mc.cores = parallel::detectCores())
}
