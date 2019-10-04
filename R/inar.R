# inar.R

#' Model training
#'
#' Computes the posterior distribution
#'
#' @return A model
#' 
#' @export
inar <- function(time_series,
                   p = 1,
                   prior = list(a_alpha = NULL,
                                a_lambda = NULL, 
								b_lambda = NULL),
                   burn_in = 10^3,
                   chain_length = 10^4,
                   random_seed = 1761,
                   verbose = TRUE)
{
    if (!length(time_series) > 0) stop("Time series must have positive length.")
    if (length(time_series) <= p) stop("Time series length must be bigger than p.")
    if (!burn_in >= 0) stop("Burn-in must be positive.")
    if (!chain_length >= 0) stop("Chain length must be positive.")
    if (!random_seed >= 0) stop("Random seed must be positive.")
    if (!is.logical(verbose)) stop("Verbose parameter must be TRUE or FALSE.")
    
    set.seed(random_seed)
    
    begin <- proc.time()
        
    if (is.null(prior[["a_lambda"]]) || is.null(prior[["b_lambda"]])) {
        prior[["a_lambda"]] <- 1
        prior[["b_lambda"]] <- 1
    }
    
    if (is.null(prior[["a_alpha"]])) prior[["a_alpha"]] <- rep(1, p)
    
    post <- .posterior_inar(time_series,
                       p,
                       prior[c("a_alpha", "a_lambda", "b_lambda")],
                       burn_in,
                       chain_length,
                       random_seed,
                       verbose)
    
    model <- list(time_series = time_series,
                  p = p,
                  prior = prior[c("a_alpha", "a_lambda", "b_lambda")],
                  burn_in = list(length = burn_in,
                                 alpha = post$burn_in$alpha,
                                 lambda = post$burn_in$lambda),
                  chain = list(length = chain_length,
                               alpha = post$chain$alpha,
                               lambda = post$chain$lambda),
                  est = list(alpha = apply(post$chain$alpha, 2, mean),
                             lambda = mean(post$chain$lambda)),
                  elapsed = proc.time() - begin)
    
    class(model) <- "inar"
    
    invisible(model)
}

predict.inar <- function(model, h = 1, replications = 10^4) {
    pred <- .predictive_distribution_inar(model, h, replications)
    list(est = .generalized_median(pred), distr = pred)
}

cross_validate_inar <- function(time_series,
                           p = 1,
                           h = 1,
                           training_epoch = NULL,
                           prior = list(a_alpha = NULL,
                                        a_lambda = NULL, b_lambda = NULL),
                           burn_in = 10^3,
                           chain_length = 10^4,
                           random_seed = 1761,
                           verbose = TRUE)
{
    if(is.null(training_epoch)) training_epoch <- round(0.7 * length(time_series))
    
    # y_hat <- numeric(length(time_series) - (training_epoch + h) + 1)
    y_hat <- matrix(NA, nrow = length(h), ncol = length(time_series) - (training_epoch + min(h)) + 1)
    diff_y_hat_obs <- matrix(NA, nrow = length(h), ncol = length(time_series) - (training_epoch + min(h)) + 1)
    
    for (j in 1:ncol(y_hat)) {
        if (verbose) cat(sprintf("Training up to epoch %d ...\n", training_epoch + j - 1))
        model <- inar(time_series[1:(training_epoch + j - 1)],
                        p = p,
                        prior = prior,
                        burn_in = burn_in,
                        chain_length = chain_length,
                        random_seed = random_seed,
                        verbose = verbose)
        
        for (i in 1:nrow(y_hat)) {
            if (training_epoch + j + i - 1  <= length(time_series)) {
                if (verbose) cat(sprintf("Predicting epoch %d (h = %d)...\n", training_epoch + j + i - 1, h[i]))
                y_hat[i, j + i - 1] <- predict(model, h = h[i])$est
                diff_y_hat_obs[i, j + i - 1] = abs(time_series[training_epoch + j + i - 1] - y_hat[i, j + i - 1])
            }
        }
    }
    
    mae <- numeric(length(h))
    mae <- rowMeans(diff_y_hat_obs, na.rm = TRUE)
    names(mae) <- paste0("MAE (h = ", h, ")")
    
    list(est = y_hat, mae = mae)
}
#' Model summaries INAR(p)
#'
#' Summarises the model
#'
#' @return A summary
#' 
#' @export
summary.inar <- function(model) {
    printf <- function(...) cat(sprintf(...))
    printf("\n========================\n")
    printf("INAR(%d) Model Summary\n", model$p)
    printf("========================\n")
    printf("Time series length: %d\n", length(model$time_series))
    printf("Prior parameters:\n")
    printf("  a_lambda = %.2f, b_lambda = %.2f\n", model$prior[["a_lambda"]], model$prior[["b_lambda"]])
    cat("  ")
    for (i in 1:model$p) {
        printf("a_%d = %.2f", i, model$prior[["a_alpha"]][i])
        if (i < model$p) cat(", ")
        else cat("\n")
    }
    printf("Effective Markov chains length: %d (%d burn-in)\n", model$chain$length, model$burn_in$length)
    printf("Some posterior means with 95%% credible intervals:\n")
    for (i in 1:model$p) {
        post_qs <- unname(quantile(model$chain$alpha[, i], probs = c(0.025, 0.975)))
        printf("  alpha_%d = %.2f  [ %.2f ; %.2f ]\n", i, mean(model$chain$alpha[, i]), post_qs[1], post_qs[2])
    }
    post_qs <- unname(quantile(model$chain$lambda, probs = c(0.025, 0.975)))
    printf("  lambda_%d = %.2f  [ %.2f ; %.2f ]\n", i, mean(model$chain$lambda), post_qs[1], post_qs[2])
    
    printf("Total simulation time: %.2f seconds\n\n", round(model$elapsed[3]))
}