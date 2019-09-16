##' treatment HR changes
##'
##' examine changes in main effect (HR, Cox regression) when covariates are
##'     added one at a time, sequentially increasing (starting from a small
##'     model) or sequentially decreasing (starting from large model)
##' @param data data frame
##' @param surv name of response 'Surv'-variable in data set
##' @param main name of main effect (binary)
##' @param terms vector of terms to adjust for (can be named)
##' @param decr.inc logical; decreasing order for sequential inclusion?
##' @param decr.exc logical; decreasing order for sequential exclusion?
##' @param rms use \code{rms::cph} instead of \code{survival::coxph}? (mainly if
##'     xtra.adj contains something specific to the rms package)
##' @export
coxreg_HRS <- function(data, surv, main, terms,
                       decr.inc = NULL, decr.exc = NULL,
                       rms = FALSE){
    ## set modeling function
    modelFNC <- if(rms){
                    ## rms::cph
                    warning("'rms' not implemented yet")
                    survival::coxph
                } else{
                    survival::coxph
                }
    ## check if terms have names
    if(is.null(names(terms))) names(terms) <- terms
    ## have data in data.frame format, so that it can keep a Surv object
    data <- as.data.frame(data, stringsAsFactors = FALSE)
    ## check surv argument
    if(length(surv) == 1){
        surv.name <- surv
        if(class(data[[surv.name]]) != 'Surv'){
            stop("'surv' is not a 'Surv' object in data")
        }
    } else if(length(surv) == 2){
        data$outcome <- survival::Surv(time = surv[1], event = surv[2])
        surv.name <- "outcome"
        bnry <- setdiff(bnry, "outcome")
        real <- setdiff(real, "outcome")
    } else stop("'surv' argument in strange form")
    ## model formula
    Sf <- paste0(surv, " ~ ", main)
    Lf <- paste0(surv, " ~ ", main, " + ",
                 paste0(terms, collapse = " + "))
    Smod <- modelFNC(formula(Sf), data = data)
    Lmod <- modelFNC(formula(Lf), data = data)
    ## update function
    upd <- function(f, s, op = " + ") paste0(f, op, s)
    ## get unadjusted effect of main
    est <- Smod$coefficients[1]
    se <-  sqrt(Smod$var[1,1])
    UNADJ <- data.frame(term = "(none)",
                        HR = exp(est),
                        HR.l = exp(est - 1.96*se),
                        HR.h = exp(est + 1.96*se))
    ## get adjusted effect of main
    est <- Lmod$coefficients[1]
    se <-  sqrt(Lmod$var[1,1])
    ADJ <- data.frame(term = "(all)",
                      HR = exp(est),
                      HR.l = exp(est - 1.96*se),
                      HR.h = exp(est + 1.96*se))
    ## get main + univariate adjustments HRs
    P <- UNADJ
    for(i in seq_along(terms)){ ## i = 1
        mtmp <- modelFNC(formula(upd(Sf, terms[i])), data = data)
        est <- mtmp$coefficients[1]
        se <-  sqrt(mtmp$var[1,1])
        df <- data.frame(term = names(terms)[i],
                         HR = exp(est),
                         HR.l = exp(est - 1.96*se),
                         HR.h = exp(est + 1.96*se))
        P <- if(is.null(P)) df else rbind(P, df)
    }
    rownames(P) <- NULL
    ## determine decr-parameters
    if(is.null(decr.inc) & is.null(decr.exc)){
        decr.inc <- UNADJ$HR - ADJ$HR > 0
        decr.exc <- !decr.inc
    }
    ## get HR from model of sequential inclusion
    Q <- UNADJ
    if(is.null(decr.inc) & Q$HR >= 1){
        decr.inc <- FALSE
    } else {
        decr.inc <- TRUE
    }
    f <- Sf
    curr.terms <- terms
    for(dummy in seq_along(terms)){ ## dummy = 1
        X <- NULL
        for(i in seq_along(curr.terms)){ ## i =1
            mtmp <- modelFNC(formula(upd(f, curr.terms[i], "+")),
                                    data = data)
            est <- mtmp$coefficients[1]
            se <-  sqrt(mtmp$var[1,1])
            df <- data.frame(term = names(curr.terms)[i],
                             HR = exp(est),
                             HR.l = exp(est - 1.96*se),
                             HR.h = exp(est + 1.96*se))
            X <- if(is.null(X)) df else rbind(X, df)
        }
        indx <- order(X$HR, decreasing = decr.inc)
        X <- X[indx, ]
        val <- curr.terms[indx][1]
        Q <- if(is.null(Q)) X[1,] else rbind(Q, X[1,])
        curr.terms <- setNames(setdiff(curr.terms, val),
                               nm = setdiff(names(curr.terms), names(val)))
        f <- upd(f, val)
    }
    rownames(Q) <- NULL
    ## get HR from model of sequential exclusion
    R <- ADJ
    if(is.null(decr.exc) & R$HR >= 1){
        decr.exc <- FALSE
    } else {
        decr.exc <- TRUE
    }
    f <- Lf
    curr.terms <- terms
    for(dummy in seq_along(terms)){ ## dummy = 1
        X <- NULL
        for(i in seq_along(curr.terms)){ ## i =1
            mtmp <- modelFNC(formula(upd(f, curr.terms[i], "-")),
                                    data = data)
            est <- mtmp$coefficients[1]
            se <-  sqrt(mtmp$var[1,1])
            df <- data.frame(term = names(curr.terms)[i],
                             HR = exp(est),
                             HR.l = exp(est - 1.96*se),
                             HR.h = exp(est + 1.96*se))
            X <- if(is.null(X)) df else rbind(X, df)
        }
        indx <- order(X$HR, decreasing = decr.exc)
        X <- X[indx, ]
        val <- curr.terms[indx][1]
        R <- if(is.null(R)) X[1,] else rbind(R, X[1,])
        curr.terms <- setNames(setdiff(curr.terms, val),
                               nm = setdiff(names(curr.terms), names(val)))
        f <- upd(f, val)
    }
    rownames(R) <- NULL
    ## return list of results
    list(
        "univariate" = P,
        "sequential_inclusion" = Q,
        "sequential_exclusion" = R
    )
}
