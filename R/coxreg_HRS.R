##' treatment HR changes
##'
##' examine  changes in  main effect  (HR, Cox  regression) when  covariates are
##'     added one at a time; in solation, sequentially increasing (starting from
##'     a small model), or sequentially decreasing (starting from large model)
##' @param data data frame
##' @param surv name of response 'Surv'-variable in data set
##' @param main name of main effect (binary)
##' @param terms vector of terms to adjust for (can be named)
##' @param uni logical; include each term's 'univariate' effect on main?
##' @param inc logical; include each  term's effect, by sequential inclusion, on
##'     main?
##' @param exc logical; include each  term's effect, by sequential exclusion, on
##'     main?
##' @param decr.inc logical; decreasing order for sequential inclusion?
##' @param decr.exc logical; decreasing order for sequential exclusion?
##' @param rms use \code{rms::cph} instead of \code{survival::coxph}? (mainly if
##'     xtra.adj contains something specific to the rms package)
##' @return a list of data frames
##' @export
coxreg_HRS <- function(data, surv, main, terms,
                       uni = TRUE, inc = TRUE, exc = TRUE,
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
    ffnc <- function(ts = NULL){
        paste0(surv, " ~ ", main,
               if(length(ts) != 0) " + " else NULL,
               paste0(ts, collapse = " + "))
    }
    Sf <- ffnc()
    Lf <- ffnc(terms)
    ## Sf <- paste0(surv, " ~ ", main) ## REMOVE
    ## Lf <- paste0(surv, " ~ ", main, " + ",
    ##              paste0(terms, collapse = " + "))
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
    if(uni){
        UNI <- UNADJ
        for(i in seq_along(terms)){ ## i = 1
            mtmp <- modelFNC(formula(upd(Sf, terms[i])), data = data)
            est <- mtmp$coefficients[1]
            se <-  sqrt(mtmp$var[1,1])
            df <- data.frame(term = names(terms)[i],
                             HR = exp(est),
                             HR.l = exp(est - 1.96*se),
                             HR.h = exp(est + 1.96*se))
            UNI <- if(is.null(UNI)) df else rbind(UNI, df)
        }
        rownames(UNI) <- NULL
    }
    ## determine decr-parameters
    if(is.null(decr.inc) & is.null(decr.exc)){
        decr.inc <- UNADJ$HR < ADJ$HR
        decr.exc <- !decr.inc
    }
    if(inc){
        ## get HR from model of sequential inclusion
        INC <- UNADJ
        if(is.null(decr.inc)) if(INC$HR >= 1){
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
            INC <- if(is.null(INC)) X[1,] else rbind(INC, X[1,])
            curr.terms <- setNames(setdiff(curr.terms, val),
                                   nm = setdiff(names(curr.terms), names(val)))
            f <- upd(f, val)
    }
    rownames(INC) <- NULL
    }
    if(exc){
        ## get HR from model of sequential exclusion
        EXC <- ADJ
        if(is.null(decr.exc)) if(EXC$HR >= 1){
                                  decr.exc <- FALSE
                              } else {
                                  decr.exc <- TRUE
                              }
        f <- Lf
        curr.terms <- terms
        for(dummy in seq_along(terms)){ ## dummy = 1
            X <- NULL
            for(i in seq_along(curr.terms)){ ## i =1
                tmp.terms <- setdiff(curr.terms, curr.terms[i])
                mtmp <- modelFNC(formula(ffnc(tmp.terms)),
                                 data = data)
                ## mtmp <- modelFNC(formula(upd(f, curr.terms[i], "-")),
                ##                  data = data) ## REMOVE
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
            EXC <- if(is.null(EXC)) X[1,] else rbind(EXC, X[1,])
            curr.terms <- setNames(setdiff(curr.terms, val),
                                   nm = setdiff(names(curr.terms), names(val)))
            ## f <- upd(f, val) ## REMOVE
        }
        rownames(EXC) <- NULL
    }
    ## return list of results
    ## list(
    ##     "univariate" = UNI,
    ##     "sequential_inclusion" = INC,
    ##     "sequential_exclusion" = EXC
    ## )
    L <- as.list(NULL)
    if(uni) L$univariate <- UNI
    if(inc) L$sequential_inclusion <- INC
    if(exc) L$sequential_exclusion <- EXC
    L
}

if(FALSE){

    ## test data for coxreg_HRS
    TMP = readRDS("ignore/hrs-list.Rds")
    data = TMP$data
    main = TMP$main
    surv = TMP$survs[1]
    terms = TMP$terms
    decr.inc = NULL
    decr.exc = NULL
    rms = FALSE

    coxreg_HRS(data = data, surv = surv, main = main, terms = terms,
               uni = TRUE, inc = TRUE, exc = TRUE,
               decr.inc = NULL, decr.exc = NULL, rms = FALSE)

    coxreg_HRS(data = data, surv = surv, main = main, terms = terms,
               uni = TRUE, inc = FALSE, exc = FALSE,
               decr.inc = NULL, decr.exc = NULL, rms = FALSE)

    coxreg_HRS(data = data, surv = surv, main = main, terms = terms,
               uni = FALSE, inc = TRUE, exc = FALSE,
               decr.inc = NULL, decr.exc = NULL, rms = FALSE)

    coxreg_HRS(data = data, surv = surv, main = main, terms = terms,
               uni = FALSE, inc = FALSE, exc = TRUE,
               decr.inc = NULL, decr.exc = NULL, rms = FALSE)

    coxreg_HRS(data = data, surv = surv, main = main, terms = terms,
               uni = FALSE, inc = FALSE, exc = FALSE,
               decr.inc = NULL, decr.exc = NULL, rms = FALSE)


}
