##' bias analysis for cox regression
##'
##' use obsSens package to get estimates  of 'true' effect given the addition of
##'     a confounder
##' @param data data frame
##' @param surv name of response 'Surv'-variable in data set
##' @param main name of main effect (binary)
##' @param bnry names of binary variables to include
##' @param real names of continuous 'normalish' variables to include
##' @param xtra.adj extra string to include in formula
##' @param bnry.manual additional binary  confounders not in model,  named list
##'     where each element is HR, proportion  at \code{main == 0}, proportion at
##'     \code{main == 1}
##' @param real.manual additional real  confounders not  in model,  named list
##'     where each element  is HR, mean at \code{main ==  0}, mean at \code{main
##'     == 1}
##' @param rms use \code{rms::cph} instead of \code{survival::coxph}? (mainly if
##'     xtra.adj contains something specific to the rms package)
##' @param verbose logical; send a message if deemed necessary?
##' @export
coxreg_bias <- function(data, surv, main,
                        bnry = NULL, real = NULL, xtra.adj = NULL,
                        bnry.manual = NULL, real.manual = NULL,
                        rms = FALSE, verbose = TRUE){
    ## must have something to work with
    if(is.null(bnry) & is.null(real) &
       is.null(bnry.manual) & is.null(real.manual)){
        stop("to much null")
    }
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
    ## limit data to what is needed ONLY if xtra.adj is NULL
    D <- if(is.null(xtra.adj)){
             subset(data, subset = TRUE,
                    select = c(surv.name, main, bnry, real))
         } else data
    ## force binary data to be 0/1
    bnry01 <- function(x){
        if(length(unique(x[!is.na(x)])) != 2) stop("some 'bnry' not binary")
        as.numeric(as.factor(x))-1
    }
    D[, c(main, bnry)] <- lapply(D[, c(main, bnry)], bnry01)
    ## guide to variables
    g0 <- data.frame(term = c(main, bnry, real, names(bnry.manual), names(real.manual)),
                     type = c("main",
                              rep("bnry", length(bnry)),
                              rep("real", length(real)),
                              rep("bnry", length(bnry.manual)),
                              rep("real", length(real.manual))),
                     manual = rep(c(0,1), c(length(c(main, bnry, real)),
                                            length(c(bnry.manual, real.manual)))),
                     stringsAsFactors = FALSE)
    ## model formula
    ftxt <- paste0(surv, " ~ ", main,
                   if(!is.null(bnry)) " + " else NULL,
                   paste(bnry, collapse = " + "),
                   if(!is.null(real)) " + " else NULL,
                   paste(real, collapse = " + "),
                   if(!is.null(xtra.adj)) " + " else NULL,
                   xtra.adj)
    if(!rms){
        M <- survival::coxph(formula(ftxt), data = D)
        sm <- survival:::summary.coxph(M)
        mod <- data.frame(term = dimnames(sm$conf.int)[[1]],
                        adjHR = sm$conf.int[, "exp(coef)"],
                        adjHR.l = sm$conf.int[, "lower .95"],
                        adjHR.u = sm$conf.int[, "upper .95"])
        rownames(mod) <- NULL
    } else {
        tmp <- D
        tmp[] <- lapply(D, function(x) if(class(x) == "Surv") NULL else x)
        if(verbose){
            .dd <<- rms::datadist(tmp)
            options(datadist = ".dd")
            mess <- paste0("\n[coxreg_bias] an effect of 'rms = TRUE' ",
                           "is that an object '.dd' is created in the ",
                           "global workspace and options(datadist='.dd') ",
                           "has been executed. (Hide this message by ",
                           "setting verbose = FALSE.)\n")
            message(mess)
        }
        M <- rms::cph(formula(ftxt), data = D, x = TRUE, y = TRUE)
        s <- rms:::summary.rms(M)
        dn <- dimnames(s)[[1]]
        indx1 <- dn != " Hazard Ratio"
        indx2 <- dimnames(s)[[2]] %in% c("Effect", "Lower 0.95", "Upper 0.95")
        d <- as.data.frame(s[!indx1, indx2])
        names(d) <- c("adjHR", "adjHR.l", "adjHR.u")
        rownames(d) <- NULL
        d$term <- dn[indx1]
        mod <- d[, c("term", "adjHR", "adjHR.l", "adjHR.u")]
    }
    ## get input stats from variables in data set
    if(!is.null(c(bnry, real))){
        br_data <- D[, c(bnry, real)]
        lmean <- function(X) unlist(lapply(X, mean, na.rm = TRUE))
        tmp <- split(br_data, f = D[[main]])
        d <- as.data.frame(lapply(tmp, lmean), check.names = FALSE)
        names(d) <- paste0('stat', names(d))
        d$term <- rownames(d)
        stat <- d[, c('term', 'stat0', 'stat1')]
    } else stat <- NULL
    ## get input stats from manually added variables
    manual <- c(bnry.manual, real.manual)
    if(!is.null(manual)){
        man.stat <- data.frame(term = names(manual),
                               stat0 = unlist(lapply(manual, function(x) x[2])),
                               stat1 = unlist(lapply(manual, function(x) x[3])))
        man.mod <- data.frame(term = names(manual),
                              adjHR = unlist(lapply(manual, function(x) x[1])),
                              adjHR.l = NA,
                              adjHR.u = NA)
    } else{
        man.stat <- man.mod <- NULL
    }
    R <- merge(g0, merge(rbind(stat, man.stat), rbind(mod, man.mod),
                         by = "term", all = TRUE), by = "term", all = TRUE)
    ## determine the changed effect of main when added confounder which is
    ## similar to the already existing covariates in distribution and HR
    R$mainHRinv.u <- R$mainHRinv.l <- R$mainHRinv <-
        R$mainHR.u <- R$mainHR.l <- R$mainHR <- rep(NA, nrow(R))
    for(i in 1:nrow(R)){ ## i = 4
        if(is.na(R$type[i]) | !R$type[i] %in% c("bnry", "real")) next
        if(R$type[i] == "bnry"){
            tmp <- obsSens::obsSensSCC(M, which = 1,
                                       g0 = R$adjHR[i],
                                       p0 = R$stat0[i], p1 = R$stat1[i],
                                       logHaz = FALSE)
            R$mainHR[i] <- tmp$beta[,,1]
            R$mainHR.l[i] <- tmp$lcl[,,1]
            R$mainHR.u[i] <- tmp$ucl[,,1]
            tmp <- obsSens::obsSensSCC(M, which = 1,
                                       g0 = 1/R$adjHR[i],
                                       p0 = R$stat0[i], p1 = R$stat1[i],
                                       logHaz = FALSE)
            R$mainHRinv[i] <- tmp$beta[,,1]
            R$mainHRinv.l[i] <- tmp$lcl[,,1]
            R$mainHRinv.u[i] <- tmp$ucl[,,1]
        }
        if(R$type[i] == "real"){
            tmp <- obsSens::obsSensSCN(M, which = 1,
                                       gamma = log(R$adjHR[i]),
                                       delta = R$stat1[i] - R$stat0[i],
                                       logHaz = FALSE)
            R$mainHR[i] <- tmp$beta[,1]
            R$mainHR.l[i] <- tmp$lcl[,1]
            R$mainHR.u[i] <- tmp$ucl[,1]
            tmp <- obsSens::obsSensSCN(M, which = 1,
                                       gamma = log(1/R$adjHR[i]),
                                       delta = R$stat1[i] - R$stat0[i],
                                       logHaz = FALSE)
            R$mainHRinv[i] <- tmp$beta[,1]
            R$mainHRinv.l[i] <- tmp$lcl[,1]
            R$mainHRinv.u[i] <- tmp$ucl[,1]
        }
    }
    attr(R, "formula") <- ftxt
    R
}

##-##' coxreg bias plot
##-##'
##-##'
##-##' @param x
##-##' @param incl.main
##-##' @param true
##-##' @return
##-##' @export
## XK include this?? is is an ass-load of dependencies
## cb_plot <- function(x, incl.main = FALSE, true = FALSE){
##     A <- filter(x, type == "main") %>%
##         transmute(term, alt = term, eff = 'Effect as-is',
##                   HR = adjHR, ci1 = adjHR.l, ci2 = adjHR.u)
##     B <- filter(x, type %in% c("bnry", "real")) %>%
##         transmute(term, alt = paste0("U as ", term), eff =  'Effect as-is',
##                   HR = mainHR, ci1 = mainHR.l, ci2 = mainHR.u)
##     C <- filter(x, type %in% c("bnry", "real")) %>%
##         transmute(term, alt = paste0("U as ", term), eff = "Inverse effect",
##                   HR = mainHRinv, ci1 = mainHRinv.l, ci2 = mainHRinv.u)
##     X <- if(incl.main) rbind(A, B, C) else rbind(B, C)
##     if(true){
##         ggplot(X, aes(x = HR, y = 1, xmin = ci1, xmax = ci2)) +
##             geom_vline(xintercept = 1, lty = 2) +
##             geom_errorbarh(height = 0.4) +
##             geom_point()+
##             scale_x_log10() +
##             scale_y_continuous(limits = c(0,2), breaks = NULL) +
##             labs(x = "Hazard ratio", y = NULL) +
##             facet_grid(alt ~ eff)
##     } else {
##         X$alt <- factor(X$alt, levels = rev(unique(X$alt)))
##         ggplot(X, aes(x = HR, y = alt, xmin = ci1, xmax = ci2)) +
##             geom_vline(xintercept = 1, lty = 2) +
##             geom_errorbarh(height = 0.4) +
##             geom_point()+
##             scale_x_log10() +
##             labs(x = "Hazard ratio", y = NULL) +
##             facet_wrap( ~ eff)
##     }
## }
