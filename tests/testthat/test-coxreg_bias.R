test_that("'coxreg_bias' works", {
    n <- 1000
    set.seed(20190829)
    data <- data.frame(
        id = paste0("id", 1001:(1000 + n)),
        r1 = round(rnorm(n, 20, 5)),
        r2 = round(rexp(n, 1/20)),
        c1 = sample(letters[1:5], size = n, replace = TRUE),
        c2 = factor(sample(LETTERS[5:3], size = n, replace = TRUE),
                    levels = LETTERS[6:3]),
        b1 = sample(LETTERS[6:7], size = n, replace = TRUE, prob = 2:3),
        b2 = rbinom(n, 1, 0.1),
        b3 = sample(c("No", "Yes"), size = n, replace = TRUE, prob = 1:2),
        b4 = sample(c(TRUE, FALSE), size = n, replace = TRUE),
        b5 = factor("one-level", levels = c("none-of-these", "one-level")),
        d1 = as.Date("2000-01-01") + rpois(n, 365),
        d2 = as.Date(floor(rexp(n, 1/3650)), origin = "1975-01-01"),
        s1 = survival::Surv(time = rnorm(n, 50, 7), event= rbinom(n, 1, 0.1)),
        s2 = survival::Surv(time = rexp(n, 1/40), event = rbinom(n, 1, 0.2)),
        stringsAsFactors = FALSE
    )

    coxreg_bias(data, surv = "s1", main = "b1", bnry = c("b2", "b3"),
                real = c("r1"), xtra.adj = "c1",
                bnry.manual = list("manbnry" = c(1.2, 0.1, 0.2)),
                real.manual = list("manreal" = c(1.05, 10, 12)),
                rms = FALSE)

    coxreg_bias(data, surv = "s1", main = "b1", bnry = c("b2", "b3"),
                real = c("r1"), xtra.adj = "c1",
                bnry.manual = list("manbnry" = c(1.2, 0.1, 0.2)),
                real.manual = list("manreal" = c(1.05, 10, 12)),
                rms = TRUE)


    surv <- "s1"
    main <- "b1"
    bnry <- c("b2", "b3")
    real <- "r1"
    xtra.adj = "c1"
    bnry.manual <- list("manbnry" = c(1.2, 0.1, 0.2))
    real.manual <- list("manreal" = c(1.05, 10, 12))
    rms = FALSE

})
