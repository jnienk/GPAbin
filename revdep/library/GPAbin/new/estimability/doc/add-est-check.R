## ----echo = FALSE, results = "hide", message = FALSE--------------------------
require("estimability")
knitr::opts_chunk$set(fig.width = 4.5, class.output = "ro")

## -----------------------------------------------------------------------------
mylm <- function(formula, data) {
    y <- data[[all.vars(formula)[1]]]
    X <- model.matrix(formula[-2], data = data)
    ch <- chol(t(X) %*% X, pivot = TRUE) |> suppressWarnings()
    rank = attr(ch, "rank")
    pivot <- attr(ch, "pivot")
    XpXinv <- chol2inv(ch, size = rank)
    coef <- rep(NA, ncol(X))
    names(coef) <- colnames(X)
    coef[pivot[1:rank]] <- XpXinv %*% (t(X[, pivot[1:rank], drop = FALSE]) %*% y)
    nonNA <- which(!is.na(coef))
    fit <- X[ , nonNA, drop = FALSE] %*% coef[nonNA]
    mse <- sum((y - fit)^2) / (length(y) - rank)
    ord <- order(pivot[1:rank])
    vcov <- mse * XpXinv[ord, ord, drop = FALSE]
    dimnames(vcov) <- list(names(coef)[nonNA], names(coef)[nonNA])
    nbasis <- estimability::nonest.basis(ch)
    terms <- terms(formula) |> delete.response()
    xlev = list()
    for (v in all.vars(terms))
        if(!is.null(lv <- levels(data[[v]])))
            xlev[[v]] <- lv
    obj <- list(terms = terms, xlev = xlev, coef = coef, vcov = vcov, nbasis = nbasis)
    class(obj) <- "mylm"
    obj
}

## -----------------------------------------------------------------------------
predict.mylm <- function(mod, newdata) {
    mf <- model.frame(mod$terms, data = newdata, xlev = mod$xlev)
    X <- model.matrix(mod$terms, data = mf)
    pred <- se <- rep(NA, nrow(X))
    estble <- estimability::is.estble(X, mod$nbasis)
    XX <- X[estble, !is.na(mod$coef), drop = FALSE]
    pred[estble] <- XX %*% mod$coef[!is.na(mod$coef)]
    se[estble] <- sqrt(diag(XX %*% mod$vcov %*% t(XX)))
    list(pred = pred, se = se)
}

## -----------------------------------------------------------------------------
# test code
warp = warpbreaks[11:40, ]
warp.mylm = mylm(breaks ~ wool*tension, warp)

new = do.call(expand.grid, warp.mylm$xlev)
cbind(new, predict(warp.mylm, newdata = new))

