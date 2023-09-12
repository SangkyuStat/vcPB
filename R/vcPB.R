#' Varying-Coefficient Disparity Decomposition Analysis for a Longitudinal Data
#'
#' The \code{vc.pb} offers Peters-Belson(PB) type of nonparametric varying-coefficient regression method which measures the disparity between a majority group
#' and a minority group for the longitudinal data.
#' @param formula a formula for the model.
#' @param group a vector within the \code{data} which is used for separating majority and minority groups.
#' @param data a data frame and data has to be included with the form of \code{data.frame}.
#' @param local_time (optional) a vector used for the local points of time variable in the kernel regression.
#' @param modifier (optional) a vector from the \code{data} which is an optional argument to add the varying term into the model. The default is \code{NULL}. If the class of the vector is given as
#' \code{integer} then, the continuous version of \code{vc.PB} is performed and if the class is \code{factor} or \code{character}, then the discrete version is proceeded. Three different sets of
#' inputs are needed for different versions.
#' @param bandwidth_M (optional) a bandwidth for the time variable used for estimating the time-varying coefficient of the majority group.
#' @param bandwidth_m (optional) a bandwidth for the time variable used for estimating the time-varying coefficient of the minority group.
#' @param bandwidth_xM (optional) a vector of \code{p} number of bandwidths for estimating the local expectations of the design matrix for the majority group.
#' @param bandwidth_xm (optional) a vector of \code{p} number of bandwidths for estimating the local expectations of the design matrix for the minority group.
#' @param bandwidth_Z_M (optional) a bandwidth for the varying variable used for estimating the time-varying coefficient of the majority group. Used only when the class of \code{modifier} is \code{integer}.
#' @param bandwidth_Z_m (optional) a bandwidth for the varying variable used for estimating the time-varying coefficient of the minority group. Used only when the class of \code{modifier} is \code{integer}.
#' @param bandwidth_Z_xM (optional) a vector of \code{p} number of bandwidths for estimating the local expectations of the design matrix related to varying variable for the majority group. Used only when the class of \code{modifier} is \code{integer}.
#' @param bandwidth_Z_xm (optional) a vector of \code{p} number of bandwidths for estimating the local expectations of the design matrix related to varying variable for the minority group. Used only when the class of \code{modifier} is \code{integer}.
#' @param detail a bool argument whether the detailed results are provided or not.
#' @param ... used for controlling the others.
#' @author Sang Kyu Lee
#' @return \code{vc.pb} returns an object of class \code{"vc.pb"}, which is a list containing
#' following components:
#' @return
#' \item{call}{a matched call.}
#' \item{overall_disparity}{overall disparity between major and minor groups.}
#' \item{explained_disparity}{explained disparity between major and minor groups, this component is given only when \code{varying} is null.}
#' \item{explained_disparity_by_X}{explained disparity from the variables without \code{modifier} variable given that the modifier variable is from the majority group, this component is given only when \code{varying} is not null.}
#' \item{explained_disparity_by_Z}{explained disparity from \code{modifier} variable, this component is given only when \code{varying} is not null.}
#' \item{unexplained_disparity}{unexplained disparity between major and minor groups.}
#' \item{times}{local time points used for kernel regression.}
#' \item{major}{a majority group label.}
#' \item{minor}{a minority group label.}
#' \item{modfier, varying.type}{the modifier variable and the type of the modifier variable, these components are given only when \code{varying} is not null.}
#' \item{bandwidths}{various corresponding bandwidths. Please see the details or the attached reference for more information.}
#' @importFrom rlist list.append
#' @importFrom KernSmooth dpill
#' @importFrom stats complete.cases dnorm glm lm model.extract model.matrix model.response na.pass model.frame quasibinomial
#' @export
vc.pb <-function(formula, group, data, modifier = NULL,
                 local_time = NULL,
                 bandwidth_M = NULL, bandwidth_m = NULL,
                 bandwidth_xM = NULL, bandwidth_xm = NULL,
                 bandwidth_Z_M = NULL, bandwidth_Z_m = NULL,
                 bandwidth_Z_xM = NULL, bandwidth_Z_xm = NULL,
                 detail = FALSE, ...)
{
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)

  m1 <- match(c("formula", "group", "data"), names(mf), 0)
  m2 <- match(c("formula", "data"), names(mf), 0)
  mf1 <- mf[c(1, m1)]
  mf1$drop.unused.levels <- TRUE
  mf1$na.action <- na.pass
  mf1[[1]] <- quote(model.frame)
  mf1$formula <- lme4::subbars(mf1$formula)

  mf1 <- eval(mf1, parent.frame())
  mt1 <- attr(mf1, "terms")

  mod <- getvc(formula)
  vf <- mod[[1]]
  tf <- mod[[2]]

  mf2 <- mf[c(1, m2)]
  mf2$drop.unused.levels <- TRUE
  mf2$na.action <- na.pass
  mf2[[1]] <- quote(model.frame)
  mf2$formula <- formula(paste("~", paste0(vf, collapse = "+")))
  mf2 <- eval(mf2, parent.frame())
  vayingX <- names(mf2)

  mf3 <- mf[c(1, m2)]
  mf3$drop.unused.levels <- TRUE
  mf3$na.action <- na.pass
  mf3[[1]] <- quote(model.frame)
  mf3$formula <- formula(paste("~", tf))
  mf3 <- eval(mf3, parent.frame())

  group <- model.extract(mf1, "group")
  if(is.null(group)) stop("group has to be defined properly.")
  disparity.var <- colnames(group)
  disparity.group <- levels(as.factor(model.extract(mf1, "group")))
  major <- disparity.group[1]
  minor <- disparity.group[2]

  time <- mf3

  yM <- model.response(mf1, "numeric")[group == disparity.group[1]]
  xM <- model.frame(mt1, mf1, contrasts.arg = NULL, xlev = NULL)[group == disparity.group[1],]
  timeM <- time[group == disparity.group[1],]
  ym <- model.response(mf1, "numeric")[group == disparity.group[2]]
  xm <- model.frame(mt1, mf1, contrasts.arg = NULL, xlev = NULL)[group == disparity.group[2],]
  timem <- time[group == disparity.group[2],]

  xM <- xM[,!(names(xM) %in% names(time))]
  xm <- xm[,!(names(xm) %in% names(time))]

  completeM <- complete.cases(xM) & complete.cases(yM) & complete.cases(timeM)
  completem <- complete.cases(xm) & complete.cases(ym) & complete.cases(timem)

  if(sum(!completeM) > 0 | sum(!completem) > 0){
    xM <- xM[completeM, ]
    xm <- xm[completem, ]
    yM <- yM[completeM]
    ym <- ym[completem]
    timeM <- timeM[completeM]
    timem <- timem[completem]
    warning("The data is not complete, the missing observations are ignored.")
  }

  xM <- xM[,-1]
  xm <- xm[,-1]

  if(is.null(local_time)) local_time <- seq(min(c(timeM, timem)), max(c(timeM, timem)),
                                            length.out = 100)

  if(is.null(modifier)){
    if(!is.null(bandwidth_Z_M)) warning("bandwidth_Z_M is ignored, since it is not required in the model.")
    if(!is.null(bandwidth_Z_m)) warning("bandwidth_Z_m is ignored, since it is not required in the model.")
    if(!is.null(bandwidth_Z_xM)) warning("bandwidth_Z_xM is ignored, since it is not required in the model.")
    if(!is.null(bandwidth_Z_xm)) warning("bandwidth_Z_xm is ignored, since it is not required in the model.")

    fitted <- time.disparity(yM = yM, xM = xM,
                             ym = ym, xm = xm,
                             time_M = timeM, time_m = timem,
                             qx = local_time,
                             varying_X = vayingX,
                             bandwidth_m = bandwidth_m,
                             bandwidth_M = bandwidth_M,
                             bandwidth_xm = bandwidth_xm,
                             bandwidth_xM = bandwidth_xM,
                             detail = F)
  } else if (!is.null(modifier)){
    # changing
    modm <- unlist(as.vector(xm[modifier]))
    modM <- unlist(as.vector(xM[modifier]))

    if((is.character(modm) & is.character(modM)) | (is.factor(modm) & is.factor(modM))){
      if(!is.null(bandwidth_Z_M)) warning("bandwidth_Z_M is ignored, since it is not required in the model.")
      if(!is.null(bandwidth_Z_m)) warning("bandwidth_Z_m is ignored, since it is not required in the model.")
      if(!is.null(bandwidth_Z_xM)) warning("bandwidth_Z_xM is ignored, since it is not required in the model.")
      if(!is.null(bandwidth_Z_xm)) warning("bandwidth_Z_xm is ignored, since it is not required in the model.")

      fitted <- time.disparity.varying.discrete(yM = yM, xM = xM,
                                                ym = ym, xm = xm,
                                                time_M = timeM, time_m = timem,
                                                qx = local_time,
                                                bandwidth_m = bandwidth_m,
                                                bandwidth_M = bandwidth_M,
                                                bandwidth_xm = bandwidth_xm,
                                                bandwidth_xM = bandwidth_xM,
                                                varying_X = vayingX,
                                                varying = modifier)
    } else if((is.integer(modm) & is.integer(modM)) | (is.numeric(modm) & is.numeric(modM))){
      fitted <- time.disparity.varying.continuous(yM = yM, xM = xM,
                                                  ym = ym, xm = xm,
                                                  time_M = timeM, time_m = timem,
                                                  qx = local_time,
                                                  varying_X = vayingX,
                                                  bandwidth1_m = bandwidth_m,
                                                  bandwidth1_M = bandwidth_M,
                                                  bandwidth1_xm = bandwidth_xm,
                                                  bandwidth1_xM = bandwidth_xM,
                                                  bandwidth2_m = bandwidth_Z_m,
                                                  bandwidth2_M = bandwidth_Z_M,
                                                  bandwidth2_xm = bandwidth_Z_xm,
                                                  bandwidth2_xM = bandwidth_Z_xM,
                                                  varying = modifier)
    } else {
      stop("The class of varying should be defined as the one of factor, character, integer and numeric.")
    }
  } else {
    stop("modifier should be defined properly.")
  }

  if(is.null(modifier)){
    res <- list(
      call = cl,
      overall_disparity = drop(fitted$result$All_Disparity),
      explained_disparity = drop(fitted$result$Explained),
      unexplained_disparity = drop(fitted$result$Unexplained),
      times = local_time,
      major = major,
      minor = minor,
      bandwidth_xM = drop(fitted$result$bandwidth_xM),
      bandwidth_xm = drop(fitted$result$bandwidth_xm),
      bandwidth_M = drop(fitted$result$bandwidth_M),
      bandwidth_m = drop(fitted$result$bandwidth_m)
    )
    if(detail) res$results = fitted$result
  } else if((!is.null(modifier) & ((is.integer(modm) & is.integer(modM)) | (is.numeric(modm) & is.numeric(modM))))){
    res <- list(
      call = cl,
      overall_disparity = drop(fitted$result$All_Disparity),
      explained_disparity_by_Z = drop(fitted$result$Explained),
      explained_disparity_by_X = drop(fitted$result$Explained_X_given_Z),
      unexplained_disparity = drop(fitted$result$Unexplained),
      times = local_time,
      modifier = modifier,
      major = major,
      minor = minor,
      varying.type = fitted$varying.type,
      bandwidth_xM = drop(fitted$result$bandwidth1_xM_seq),
      bandwidth_xm = drop(fitted$result$bandwidth1_xm_seq),
      bandwidth_M = drop(fitted$result$bandwidth1_M),
      bandwidth_m = drop(fitted$result$bandwidth1_m),
      bandwidth_Z_xM = drop(fitted$result$bandwidth2_xM_seq),
      bandwidth_Z_xm = drop(fitted$result$bandwidth2_xm_seq),
      bandwidth_Z_M = drop(fitted$result$bandwidth2_M),
      bandwidth_Z_m = drop(fitted$result$bandwidth2_m)
    )
    if(detail) res$results = fitted$result
  } else {
    varying_combination_idx = fitted$varying_combination_idx
    res <- list(
      call = cl,
      overall_disparity = drop(fitted[[varying_combination_idx]]$All_Disparity),
      explained_disparity_by_Z = drop(fitted[[varying_combination_idx]]$Explained),
      explained_disparity_by_X = drop(fitted[[varying_combination_idx]]$Explained_X_given_Z),
      unexplained_disparity = drop(fitted[[varying_combination_idx]]$Unexplained),
      times = local_time,
      modifier = modifier,
      major = major,
      minor = minor,
      disp.name = fitted$varying_combination,
      varying_combination_idx = varying_combination_idx,
      varying.type = fitted$varying.type,
      bandwidth_xM = drop(fitted[[varying_combination_idx]]$bandwidth_xM_seq),
      bandwidth_xm = drop(fitted[[varying_combination_idx]]$bandwidth_xm_seq),
      bandwidth_M = drop(fitted[[varying_combination_idx]]$bandwidth_M),
      bandwidth_m = drop(fitted[[varying_combination_idx]]$bandwidth_m)
    )
    if(detail) res$results = fitted[[varying_combination_idx]]
  }
  class(res) = "vc.pb"
  res
}

#' @method print vc.pb
#' @importFrom stats quantile
#' @export
print.vc.pb <- function(x, digits = max(3, getOption("digits") - 3), ...){
  cat("\nCall:\n\n")
  time.quantiles <- seq(0.1, 0.9, length.out = 5)
  print(x$call)
  cat("\nGroup Variables:\n\n")
  if(is.null(x$modifier)){
    cat("Major:", x$major,"\n")
    cat("Minor:", x$minor,"\n")
  } else {
    if(x$varying.type == "discrete"){
      highfreq <- strsplit(x$disp.name, split = " ")
      cat("Major:", x$major, "(", highfreq[[1]][1], ")", "\n")
      cat("Minor:", x$minor, "(", highfreq[[1]][2], ")", "\n")
    } else {
      cat("Major:", x$major,"\n")
      cat("Minor:", x$minor,"\n")
    }
  }

  time.quantiles.value <- stats::quantile(x$times, probs = time.quantiles)
  idx <- NULL
  for(i in 1:length(time.quantiles.value)){
    idx <- c(idx, which.min(x$times < time.quantiles.value[i]))
  }

  if(!is.null(x$modifier)){
    result.mat <- matrix(0, ncol = length(idx), nrow = 4)
    colnames(result.mat) <- paste0("t = ", round(x$times[idx],  digits = digits))
    # rownames(result.mat) <- c("       Overall Disparity: ", "Explained Disparity by Z: ",
    #                           "Explained Disparity by X: ", "   Unexplained Disparity: ")
    rownames(result.mat) <- c("                  Overall Disparity: ",
                              "Explained Disparity by the modifier: ",
                              "           Explained Disparity by X: ",
                              "              Unexplained Disparity: ")
    cat("\nDisparity Summary Results:\n\n")
    # cat("\nDisparity Summary Results:\n")
    # cat("\n                         ", paste0("t = ", round(x$times[idx],  digits = digits)))
    # cat("\n                  Overall Disparity: ")
    # cat(round(x$overall_disparity[idx], digits = digits))
    # cat("\nExplained Disparity by the modifier: ")
    # cat(round(x$explained_disparity_by_Z[idx], digits = digits))
    # cat("\n           Explained Disparity by X: ")
    # cat(round(x$explained_disparity_by_X[idx], digits = digits))
    # cat("\n              Unexplained Disparity: ")
    # cat(round(x$unexplained_disparity[idx], digits = digits))
    result.mat[1,] <- round(x$overall_disparity[idx], digits = digits)
    result.mat[2,] <- round(x$explained_disparity_by_Z[idx], digits = digits)
    result.mat[3,] <- round(x$explained_disparity_by_X[idx], digits = digits)
    result.mat[4,] <- round(x$unexplained_disparity[idx], digits = digits)
  } else {
    result.mat <- matrix(0, ncol = length(idx), nrow = 3)
    colnames(result.mat) <- paste0("t = ", round(x$times[idx],  digits = digits))
    rownames(result.mat) <- c("    Overall Disparity: ",
                              "  Explained Disparity: ",
                              "Unexplained Disparity: ")
    cat("\nDisparity Summary Results:\n\n")
    # cat("\n                      ", paste0("t = ", round(x$times[idx],  digits = digits)))
    # cat("\n    Overall Disparity: ")
    # cat(round(x$overall_disparity[idx], digits = digits))
    # cat("\n  Explained Disparity: ")
    # cat(round(x$explained_disparity[idx], digits = digits))
    # cat("\nUnexplained Disparity: ")
    # cat(round(x$unexplained_disparity[idx], digits = digits), "\n\n")
    result.mat[1,] <- round(x$overall_disparity[idx], digits = digits)
    result.mat[2,] <- round(x$explained_disparity[idx], digits = digits)
    result.mat[3,] <- round(x$unexplained_disparity[idx], digits = digits)
  }
  print(result.mat)
  cat("\n")
}

#' @importFrom lme4 expandDoubleVerts
#' @importFrom methods is
sb = function(term){
  fb <- function(term) {
    if (is.name(term) || !is.language(term))
      return(NULL)
    if (term[[1]] == as.name("("))
      return(fb(term[[2]]))
    stopifnot(is.call(term))
    if (term[[1]] == as.name("|"))
      return(term)
    if (length(term) == 2)
      return(fb(term[[2]]))
    c(fb(term[[2]]), fb(term[[3]]))
  }
  modterm <- lme4::expandDoubleVerts(if (is(term, "formula"))
    term[[length(term)]]
    else term)
  strsplit(as.character(fb(modterm)), split = " | ")
}

getvc <- function(f){
  bars <- sb(f)
  res.varying <- NULL
  res.time <- bars[[1]][3]
  for(j in 1:length(bars)){
    res.varying <- c(res.varying, bars[[j]][1])
    if(bars[[j]][3] != res.time) stop("time variables have to be the same for the varying variables")
  }
  res = list(res.varying = res.varying, res.time = res.time)
  res
}

time.disparity.varying.discrete = function(yM, xM, ym, xm, varying, varying_X,
                                           time_M, time_m,
                                           qx = seq(min(c(time_M, time_m)), max(c(time_M, time_m)),
                                                    length.out = 100),
                                           bandwidth_M = NULL, bandwidth_m = NULL,
                                           bandwidth_xM = NULL, bandwidth_xm = NULL){
  varying_idx = which(names(xM) %in% varying)
  varying_lbl = levels(as.factor(xM[,names(xM) %in% varying]))
  l = length(varying_lbl)
  varying_list_M = list()
  varying_list_m = list()
  for(q in 1:l){
    varying_list_M[[q]] = cbind(yM[xM[,varying_idx] == varying_lbl[q]],
                                time_M[xM[,varying_idx] == varying_lbl[q]],
                                xM[xM[,varying_idx] == varying_lbl[q], which(!(names(xM) %in% varying))])
    varying_list_m[[q]] = cbind(ym[xm[,varying_idx] == varying_lbl[q]],
                                time_m[xm[,varying_idx] == varying_lbl[q]],
                                xm[xm[,varying_idx] == varying_lbl[q], which(!(names(xm) %in% varying))])
  }
  varying_lbl_M = varying_lbl[which.max(unlist(lapply(lapply(varying_list_M, dim), function(x) x[1])))]
  varying_lbl_m = varying_lbl[which.max(unlist(lapply(lapply(varying_list_m, dim), function(x) x[1])))]
  varying_lbl_M_idx = which.max(unlist(lapply(lapply(varying_list_M, dim), function(x) x[1])))
  varying_lbl_m_idx = which.max(unlist(lapply(lapply(varying_list_m, dim), function(x) x[1])))
  res_list = list()
  for(k in 1:l){
    xM_temp = varying_list_M[[k]][,-c(1,2)]
    yM_temp = varying_list_M[[k]][,1]
    time_M_temp = varying_list_M[[k]][,2]
    if(k == 1){
      p = ncol(xM_temp)

      charac = sapply(xM_temp, is.character)
      fac = sapply(xM_temp, is.factor)
      varying_coef_bool = names(xM_temp) %in% varying_X
      xM_model = list()
      xm_model = list()
      cat_idx = NULL
      idx = 1
      for(j in 1:p){
        # the model.matrix function does not work if there is only one level
        # in the variable, however, the estimation is nonsense for this case
        # so I don't have to worry about this
        xM_model[[j]] = model.matrix(~., data = as.data.frame(xM_temp[,j]))
        if((charac[j] == T | fac[j] == T) & (dim(xM_model[[j]])[2] == 2) & varying_coef_bool[j]){
          cat_idx = c(cat_idx, (idx + 1))
          idx = idx + 1
        } else if((charac[j] == T | fac[j] == T) & (dim(xM_model[[j]])[2] > 2) & varying_coef_bool[j]){
          cat_idx = c(cat_idx, (idx + 1:(dim(xM_model[[j]])[2]-1)))
          idx = max(idx + 1:(dim(xM_model[[j]])[2]-1))
        } else{
          idx = idx + 1
        }
      }
    }
    xM_model = model.matrix(~., data = as.data.frame(xM_temp))
    varying_coef_model_idx = which(colnames(xM_model) %in% varying_X)
    xm_temp = varying_list_m[[k]][,-c(1,2)]
    ym_temp = varying_list_m[[k]][,1]
    time_m_temp = varying_list_m[[k]][,2]
    xm_model = model.matrix(~., data = as.data.frame(xm_temp))

    pred_xM = matrix(0, nrow = length(qx), ncol = dim(xm_model)[2] - 1)
    pred_xm = matrix(0, nrow = length(qx), ncol = dim(xm_model)[2] - 1)
    bandwidth_xM_seq <- NULL
    bandwidth_xm_seq <- NULL
    for(i in 1:length(qx)){
      for(j in 2:dim(xM_model)[2]){
        if(j %in% varying_coef_model_idx){
          if(i == 1){
            bandwidth_xM_temp <- ifelse(is.null(bandwidth_xM), 3*KernSmooth::dpill(x = time_M_temp, y = xM_model[,j]), bandwidth_xM[which(varying_coef_model_idx %in% j)])
            bandwidth_xm_temp <- ifelse(is.null(bandwidth_xm), 3*KernSmooth::dpill(x = time_m_temp, y = xm_model[,j]), bandwidth_xm[which(varying_coef_model_idx %in% j)])

            bandwidth_xM_seq <- c(bandwidth_xM_temp, bandwidth_xM_seq)
            bandwidth_xm_seq <- c(bandwidth_xm_temp, bandwidth_xm_seq)
          }
          kernel_M = dnorm((time_M_temp - qx[i])/bandwidth_xM_temp)
          kernel_m = dnorm((time_m_temp - qx[i])/bandwidth_xm_temp)
        }

        if(j %in% cat_idx){
          pred_xM[i,(j-1)] = mean(glm(xM_model[,j] ~ 1, family = quasibinomial(link = "logit"), weights=kernel_M)$fitted.values)
          pred_xm[i,(j-1)] = mean(glm(xm_model[,j] ~ 1, family = quasibinomial(link = "logit"), weights=kernel_m)$fitted.values)
        } else if(j %in% varying_coef_model_idx){
          pred_xM[i,(j-1)] = mean(lm(xM_model[,j] ~ 1, weights=kernel_M)$fitted.values)
          pred_xm[i,(j-1)] = mean(lm(xm_model[,j] ~ 1, weights=kernel_m)$fitted.values)
        } else {
          pred_xM[i,(j-1)] = mean(xM_model[,j])
          pred_xm[i,(j-1)] = mean(xm_model[,j])
        }
        # if(j %in% cat_idx){
        #   pred_xM[i,(j-1)] = mean(glm(xM_model[,j] ~ 1, family = quasibinomial(link = "logit"), weights=kernel_M)$fitted.values)
        #   pred_xm[i,(j-1)] = mean(glm(xm_model[,j] ~ 1, family = quasibinomial(link = "logit"), weights=kernel_m)$fitted.values)
        # } else {
        #   pred_xM[i,(j-1)] = mean(lm(xM_model[,j] ~ 1, weights=kernel_M)$fitted.values)
        #   pred_xm[i,(j-1)] = mean(lm(xm_model[,j] ~ 1, weights=kernel_m)$fitted.values)
        # }
      }
    }
    hatcoeffM = matrix(0,length(qx), dim(xM_model)[2])
    hatcoeffm = matrix(0,length(qx), dim(xm_model)[2])

    bandwidth_M_temp <- ifelse(is.null(bandwidth_M), 3*KernSmooth::dpill(x = time_M_temp, y = as.numeric(yM_temp)), bandwidth_M)
    bandwidth_m_temp <- ifelse(is.null(bandwidth_m), 3*KernSmooth::dpill(x = time_m_temp, y = as.numeric(ym_temp)), bandwidth_m)
    for(i in 1:length(qx))
    {
      kernel_M = dnorm((time_M_temp - qx[i])/bandwidth_M_temp)
      kernel_m = dnorm((time_m_temp - qx[i])/bandwidth_m_temp)

      hatcoeffM[i,] = lm(yM ~. , data = data.frame(yM = as.numeric(yM_temp), xM_temp),
                         weights=kernel_M)$coefficient
      hatcoeffm[i,] = lm(ym ~. , data = data.frame(ym = as.numeric(ym_temp), xm_temp),
                         weights=kernel_m)$coefficient
    }
    hatcoeff1M = hatcoeffM
    hatcoeff1m = hatcoeffm
    pred_X1M = t(pred_xM)
    pred_X1m = t(pred_xm) # z^M with coefficient minor part

    for(h in 1:l){
      xm_temp = varying_list_m[[h]][,-c(1,2)]
      ym_temp = varying_list_m[[h]][,1]
      time_m_temp = varying_list_m[[h]][,2]

      xm_model = model.matrix(~., data = as.data.frame(xm_temp))

      pred_xm = matrix(0, nrow = length(qx), ncol = dim(xm_model)[2] - 1)
      for(i in 1:length(qx)){
        for(j in 2:dim(xm_model)[2]){
          if(j %in% varying_coef_model_idx){
            if(i == 1){
              bandwidth_xM_temp <- ifelse(is.null(bandwidth_xM), 3*KernSmooth::dpill(x = time_M_temp, y = xM_model[,j]), bandwidth_xM[which(varying_coef_model_idx %in% j)])
              bandwidth_xm_temp <- ifelse(is.null(bandwidth_xm), 3*KernSmooth::dpill(x = time_m_temp, y = xm_model[,j]), bandwidth_xm[which(varying_coef_model_idx %in% j)])
            }
            kernel_mm = dnorm((time_m_temp - qx[i])/ifelse(h == q, bandwidth_xM_temp, bandwidth_xm_temp))
          }

          if(j %in% cat_idx){
            pred_xm[i,(j-1)] = mean(glm(xm_model[,j] ~ 1, weights = kernel_mm, family = quasibinomial(link = "logit"))$fitted.values)
          } else if(j %in% varying_coef_model_idx){
            pred_xm[i,(j-1)] = mean(lm(xm_model[,j] ~ 1, weights = kernel_mm)$fitted.values)
          } else {
            pred_xm[i,(j-1)] = mean(xm_model[,j])
          }
        }
        # for(j in 2:dim(xm_model)[2]){
        #   if(i == 1){
        #     bandwidth_xM_temp <- ifelse(is.null(bandwidth_xM), 2*KernSmooth::dpill(x = time_M_temp, y = xM_model[,j]), bandwidth_xM[j-1])
        #     bandwidth_xm_temp <- ifelse(is.null(bandwidth_xm), 2*KernSmooth::dpill(x = time_m_temp, y = xm_model[,j]), bandwidth_xm[j-1])
        #   }
        #
        #   kernel_mm = dnorm((time_m_temp - qx[i])/ifelse(h == q, bandwidth_xM_temp, bandwidth_xm_temp))
        #   if(j %in% cat_idx){
        #     pred_xm[i,(j-1)] = mean(glm(xm_model[,j] ~ 1, weights = kernel_mm, family = quasibinomial(link = "logit"))$fitted.values)
        #   } else {
        #     pred_xm[i,(j-1)] = mean(lm(xm_model[,j] ~ 1, weights = kernel_mm)$fitted.values)
        #   }
        # }
      }
      hatcoeffm = matrix(0,length(qx), dim(xm_model)[2])

      for(i in 1:length(qx))
      {
        kernel_m = dnorm((time_m_temp - qx[i])/ifelse(h == q, bandwidth_M_temp, bandwidth_m_temp))
        hatcoeffm[i,] = lm(ym ~. , data = data.frame(ym = as.numeric(ym_temp), xm_temp),
                           weights=kernel_m)$coefficient
      }
      hatcoeff2m = hatcoeffm
      pred_X2m = t(pred_xm) # z^m with coefficient minor part

      result = list()
      # eq3 = hatcoeff1M[,1] - hatcoeff1m[,1]
      # eq4 = (hatcoeff1M[,-1] - hatcoeff1m[,-1]) %*% (colMeans(pred_X1M))
      # eq5 = hatcoeff1m[,-1]%*%(colMeans(pred_X1M) - colMeans(pred_X1m))
      # eq6 = hatcoeff1m[,1] - hatcoeff2m[,1]
      # eq7 = hatcoeff1m[,-1]%*%colMeans(pred_X1m) - hatcoeff2m[,-1]%*%colMeans(pred_X2m)
      eq3 = hatcoeff1M[,1] - hatcoeff1m[,1]
      eq4 = diag((hatcoeff1M[,-1] - hatcoeff1m[,-1]) %*% (pred_X1M))
      eq5 = diag(hatcoeff1m[,-1]%*%((pred_X1M) - (pred_X1m)))
      eq6 = hatcoeff1m[,1] - hatcoeff2m[,1]
      eq7 = diag(hatcoeff1m[,-1]%*%(pred_X1m)) - diag(hatcoeff2m[,-1]%*%(pred_X2m))

      result$Unexplained = eq3 + eq4
      result$Explained_X_given_Z = eq5
      result$Direct_Explained = eq6
      result$Indirect_Explained = eq7
      result$Explained = eq6 + eq7
      result$All_Disparity = eq3 + eq4 + eq5 + eq6 + eq7
      result$hatcoeffmm = hatcoeff1m
      result$hatcoeffMM = hatcoeff1M
      result$hatcoeffmM = hatcoeff2m
      result$time = qx
      result$varying_coef_model_idx = varying_coef_model_idx
      result$varying_coef_bool = varying_coef_bool
      result$bandwidth_xM_seq = bandwidth_xM_seq
      result$bandwidth_xm_seq = bandwidth_xm_seq
      result$bandwidth_M = bandwidth_M_temp
      result$bandwidth_m = bandwidth_m_temp

      res_list = list.append(res_list, result)
    }
  }
  obj.names = as.vector(t(outer(varying_lbl, varying_lbl, FUN = paste)))
  names(res_list) = obj.names
  res_list$yM = yM
  res_list$ym = ym
  res_list$xM = xM
  res_list$xm = xm
  res_list$varying = varying
  res_list$varying_combination = paste(varying_lbl_M, varying_lbl_m)
  res_list$varying_combination_idx = (varying_lbl_M_idx-1)*l + varying_lbl_m_idx
  res_list$time_M = time_M
  res_list$time_m = time_m
  res_list$qx = qx
  res_list$bandwidth_M = bandwidth_M
  res_list$bandwidth_m = bandwidth_m
  res_list$bandwidth_xM = bandwidth_xM
  res_list$bandwidth_xm = bandwidth_xm
  res_list$obj.names = obj.names
  res_list$varying.type = "discrete"
  res_list
}

time.disparity.varying.continuous = function(yM, xM, ym, xm, varying, varying_X,
                                             time_M, time_m,
                                             qx = seq(min(c(time_M, time_m)), max(c(time_M, time_m)),
                                                      length.out = 100),
                                             bandwidth1_m = NULL,
                                             bandwidth1_M = NULL,
                                             bandwidth1_xm = NULL,
                                             bandwidth1_xM = NULL,
                                             bandwidth2_m = NULL,
                                             bandwidth2_M = NULL,
                                             bandwidth2_xm = NULL,
                                             bandwidth2_xM = NULL,
                                             varying_M_mean = NULL,
                                             varying_m_mean = NULL){

  varying_idx = which(names(xM) %in% varying)
  varying_list_M = cbind(yM,
                         time_M,
                         xM[, which(!(names(xM) %in% varying))])
  varying_list_m = cbind(ym,
                         time_m,
                         xm[, which(!(names(xm) %in% varying))])
  res_list = list()

  xM_temp = varying_list_M[,-c(1,2)]
  yM_temp = varying_list_M[,1]
  time_M_temp = varying_list_M[,2]
  varying_M_temp = xM[, which(names(xM) %in% varying)]
  if(is.null(varying_M_mean)){
    varying_M_temp_mean = mean(varying_M_temp)
  } else {
    varying_M_temp_mean = varying_M_mean
  }
  p = ncol(xM_temp)

  charac = sapply(xM_temp, is.character)
  fac = sapply(xM_temp, is.factor)
  varying_coef_bool = names(xM_temp) %in% varying_X
  xM_model = list()
  xm_model = list()
  cat_idx = NULL
  idx = 1
  for(j in 1:p){
    # the model.matrix function does not work if there is only one level
    # in the variable, however, the estimation is nonsense for this case
    # so I don't have to worry about this
    xM_model[[j]] = model.matrix(~., data = as.data.frame(xM_temp[,j]))
    if((charac[j] == T | fac[j] == T) & (dim(xM_model[[j]])[2] == 2) & varying_coef_bool[j]){
      cat_idx = c(cat_idx, (idx + 1))
      idx = idx + 1
    } else if((charac[j] == T | fac[j] == T) & (dim(xM_model[[j]])[2] > 2) & varying_coef_bool[j]){
      cat_idx = c(cat_idx, (idx + 1:(dim(xM_model[[j]])[2]-1)))
      idx = max(idx + 1:(dim(xM_model[[j]])[2]-1))
    } else{
      idx = idx + 1
    }
  }

  xM_model = model.matrix(~., data = as.data.frame(xM_temp))
  varying_coef_model_idx = which(colnames(xM_model) %in% varying_X)
  xm_temp = varying_list_m[,-c(1,2)]
  ym_temp = varying_list_m[,1]
  time_m_temp = varying_list_m[,2]
  varying_m_temp = xm[, which(names(xm) %in% varying)]
  if(is.null(varying_m_mean)){
    varying_m_temp_mean = mean(varying_m_temp)
  } else {
    varying_m_temp_mean = varying_m_mean
  }
  xm_model = model.matrix(~., data = as.data.frame(xm_temp))

  pred_xM = matrix(0, nrow = length(qx), ncol = dim(xM_model)[2] - 1)
  pred_xm = matrix(0, nrow = length(qx), ncol = dim(xm_model)[2] - 1)

  bandwidth1_xM_seq <- NULL
  bandwidth1_xm_seq <- NULL
  bandwidth2_xM_seq <- NULL
  bandwidth2_xm_seq <- NULL

  for(i in 1:length(qx)){
    for(j in 2:dim(xM_model)[2]){
      if(j %in% varying_coef_model_idx){
        if(i == 1){
          bandwidth1_xM_temp <- ifelse(is.null(bandwidth1_xM), 3*KernSmooth::dpill(x = time_M_temp, y = xM_model[,j]), bandwidth1_xM[which(varying_coef_model_idx %in% j)])
          bandwidth1_xm_temp <- ifelse(is.null(bandwidth1_xm), 3*KernSmooth::dpill(x = time_m_temp, y = xm_model[,j]), bandwidth1_xm[which(varying_coef_model_idx %in% j)])
          bandwidth2_xM_temp <- ifelse(is.null(bandwidth2_xM), 3*KernSmooth::dpill(x = varying_M_temp, y = xM_model[,j]), bandwidth2_xM[which(varying_coef_model_idx %in% j)])
          bandwidth2_xm_temp <- ifelse(is.null(bandwidth2_xm), 3*KernSmooth::dpill(x = varying_m_temp, y = xm_model[,j]), bandwidth2_xm[which(varying_coef_model_idx %in% j)])

          bandwidth1_xM_seq <- c(bandwidth1_xM_temp, bandwidth1_xM_seq)
          bandwidth1_xm_seq <- c(bandwidth1_xm_temp, bandwidth1_xm_seq)
          bandwidth2_xM_seq <- c(bandwidth2_xM_temp, bandwidth2_xM_seq)
          bandwidth2_xm_seq <- c(bandwidth2_xm_temp, bandwidth2_xm_seq)
        }

        kernel_M = dnorm((time_M_temp - qx[i])/bandwidth1_xM_temp) * dnorm((varying_M_temp - varying_M_temp_mean)/bandwidth2_xM_temp)
        kernel_m = dnorm((time_m_temp - qx[i])/bandwidth1_xm_temp) * dnorm((varying_m_temp - varying_m_temp_mean)/bandwidth2_xm_temp)

      }

      if(j %in% cat_idx){
        pred_xM[i,(j-1)] = mean(glm(xM_model[,j] ~ 1, family = quasibinomial(link = "logit"), weights=kernel_M)$fitted.values)
        pred_xm[i,(j-1)] = mean(glm(xm_model[,j] ~ 1, family = quasibinomial(link = "logit"), weights=kernel_m)$fitted.values)
      } else if(j %in% varying_coef_model_idx){
        pred_xM[i,(j-1)] = mean(lm(xM_model[,j] ~ 1, weights=kernel_M)$fitted.values)
        pred_xm[i,(j-1)] = mean(lm(xm_model[,j] ~ 1, weights=kernel_m)$fitted.values)
      } else {
        pred_xM[i,(j-1)] = mean(xM_model[,j])
        pred_xm[i,(j-1)] = mean(xm_model[,j])
      }
    }
  }
  hatcoeffM = matrix(0,length(qx), dim(xM_model)[2])
  hatcoeffm = matrix(0,length(qx), dim(xm_model)[2])

  bandwidth1_M_temp <- ifelse(is.null(bandwidth1_M), 3*KernSmooth::dpill(x = time_M_temp, y = as.numeric(yM_temp)), bandwidth1_M)
  bandwidth1_m_temp <- ifelse(is.null(bandwidth1_m), 3*KernSmooth::dpill(x = time_m_temp, y = as.numeric(ym_temp)), bandwidth1_m)
  bandwidth2_M_temp <- ifelse(is.null(bandwidth2_M), 3*KernSmooth::dpill(x = varying_M_temp, y = as.numeric(yM_temp)), bandwidth2_M)
  bandwidth2_m_temp <- ifelse(is.null(bandwidth2_m), 3*KernSmooth::dpill(x = varying_m_temp, y = as.numeric(ym_temp)), bandwidth2_m)
  for(i in 1:length(qx))
  { # kernel part needs to be changed 3:42 / (08/21)
    kernel_M = dnorm((time_M_temp - qx[i])/bandwidth1_M_temp) * dnorm((varying_M_temp - varying_M_temp_mean)/bandwidth2_M_temp)
    kernel_m = dnorm((time_m_temp - qx[i])/bandwidth1_m_temp) * dnorm((varying_m_temp - varying_m_temp_mean)/bandwidth2_m_temp)

    hatcoeffM[i,] = lm(yM ~. , data = data.frame(yM = as.numeric(yM_temp), xM_temp),
                       weights=kernel_M)$coefficient
    hatcoeffm[i,] = lm(ym ~. , data = data.frame(ym = as.numeric(ym_temp), xm_temp),
                       weights=kernel_m)$coefficient
  }
  hatcoeff1M = hatcoeffM
  hatcoeff1m = hatcoeffm
  pred_X1M = t(pred_xM)
  pred_X1m = t(pred_xm)


  pred_xm = matrix(0, nrow = length(qx), ncol = dim(xm_model)[2] - 1)
  for(i in 1:length(qx)){
    # kernel part needs to be changed 3:42 / (08/21)
    for(j in 2:dim(xm_model)[2]){
      if(j %in% varying_coef_model_idx){
        if(i == 1){
          bandwidth1_xm_temp <- ifelse(is.null(bandwidth1_xm), 3*KernSmooth::dpill(x = time_m_temp, y = xm_model[,j]), bandwidth1_xm[which(varying_coef_model_idx %in% j)])
          bandwidth2_xm_temp <- ifelse(is.null(bandwidth2_xm), 3*KernSmooth::dpill(x = varying_m_temp, y = xm_model[,j]), bandwidth2_xm[which(varying_coef_model_idx %in% j)])
        }
        kernel_mm = dnorm((time_m_temp - qx[i])/bandwidth1_xm_temp) * dnorm((varying_m_temp - varying_M_temp_mean)/bandwidth2_xm_temp)
      }

      if(j %in% cat_idx){
        pred_xm[i,(j-1)] = mean(glm(xm_model[,j] ~ 1, weights = kernel_mm, family = quasibinomial(link = "logit"))$fitted.values)
      } else if(j %in% varying_coef_model_idx) {
        pred_xm[i,(j-1)] = mean(lm(xm_model[,j] ~ 1, weights = kernel_mm)$fitted.values)
      } else {
        pred_xm[i,(j-1)] = mean(xm_model[,j])
      }
    }
  }
  hatcoeffm = matrix(0,length(qx), dim(xm_model)[2])

  for(i in 1:length(qx))
  {
    # kernel part needs to be changed 3:42 / (08/21)
    kernel_m = dnorm((time_m_temp - qx[i])/bandwidth1_m_temp) * dnorm((varying_m_temp - varying_M_temp_mean)/bandwidth2_m_temp)
    hatcoeffm[i,] = lm(ym ~. , data = data.frame(ym = as.numeric(ym_temp), xm_temp),
                       weights=kernel_m)$coefficient
  }
  hatcoeff2m = hatcoeffm
  pred_X2m = t(pred_xm)

  result = list()
  eq3 = hatcoeff1M[,1] - hatcoeff1m[,1]
  eq4 = diag((hatcoeff1M[,-1] - hatcoeff1m[,-1]) %*% (pred_X1M))
  eq5 = diag(hatcoeff1m[,-1]%*%((pred_X1M) - (pred_X1m)))
  eq6 = hatcoeff1m[,1] - hatcoeff2m[,1]
  eq7 = diag(hatcoeff1m[,-1]%*%(pred_X1m)) - diag(hatcoeff2m[,-1]%*%(pred_X2m))

  result$Unexplained = eq3 + eq4
  result$Explained_X_given_Z = eq5
  result$Direct_Explained = eq6
  result$Indirect_Explained = eq7
  result$Explained = eq6 + eq7
  result$All_Disparity = eq3 + eq4 + eq5 + eq6 + eq7
  result$time = qx
  result$bandwidth1_xM_seq = bandwidth1_xM_seq
  result$bandwidth1_xm_seq = bandwidth1_xm_seq
  result$bandwidth1_M = bandwidth1_M_temp
  result$bandwidth1_m = bandwidth1_m_temp
  result$bandwidth2_xM_seq = bandwidth2_xM_seq
  result$bandwidth2_xm_seq = bandwidth2_xm_seq
  result$bandwidth2_M = bandwidth2_M_temp
  result$bandwidth2_m = bandwidth2_m_temp

  res_list$result = result

  res_list$yM = yM
  res_list$ym = ym
  res_list$xM = xM
  res_list$xm = xm
  res_list$varying = varying
  res_list$time_M = time_M
  res_list$time_m = time_m
  res_list$qx = qx
  res_list$varying_M_mean = varying_M_temp_mean
  res_list$varying_m_mean = varying_m_temp_mean
  res_list$varying.type = "continuous"
  res_list
}

time.disparity = function(yM, xM, ym, xm,
                          time_M, time_m,
                          varying = NULL, varying_X,
                          qx = seq(min(c(time_M, time_m)), max(c(time_M, time_m)),
                                   length.out = 100),
                          bandwidth_M = NULL, bandwidth_m = NULL,
                          bandwidth_xM = NULL, bandwidth_xm = NULL,
                          detail = T){
  varying_list_M = cbind(yM, time_M, xM)
  varying_list_m = cbind(ym, time_m, xm)
  res_list = list()

  xM_temp = varying_list_M[,-c(1,2)]
  yM_temp = varying_list_M[,1]
  time_M_temp = varying_list_M[,2]
  if(!is.null(varying)){
    varying_idx = which(names(xM) %in% varying) # not including the intercept
  }
  p = ncol(xM_temp)

  charac = sapply(xM_temp, is.character)
  fac = sapply(xM_temp, is.factor)
  varying_coef_bool = names(xM_temp) %in% varying_X
  xM_model = list()
  xm_model = list()
  cat_idx = NULL
  idx = 1
  for(j in 1:p){
    # the model.matrix function does not work if there is only one level
    # in the variable, however, the estimation is nonsense for this case
    # so I don't have to worry about this
    xM_model[[j]] = model.matrix(~., data = as.data.frame(xM_temp[,j]))
    if((charac[j] == T | fac[j] == T) & (dim(xM_model[[j]])[2] == 2) & varying_coef_bool[j]){
      if(!is.null(varying)){
        if(j == varying_idx){
          varying_idx = idx
        }
      }
      cat_idx = c(cat_idx, (idx + 1))
      idx = idx + 1
    } else if((charac[j] == T | fac[j] == T) & (dim(xM_model[[j]])[2] > 2) & varying_coef_bool[j]){
      if(!is.null(varying)){
        if(j == varying_idx){
          varying_idx = idx + 1:(dim(xM_model[[j]])[2]-1) - 1
        }
      }
      cat_idx = c(cat_idx, (idx + 1:(dim(xM_model[[j]])[2]-1)))
      idx = max(idx + 1:(dim(xM_model[[j]])[2]-1))
    } else{
      idx = idx + 1
    }
  }

  xM_model = model.matrix(~., data = as.data.frame(xM_temp))
  varying_coef_model_idx = which(colnames(xM_model) %in% varying_X)
  xm_temp = varying_list_m[,-c(1,2)]
  ym_temp = varying_list_m[,1]
  time_m_temp = varying_list_m[,2]

  xm_model = model.matrix(~., data = as.data.frame(xm_temp))

  pred_xM = matrix(0, nrow = length(qx), ncol = dim(xM_model)[2] - 1)
  pred_xm = matrix(0, nrow = length(qx), ncol = dim(xm_model)[2] - 1)

  bandwidth_xM_seq <- NULL
  bandwidth_xm_seq <- NULL
  for(i in 1:length(qx)){
    for(j in 2:dim(xM_model)[2]){
      if(j %in% varying_coef_model_idx){
        if(i == 1){
          bandwidth_xM_temp <- ifelse(is.null(bandwidth_xM), 3*KernSmooth::dpill(x = time_M_temp, y = xM_model[,j]), bandwidth_xM[which(varying_coef_model_idx %in% j)])
          bandwidth_xm_temp <- ifelse(is.null(bandwidth_xm), 3*KernSmooth::dpill(x = time_m_temp, y = xm_model[,j]), bandwidth_xm[which(varying_coef_model_idx %in% j)])

          bandwidth_xM_seq <- c(bandwidth_xM_temp, bandwidth_xM_seq)
          bandwidth_xm_seq <- c(bandwidth_xm_temp, bandwidth_xm_seq)
        }

        kernel_M = dnorm((time_M_temp - qx[i])/bandwidth_xM_temp)
        kernel_m = dnorm((time_m_temp - qx[i])/bandwidth_xm_temp)
      }

      if(j %in% cat_idx){
        pred_xM[i,(j-1)] = mean(glm(xM_model[,j] ~ 1, family = quasibinomial(link = "logit"), weights=kernel_M)$fitted.values)
        pred_xm[i,(j-1)] = mean(glm(xm_model[,j] ~ 1, family = quasibinomial(link = "logit"), weights=kernel_m)$fitted.values)
      } else if(j %in% varying_coef_model_idx) {
        pred_xM[i,(j-1)] = mean(lm(xM_model[,j] ~ 1, weights=kernel_M)$fitted.values)
        pred_xm[i,(j-1)] = mean(lm(xm_model[,j] ~ 1, weights=kernel_m)$fitted.values)
      } else {
        pred_xM[i,(j-1)] = mean(xM_model[,j])
        pred_xm[i,(j-1)] = mean(xm_model[,j])
      }
    }
  }
  hatcoeffM = matrix(0,length(qx), dim(xM_model)[2])
  hatcoeffm = matrix(0,length(qx), dim(xm_model)[2])

  bandwidth_M_temp <- ifelse(is.null(bandwidth_M), 3*KernSmooth::dpill(x = time_M_temp, y = as.numeric(yM_temp)), bandwidth_M)
  bandwidth_m_temp <- ifelse(is.null(bandwidth_m), 3*KernSmooth::dpill(x = time_m_temp, y = as.numeric(ym_temp)), bandwidth_m)
  for(i in 1:length(qx))
  {
    kernel_M = dnorm((time_M_temp - qx[i])/bandwidth_M_temp)
    kernel_m = dnorm((time_m_temp - qx[i])/bandwidth_m_temp)

    hatcoeffM[i,] = lm(yM ~. , data = data.frame(yM = as.numeric(yM_temp), xM_temp),
                       weights=kernel_M)$coefficient
    hatcoeffm[i,] = lm(ym ~. , data = data.frame(ym = as.numeric(ym_temp), xm_temp),
                       weights=kernel_m)$coefficient
  }
  hatcoeff1M = hatcoeffM
  hatcoeff1m = hatcoeffm
  pred_X1M = t(pred_xM)
  pred_X1m = t(pred_xm)

  result = list()
  eq3 = hatcoeff1M[,1] - hatcoeff1m[,1]
  eq4 = diag((hatcoeff1M[,-1] - hatcoeff1m[,-1]) %*% (pred_X1M))
  if(detail & !is.null(varying)){
    eq5 = diag(hatcoeff1m[,-c(1, varying_idx + 1)]%*%((pred_X1M[-varying_idx,]) - (pred_X1m[-varying_idx,])))

    if(length(varying_idx == 1)){
      eq6 = diag(hatcoeff1m[,(varying_idx + 1)]%*%t((pred_X1M[varying_idx,]) - (pred_X1m[varying_idx,])))
    } else if(length(varying_idx > 1)){
      eq6 = diag(hatcoeff1m[,(varying_idx + 1)]%*%((pred_X1M[varying_idx,]) - (pred_X1m[varying_idx,])))
    }
  } else {
    eq5 = diag(hatcoeff1m[,-c(1)]%*%((pred_X1M) - (pred_X1m)))
  }

  result$Unexplained = eq3 + eq4
  result$Explained = eq5
  if(detail){
    result$Explained_by_Z = eq6
    result$All_Disparity = eq3 + eq4 + eq5 + eq6
  } else {
    result$All_Disparity = eq3 + eq4 + eq5
  }

  result$time = qx
  result$bandwidth_xM_seq = bandwidth_xM_seq
  result$bandwidth_xm_seq = bandwidth_xm_seq
  result$bandwidth_M = bandwidth_M_temp
  result$bandwidth_m = bandwidth_m_temp

  res_list$result = result

  res_list$yM = yM
  res_list$ym = ym
  res_list$xM = xM
  res_list$xm = xm
  res_list$time_M = time_M
  res_list$time_m = time_m
  res_list$qx = qx
  res_list$varying = varying
  res_list$bandwidth_m = bandwidth_m
  res_list$bandwidth_M = bandwidth_M
  res_list$bandwidth_xm = bandwidth_xm
  res_list$bandwidth_xM = bandwidth_xM
  res_list$hatcoeff1m = hatcoeff1m
  res_list$pred_X1M = pred_X1M
  res_list$pred_X1m = pred_X1m
  res_list$names = colnames(xM_model)
  res_list
}
