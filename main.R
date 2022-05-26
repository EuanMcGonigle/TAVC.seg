require(mosum)
require(breakfast)

robust.tavc.est = function(x, L, M = floor(2.5*sqrt(length(x))), b.max = min(floor(L/2),floor(M/2))-1, global.est = TRUE,
                            window.len = 10*(min(floor(L/2),floor(M/2))), v.est.method = c("trimmed.mean","median","MAD","median.pairwise")[1]){

  #main function for robustly estimating the TAVC at scale L of a univariate time series.
  # Parameters needed for the function are:
  #
  # x - the data
  # L - the scale at which to estimate the TAVC. L=2G, the block size used in estimation.
  # M - the maximum estimation scale. If L is set such that L>M, then the TAVC is estimated at scale M.
  # b.max - the number of starting points to use for the global TAVC estimator. Helps with robustness.
  # global.est - logical variable, set TRUE to calculate the TAVC estimator assuming stationarity of the error process. Set FALSE to
  #             compute a time-varying estimator of the TAVC.
  # window.len - if global.est = FALSE, then window.len sets the size of window over which the local TAVC is estimated.
  # v.est.method - method of plug-in estimation for the variance parameter v. Recommended to choose trimmed.mean or median.

  # For applying the bottom-up MOSUM or WBS2 algorithms with the robust TAVC estimation procedure, use the functions mosum.tavc and WBS2.tavc
  # NOTE: to apply the robust TAVC estimation with bottom-up MOSUM or WBS2, the code requires the R packages mosum and breakfast


  if (sum(is.na(x))!=0){
    stop("Data contains mising values.")
  }
  if (!is.atomic(x)){
    stop("Data is not atomic")
  }

  data.len = length(x)

  G = floor(L/2)
  N1 = floor((data.len-G)/G)
  Smax = floor(M/2)

  if (b.max>=G){
    warning("b.max is set too large. Setting b.max = G-1 instead.")
    b.max = G-1
  }


  if(global.est==TRUE){
    if(G>=Smax){
      G = Smax
      x.v = mosum(x, G = Smax, var.custom = rep(1,data.len), var.est.method = "custom", boundary.extension = FALSE)$stat^2
    } else{
      x.v = mosum(x, G = G, var.custom = rep(1,data.len), var.est.method = "custom",boundary.extension = FALSE)$stat^2
    }

    var.ests = rep(0,b.max)

    for (l in 1:b.max){
      l.seq = seq(from = G+l-1, by = G, length = N1)
      l.seq = l.seq[l.seq<=data.len-G]

      temp = x.v[l.seq]

      if(v.est.method =="trimmed.mean"){
        v = mean(temp, trim=0.25)
      } else if (v.est.method=="median"){
        v = 2.125365*median(temp)
      } else if (v.est.method=="MAD"){
        v = 3.373791*mad(temp, cons=1)
      } else if (v.est.method=="median.pairwise"){
        d = as.matrix(dist(temp, diag=TRUE,upper=TRUE))
        v = median(apply(d,1,median))*3.211
      } else {
        stop("Choice of v.est.method invalid.")
      }

      a = sqrt(G/data.len)/v

      results =  try(uniroot(robustfun2, x = temp, alpha = a, interval = c(-10,max(temp))),silent = TRUE)
      if(class(results)=="try-error"){
        var.ests[l]=NA
      } else{
        var.ests[l] = results$root
      }

    }

    final.var.est = rep(median(var.ests, na.rm = TRUE), data.len)

  } else{

    if(G>=Smax){
      G=Smax
      x.v = mosum(x, G=Smax, var.custom = rep(1,data.len), var.est.method = "custom",
                  boundary.extension = TRUE)$stat^2

    } else{
      x.v = mosum(x, G=G, var.custom = rep(1,data.len), var.est.method = "custom",
                  boundary.extension = TRUE)$stat^2
    }

    final.var.est = rep(0,data.len)

    N2 = floor((window.len-G)/G)

    time = 1
    l.seq = seq(from = time, by = G, length = N2)
    l.seq = l.seq[l.seq<=data.len]

    while(max(l.seq)<data.len){

      l.seq = l.seq + 1
      temp = x.v[l.seq]

      if(v.est.method =="trimmed.mean"){
        v = mean(temp, trim=0.25)
      } else if (v.est.method=="median"){
        v = 2.125365*median(temp)
      } else if (v.est.method=="MAD"){
        v = v = 3.373791*mad(temp, cons=1)
      } else if (v.est.method=="median.pairwise"){
        d = as.matrix(dist(temp, diag=TRUE,upper=TRUE))
        v = median(apply(d,1,median))*3.211
      } else {
        stop("Choice of v.est.method invalid.")
      }
      a = sqrt(G/window.len)/v

      results =  try(uniroot(robustfun2, x = temp, alpha = a, interval = c(-10,max(temp))),silent = TRUE)
      if(class(results)=="try-error"){
        final.var.est[median(l.seq)]=NA
      } else{
        final.var.est[median(l.seq)] = results$root
      }

      time=time+1

    }

    final.var.est2 = final.var.est[final.var.est!=0]

    #replace NAs if any are present:

    v.nas = which(is.na(final.var.est2))

    if(length(v.nas>0)){
      for(i in 1:length(v.nas)){
        final.var.est2[v.nas[i]] = final.var.est2[which.min(abs(v.nas[i]-which(!is.na(final.var.est2))))]
      }
    }

    suppressWarnings(med.var.est <- runmed(final.var.est2, k = b.max, endrule = "constant"))

    l.m = length(med.var.est)

    final.var.est = c(rep(med.var.est[1],ceiling((data.len-l.m)/2)),med.var.est,rep(med.var.est[length(med.var.est)],floor((data.len-l.m)/2)))

  }

  return(final.var.est)

}


mosum.tavc = function (x, G = bandwidths.default(length(x), G.min = max(20, ceiling(0.05 * length(x)))),
                       threshold = c("critical.value","custom")[1], alpha = 0.1, threshold.function = NULL,
                       criterion = c("eta", "epsilon")[1], eta = 0.4, epsilon = 0.2, global.est = TRUE,
                       v.est.method = c("trimmed.mean","median","MAD","median.pairwise")[1],
                       boundary.extension = FALSE, do.confint = FALSE, level = 0.05, N_reps = 1000,
                       M = floor(sqrt(2.5*length(x))), b.max = min(G,floor(M/2))-1,  N2 = 5, plot = FALSE, ...)
{

  # Multiscale MOSUM procedure, using the TAVC estimator. All parameters inherited from robust.tavc.est and mosum::multiscale.bottomUp, except:
  #
  # N2 - Used instead of window.len, controls the size of window for the local TAVC estimator, needed if global.est = FALSE.
  #     The window size is given by window.len = N_3*L.

  n <- length(x)

  TAVC.output = rep(list(rep(0,n)), length(G))
  test.stats = list()

  Smax = floor(M/2)

  L.values = c(2*G[G<Smax],M)

  num.tavc.ests = length(G[G<Smax])
  if(num.tavc.ests>0){
    for(i in 1:num.tavc.ests){
      TAVC.output[[i]] = robust.tavc.est(x = x, L = 2*G[i], M = M, b.max = G[i]-1, global.est = global.est,
                                         window.len = 2*G[i]*N2, v.est.method = v.est.method)
    }
  }

  TAVC.output[[(num.tavc.ests+1)]] = robust.tavc.est(x = x, L = M, M = M, b.max = Smax-1, global.est = global.est,
                                                     window.len = M*N2, v.est.method = v.est.method)
  k = num.tavc.ests+1
  while(k<length(G)){
    TAVC.output[[(k+1)]] = TAVC.output[[k]]
    k=k+1
  }


  if (class(G) == "integer" || class(G) == "numeric") {
    grid <- mosum:::multiscale.grid(G, method = "concatenate")
  }
  else if (class(G) == "multiscale.grid") {
    if (any(apply(G$grid, 1, diff) != 0)) {
      stop("Expecting a grid of symmetric bandwidths")
    }
    grid <- G
  }
  else stop("Expecting a vector of numbers")
  abs.bandwidth <- all(grid$grid >= 1)
  if (abs.bandwidth) {
    GRID_THRESH <- max(20, 0.05 * n)
  }
  else {
    GRID_THRESH <- 0.05
  }
  if (threshold == "critical.value" & min(grid$grid) < GRID_THRESH) {
    warning("Smallest bandwidth in grid is relatively small (in comparison to n), \n increase the smallest bandwidth or use multiscale.localPrune instead")
  }
  if (threshold != "critical.value" && threshold != "custom") {
    stop("threshold must be either 'critical.value' or 'custom'")
  }
  stopifnot(alpha >= 0 && alpha <= 1)
  stopifnot(eta <= 1 & eta > 0)
  stopifnot(!do.confint || N_reps > 0)
  cpts.complete <- numeric(0)
  bandwidths.complete <- integer(0)
  pValues.complete <- numeric(0)
  jumps.complete <- numeric(0)
  for (i in seq_len(nrow(grid$grid))) {
    G <- grid$grid[[i, 1]]

    if (threshold == "critical.value") {
      m = mosum(x, G = G, ..., threshold = "critical.value", boundary.extension = boundary.extension, alpha = alpha,
                criterion = criterion, eta = eta, epsilon = epsilon, var.est.method = "custom", var.custom = TAVC.output[[i]])

      test.stats[[i]] = m$stat

    }
    else {
      threshold_val <- threshold.function(G, n, alpha)
      m = mosum(x, G = G, ..., threshold = "custom", threshold.custom = threshold_val,boundary.extension = boundary.extension,
                alpha = alpha, criterion = criterion, eta = eta, epsilon = epsilon,var.est.method = "custom", var.custom = TAVC.output[[i]])
    }
    if (!abs.bandwidth)
      G <- floor(G * n)
    cpts <- m$cpts
    cpts.complete <- c(cpts.complete, cpts)
    bandwidths.complete <- c(bandwidths.complete, rep(G,
                                                      length(cpts)))
    pValues.complete <- c(pValues.complete, mosum:::mosum.pValue(m$stat[cpts],
                                                                 n, G))
    jumps.complete <- c(jumps.complete, m$stat[cpts] * sqrt(2/G))
  }
  points <- numeric(0)
  bandwidths <- numeric(0)
  pValues <- numeric(0)
  jumps <- numeric(0)
  cptsInOrder <- seq_len(length(cpts.complete))
  for (i in cptsInOrder) {
    p <- cpts.complete[[i]]
    G <- bandwidths.complete[[i]]
    pVal <- pValues.complete[[i]]
    jmp <- jumps.complete[[i]]
    if (suppressWarnings(min(abs(p - points))) >= eta * G) {
      points <- c(points, p)
      bandwidths <- c(bandwidths, G)
      pValues <- c(pValues, pVal)
      jumps <- c(jumps, jmp)
    }
  }
  cpts.merged <- data.frame(cpts = points, G.left = bandwidths,
                            G.right = bandwidths, p.value = pValues, jump = jumps)
  cpts <- cpts.merged[order(cpts.merged$cpts), ]
  G <- as.vector(grid$grid[, 1])
  if (!abs.bandwidth)
    G <- floor(n * G)
  ret <- structure(list(x = x, cpts = as.numeric(cpts[, 1]),
                        cpts.info = cpts, pooled.cpts = sort(unique(cpts.complete)),
                        G = G, alpha = alpha, threshold = threshold, threshold.function = threshold.function,
                        criterion = "eta", eta = eta, do.confint = FALSE, ci = NA, test.stats = test.stats,
                        TAVC.ests = list(L.vals = L.values, ests = t(matrix(unlist(TAVC.output), ncol = length(TAVC.output))))),
                   class = "multiscale.cpts")
  if (do.confint) {
    ret$ci <- confint.multiscale.cpts(ret, level = level,
                                      N_reps = N_reps)
    ret$do.confint <- TRUE
  }
  ret
}


wbs2.tavc = function(x, R = 100, M = floor(2.5*sqrt(length(x))), v.est.method = c("trimmed.mean","median","MAD","median.pairwise")[1],
                     global.est = TRUE, N2 = 5, min.int.len = floor(0.05*length(x)), th.const = 1.3){

  # WBS2 algorithm using the robust TAVC estimator. All parameters inherited from robust.tavc.est, except:
  #
  # N2 - Used instead of window.len, controls the size of window for the local TAVC estimator, needed if global.est = FALSE.
  #     The window size is given by window.len = N_3*L.
  # R - the number of systematic intervals to use at each iteration of the WBS2 algorithm.
  # min.int.len - the minimum interval length considered for the WBS2 algorithm.
  # th.const - the threshold constant used by WBS2 to test for change points.

  solutions.nested <- TRUE
  solution.set <- list()

  min.len = floor(min.int.len/2)

  Smax = floor(M/2)


  if(min.len>Smax){
    tavc.ests = t(as.matrix(robust.tavc.est(x = x, L = M, M = M, b.max = Smax-1, global.est = global.est,
                                            window.len = M*N2, v.est.method = v.est.method)))
    L.vals = M
  } else{
    tavc.ests = matrix(0,nrow = Smax-min.len+1,ncol = length(x))
    L.vals = 2*(min.len:Smax)
    for(i in 1:(Smax-min.len+1)){
      tavc.ests[i,] = robust.tavc.est(x, L = 2*(min.len+i-1), M = M, b.max = min.len+i-2 , v.est.method = v.est.method,
                                      global.est = global.est, window.len = 2*(min.len+i-1)*N2)
    }
  }


  wbs2.tavc.int = function(x, R = 100, min.len = 0, Smax, tavc.ests){

    n <- length(x)
    if (n <= 2 * min.len + 1) {
      return(matrix(NA, 4, 0))}
    else {
      cpt <- t(systematic.cusums.tavc(x, R, min.len, Smax, tavc.ests)$max.val)
      return(cbind(cpt, wbs2.tavc.int(x = x[1:cpt[2]], R, min.len, Smax, tavc.ests[,1:cpt[2],drop=FALSE]),
                   wbs2.tavc.int(x = x[(cpt[2] + 1):n], R, min.len, Smax,  tavc.ests[,(cpt[2] + 1):n,drop=FALSE]) + c(rep(cpt[2], 3), 0)))
    }

  }

  sol.wbs2 = wbs2.tavc.int(x, R = R, min.len = min.len, Smax=Smax, tavc.ests = tavc.ests)

  ord <- order(abs(sol.wbs2[4,]), decreasing = T)
  sorted.cusums  <- abs(t(sol.wbs2)[ord, , drop = F])
  solution.path <- sorted.cusums[, 2]

  temp = sorted.cusums[,2]
  sorted.cusums[,2] = sorted.cusums[,3]
  sorted.cusums[,3] = temp

  ret = list(solutions.nested = solutions.nested, solution.path = solution.path,
             solution.set = solution.set, x = x, R = R, cands = sorted.cusums,
             method = "wbs2", tavc.ests = list(L.vals = L.vals, ests = tavc.ests))

  class(ret) <- "cptpath"

  cpt.obj = breakfast::model.thresh(ret, sigma = 1, th_const = th.const)

  cpt.obj


}





robustfun = function(x)
{
  if (x >= 1)
  {
    y = log(2)

  }

  if ((x>= 0)&&(x<1))
  {
    y = - log(1-x+x^2/2)

  }
  if ((x<=0)&& (x> -1))
  {
    y = log(1+x+x^2/2)

  }
  if  (x<= -1)
  {
    y = -log(2)

  }
  return(y)

}

robustfun2 = function(x,c,alpha)
{

  y = x
  for (i in (1:length(x)))
  {
    y[i] = alpha^(-1)*robustfun(alpha*(x[i]-c))
  }
  ysum = sum(y)/length(x)
  return(ysum)
}

systematic.cusums.tavc = function(x, R, min.len, Smax, tavc.ests){

  y <- c(0, cumsum(x))
  n <- length(x)
  R <- min(R, (n - 1)*n/2)

  ind <- breakfast:::grid.intervals(n, R)
  R <- dim(ind)[2]
  res <- matrix(0, R, 4)

  interval.lengths = ind[2,]-ind[1,]+1
  if(min.len>Smax){
    tavc.lens = Smax
  } else{
    tavc.lens = min.len:Smax
  }

  ind.tavc.est = rep(0,dim(ind)[2])

  for (i in 1:length(ind.tavc.est)){
    ind.tavc.est[i] = which.min(abs(2*tavc.lens-interval.lengths[i]))
  }

  res[,1:2] <- t(ind)
  res[,3:4] <- t(apply(rbind(ind,ind.tavc.est), 2, max.cusum.tavc, y, min.len, tavc.ests))

  max.ind <- which.max(abs(res[,4]))
  max.val <- res[max.ind, c(1, 3, 2, 4), drop = FALSE]

  list(res = res, max.val = max.val, R.eff = R)

}

max.cusum.tavc = function(ind, y, min.len, tavc.ests) {

  m <- ind[2] - ind[1] + 1

  if(m > 2 * min.len + 1){
    z <- y[(ind[1] + 1):(ind[2] + 1)] - y[ind[1]]
    ip <- sqrt(((m - 1):1) / m / (1:(m - 1))) * z[1:(m - 1)] - sqrt((1:(m - 1))/m/((m - 1):1)) * (z[m] - z[1:(m - 1)])
    ip = ip/sqrt(tavc.ests[ind[3],(ind[1]):(ind[2]-1)])
    ip.max <- which.max(abs(ip[(min.len + 1):(m - min.len - 1)])) + min.len
    return(c(ip.max + ind[1] - 1, abs(ip[ip.max])))
  } else{
    return(c(ind[1], 0))
  }

}











