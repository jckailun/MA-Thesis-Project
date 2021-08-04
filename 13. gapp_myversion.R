gapp = function (x, FUNcluster=kmd, K.max, B = 100, d.power = 1, spaceH0 = c("scaledPCA", 
                                                                  "original"), verbose = interactive(), ...) 
{
  stopifnot(is.function(FUNcluster), length(dim(x)) == 2, 
            K.max >= 2, (n <- nrow(x)) >= 1, ncol(x) >= 1)
  if (B != (B. <- as.integer(B)) || (B <- B.) <= 0) 
    stop("'B' has to be a positive integer")
  cl. <- match.call()
  if (is.data.frame(x)) 
    x <- as.matrix(x)
  ii <- seq_len(n)
  W.k <- function(X, kk) {
    if (kk > 1) {
    clus <- kmd(X, kk) } else {clus = rep.int(1L, nrow(X))}


    0.5 * sum(vapply(split(ii, clus), function(I) {
      xs <- X[I, , drop = FALSE]
      sum(xs^d.power/nrow(xs))
      #sum(dist(xs)^d.power/nrow(xs))
    }, 0))
  }
  logW <- E.logW <- SE.sim <- numeric(K.max)
  if (verbose) 
    cat("Clustering k = 1,2,..., K.max (= ", K.max, "): .. ", 
        sep = "")
  for (k in 1:K.max) logW[k] <- log(W.k(x, k))
  if (verbose) 
    cat("done\n")
  spaceH0 <- match.arg(spaceH0)
  xs <- scale(x, center = TRUE, scale = FALSE)
  m.x <- rep(attr(xs, "scaled:center"), each = n)
  switch(spaceH0, scaledPCA = {
    V.sx <- svd(xs, nu = 0)$v
    xs <- xs %*% V.sx
  }, original = {
  }, stop("invalid 'spaceH0':", spaceH0))
  rng.x1 <- apply(xs, 2L, range)
  logWks <- matrix(0, B, K.max)
  if (verbose) 
    cat("Bootstrapping, b = 1,2,..., B (= ", B, ")  [one \".\" per sample]:\n", 
        sep = "")
  for (b in 1:B) {
    z1 <- apply(rng.x1, 2, function(M, nn) runif(nn, min = M[1], 
                                                 max = M[2]), nn = n)
    z <- switch(spaceH0, scaledPCA = tcrossprod(z1, V.sx), 
                original = z1) + m.x
    for (k in 1:K.max) {
      logWks[b, k] <- log(W.k(z, k))
    }
    if (verbose) 
      cat(".", if (b%%50 == 0) 
        paste(b, "\n"))
  } 
  if (verbose && (B%%50 != 0)) 
    cat("", B, "\n")
  E.logW <- colMeans(logWks)
  SE.sim <- sqrt((1 + 1/B) * apply(logWks, 2, var))
  structure(class = "clusGap", list(Tab = cbind(logW, E.logW, 
                                                gap = E.logW - logW, SE.sim), call = cl., spaceH0 = spaceH0, 
                                    n = n, B = B, FUNcluster = FUNcluster))
}
