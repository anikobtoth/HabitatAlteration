library("sp")

##### DATA MANIPULATION #####
### Change character vectors in df to factors
tofac <- function(df){
  df[,map_lgl(df, is.character)] <- map(df[,map_lgl(df, is.character)], factor)
  return(df)
}
  
### Abundance to presence-absence 
# matrix list
tobinary <- function(PA.LIST){
  binary <- lapply(PA.LIST, function(x) {
    x <- x/x
    x[is.na(x)] <- 0
    return(x)
  })
  return(binary)
}

# single matrix
tobinary.single <- function(x) {
  x <- x/x
  x[is.na(x)] <- 0
  return(x)}

# First column to rownames
namerows <- function(table){
  rownames(table) <- table[,1]
  table <- table[,2:ncol(table)]
  return(table)
}

# remove empty rows and columns in 1 matrix
# or remove rows and columns with too few observations
clean.empty <- function(x, mincol = 1, minrow = 1){
  x <- x[which(rowSums(x) > minrow-1),]
  x <- x[,which(colSums(x) > mincol-1)]
  return(x)
}

# triangular distance matrix to long format
dist2edgelist <- function(z, sppDat){  #edge list with link types attached
  k3 = as.matrix(z)
  dimnames(k3) <- list(rownames(sppDat), rownames(sppDat)) 
  
  xy <- t(combn(colnames(k3), 2))
  k3 <- data.frame(xy, dist=k3[xy], stringsAsFactors = F)
  
  k3 <- data.frame(k3, id = paste(k3$X1, k3$X2, sep = "-"), stringsAsFactors = F)
  colnames(k3) <- c("Sp1", "Sp2", "Z.Score", "id")
  
    return(k3)
}

# Biogeographic matching algorithm (see supplement)
matchbiogeo <- function(coords1, coords2) {
  if(min(nrow(coords1), nrow(coords2))==0) {return(0)
  }else{
    dists <- apply(coords1, 1, function(x) spDistsN1(pts = as.matrix(dplyr::select(coords2, longitude, latitude)), pt = as.numeric(x[3:2]), longlat = TRUE))
    mindist <- min(max(apply(dists,1,min)), max(apply(dists, 2, min))) # take the maximum distance of the closest point of other type for each type of point, take the lesser of these.
    
    keep1 <- coords1$siteid[which(apply(dists, 2, min)<=mindist)]
    keep2 <- coords2$siteid[which(apply(dists, 1, min)<=mindist)]
    
    return(list(keep1, keep2)) }
  }

##### ANALYSES ######
# FETmP
simpairs <- function(x){ #simpairs function, simpairs only out
  samples = ncol(x)  #S
   z = matrix(nrow=nrow(x),ncol=nrow(x),data=0)
  occs = array()
 
  #convert to P/A. Occs = rowsums of PA matrix.
  x <- x/x
  x[is.na(x)] <- 0
  occs <- rowSums(x)
  
  #SimPairs Algorithm
  for (i in 2:nrow(x))  {
    for (j in 1:(i-1))
    {
      a = length(which(x[i,] > 0 & x[j,] > 0)) # B
      
      #simpairs
      for (k in 0:a)
        z[i,j] = z[i,j] + choose(occs[j] , k) * choose(samples - occs[j] , occs[i] - k) / choose(samples , occs[i])
      z[i,j] = z[i,j] - choose(occs[j] , a) * choose(samples - occs[j] , occs[i] - a) / choose(samples , occs[i]) / 2
      z[i,j] = qnorm(z[i,j])
      z[j,i] = z[i,j]
    }
  }
  print("check")
  return(as.dist(z, diag = F, upper = F))
}

resamp <- function(PA, sites = 50, reps = 100){
  samp <- list()
  for(j in 1:reps){
    samp[[j]] <- PA[, sample(1:ncol(PA), sites, replace = F)]
    samp[[j]] <- samp[[j]][which(rowSums(samp[[j]]) > 0),]
  }
  return(samp)
}

#### Forbes similarity ###
# formats beta diversity results and organises them by altered and unaltered site pairings. 
beta.types <- function(PAn, unalt_sites){
  beta <- map(PAn, forbesMatrix) %>% map(as.dist, upper = F) %>% map2(.y = map(PAn, ~t(.)), dist2edgelist)
  beta <- bind_rows(beta, .id = 'taxon')
  beta$unalt1 <- beta$Sp1 %in% unalt_sites
  beta$unalt2 <- beta$Sp2 %in% unalt_sites
  beta$unalt.pair <- paste(beta$unalt1, beta$unalt2, sep = "-")
  beta$unalt.pair[beta$unalt.pair == "TRUE-TRUE"] <- "Unaltered-Unaltered"
  beta$unalt.pair[beta$unalt.pair == "FALSE-FALSE"] <- "Altered-Altered"
  beta$unalt.pair[beta$unalt.pair == "FALSE-TRUE"] <- "Unaltered-Altered"
  beta$unalt.pair[beta$unalt.pair == "TRUE-FALSE"] <- "Unaltered-Altered"
  #ggplot(beta[beta$taxon == "bat",], aes(x = Z.Score, col = alt.pair)) + geom_density(size = 1)
  return(beta)
}

# Squares richness estimator, Alroy 2018
squares<-function(n) {
  n <- n[n>0] # removes any non-sampled species from calculation
  S <- length(n)
  N <- sum(n)
  s1 <- length(which(n == 1))
  if (s1 == S)
    return(NA)
  return(S + s1^2 * sum(n^2) / (N^2 - S * s1))
}

rscale<-function(n,scale=2) {
  s <- ceiling(squares(n))
  q <- ceiling(s / scale)
  m <- n
  m[m > 2] <- 3
  rarefy(m,q) / (1 - exp(lchoose(3 * s - 3,q) - lchoose(3 * s,q)))
}

### Chao1 richness estimator
chao1 <- function(n) {
  n <- n[n>0]
  S <- length(n)
  s1 <- length(which(n == 1))
  s2 <- length(which(n == 2))
  if (s2 == 0)
    return(NA)
  return(S + (s1^2 / (2*s2)))
}

#  Emulate ggplot default colours
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

### matrix fill
matfill <- function(m){
  m <- m/m
  m[is.na(m)] <- 0
  return(sum(m)/(nrow(m)*ncol(m)))
}

### return percent and mean aggregations and segregations, for easier use in dplyr pipes.####
percpos <- function(x)
  length(which(x>0))/length(x)

percneg <- function(x)
  length(which(x<0))/length(x)

meanpos <- function(x)
  mean(x[which(x>0)])

meanneg <- function(x)
  mean(x[which(x<0)])

# identify pairs as aggregations or segregations ####
posnegzero <- function(x){
  out <- x > 0
  out[which(x>0)] <- "Aggregation"
  out[which(x<0)] <- "Segregation"
  out[which(x==0)] <- "ZERO"
  out
}

# Percent occupancy by guild
percent.occupancy.by.guild <- function(t, gld, taxon){
  pr <- data.frame(row.names = colnames(t))
  for(i in 1:length(gld))  {
    temp <- t[rownames(t) %in% rownames(spp[which(spp$guild == gld[i] & spp$life.form == taxon),]), ] 
    temp <- colSums(temp)/nrow(temp)
    pr[,i] <- temp  
  }
  
  colnames(pr) <- gld
  pr$altered_habitat <- sitedat$altered_habitat[sitedat$siteid %in% rownames(pr)]
  
  prm <- melt(pr, id.vars = "altered_habitat")
  prm$altered_habitat <- factor(prm$altered_habitat, levels = c("cropland", "disturbed forest", "fragment", "pasture", "plantation", "rural", "secondary forest", "suburban", "urban"))
  
  return(prm)
}

