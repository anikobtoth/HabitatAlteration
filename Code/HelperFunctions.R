library(sp)
library(BEST)
library(BayesianFirstAid)
library(rlang)

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
    
    keep1 <- coords1[which(apply(dists, 2, min)<=mindist), 1]
    keep2 <- coords2[which(apply(dists, 1, min)<=mindist), 1]
    
    return(list(keep1, keep2)) }
  }

# Prepare for plotting - calculate observed to expected ratios
obsDexp <- function(d, split.var, data.var, ...){
  group.vars <- enquos(...)
  xvar<- quo(!! sym(paste0(data.var, ".x"))) 
  yvar<- quo(!! sym(paste0(data.var, ".y"))) 
  
  quos_text <- function(qs) {
    unlist(lapply(seq_along(qs), function(i) quo_text(qs[[i]])))}
  
  d <- d %>% filter(!is.na(diet.match)) %>% split(.[split.var])
  m <- full_join(d[[1]], d[[2]], by = quos_text(group.vars))
 
  obsexp <- m %>% select(!!!group.vars, subsample = subsample.x, expected = !!xvar, observed = !!yvar) %>% 
    mutate(obsDexp = observed/expected) %>% 
    separate(Taxon_status, c("taxon", "status"), sep = "_") %>%
    select(-observed, -expected) %>% 
    spread(status, value = obsDexp)
  return(obsexp)
}


obsDexp2 <- function(d, split.var, data.var, ...){
  group.vars <- enquos(...)
  xvar<- quo(!! sym(paste0(data.var, ".x"))) 
  yvar<- quo(!! sym(paste0(data.var, ".y"))) 
  
  quos_text <- function(qs) {
    unlist(lapply(seq_along(qs), function(i) quo_text(qs[[i]])))}
  
  d <- d %>% filter(!is.na(diet.match)) %>% split(.[split.var])
  m <- full_join(d[[1]], d[[2]], by = quos_text(group.vars))
  
  obsexp <- m %>% select(!!!group.vars, subsample = subsample.x, expected = !!xvar, observed = !!yvar) %>% 
    mutate(obsDexp = observed/expected) %>% 
    separate(Taxon_status, c("taxon", "status"), sep = "_") %>%
    select(-observed, -expected) %>% 
    spread(diet.match, value = obsDexp) %>% 
    mutate(diff = Same - Different) %>% 
    select(-Same, -Different) %>%
    spread(status, value = diff)
  return(obsexp)
}


cont_table <- function(x){ #simpairs function, simpairs only out
  samples = ncol(x)  #S
  a = matrix(nrow=nrow(x),ncol=nrow(x),data=0)
  occs = array()
  
  #Calculate overlap
  for (i in 2:nrow(x))  {
    for (j in 1:(i-1))
    {
      a[i,j] = length(which(x[i,] > 0 & x[j,] > 0)) # B
    }
  }

  a <- as.dist(a, diag = F, upper = F)
  
  l <- dist2edgelist(a, x)
  s <- rowSums(x) %>% data.frame()
  
  t <- merge(l, s, by.x = "Sp1", by.y = 0)
  t <- merge(t, s, by.x = "Sp2", by.y = 0)
  t <- t %>% select(Sp1, Sp2, ..x, ..y, Z.Score)
  t$samples <- ncol(x)
  names(t) <- c("Sp1", "Sp2", "presSp1", "presSp2", "presBoth", "samples")
  t$absSp1 <- t$samples - t$presSp1
  t$absSp2 <- t$samples - t$presSp2
  
  t$presSp1absSp2 <- t$presSp1- t$presBoth
  t$presSp2absSp1 <- t$presSp2- t$presBoth
  
  t$absBoth <- t$samples - t$presBoth - t$presSp1absSp2 - t$presSp2absSp1
  
  return(t)
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

## Single run of FETmP
FETmP <- function(Talt, Tunalt, altered, unaltered){
  occurrences <- altered + unaltered
  p <- choose(Talt, 0:altered) * choose(Tunalt, occurrences:unaltered) / choose(Talt+Tunalt, occurrences)
  return(sum(p)-0.5*last(p))
  
}

## FETmP along vectors
FETmP_ <- function(Talt, Tunalt, altered, unaltered){
  out <- numeric()
  for(i in seq_along(Talt)){
    out[i] <- FETmP(Talt[i], Tunalt[i], altered[i], unaltered[i])
  }
  return(out)
}

#### Forbes similarity ###
# formats beta diversity results and organises them by altered and unaltered site pairings. 
beta.types <- function(PAn, unalt_sites){
  beta <- map(PAn, ochiaiMatrix) %>% map(as.dist, upper = F) %>% map2(.y = map(PAn, ~t(.)), dist2edgelist)
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
    s <- ceiling(cJ1(n))
    q <- ceiling(s / scale)
    m <- n
    m[m > 2] <- 3
    rarefy(m,q) / (1 - exp(lchoose(3 * s - 3,q) - lchoose(3 * s,q)))
  }
  
## From John Alroy (2020)####

### cj1 to replace squares

cJ1rich<-function(n)	{
  n <- n[n>0]
  S <- length(n)
  s1 <- sum(n == 1)
  if (s1 == S)
    return(NA)
  if (s1 == 0)
    return(S)
  l <- log(sum(n) / s1)
  return((S + s1) / (1 - exp(-l) + l * exp(-l)))
}

## Ochiai to replace the forbes index ####
ochiai<-function(x,y)	{
  if (is.numeric(x) && is.numeric(y) && min(x) == 0 && min(y) == 0 && length(x) == length(y))	{
    a <- length(which((x * y) > 0))
    b <- length(which(x > 0)) - a
    c <- length(which(y > 0)) - a
  } else	{
    a <- length(na.omit(match(x,y)))
    b <- length(x) - a
    c <- length(y) - a
  }
  return(a / ((a + b) * (a + c))^0.5)
}

ochiaiMatrix<-function(x)	{
  x[is.na(x)] <- 0
  m <- matrix(nrow=ncol(x),ncol=ncol(x))
  for (i in 1:ncol(x))
    for (j in 1:ncol(x))
      m[i,j] <- ochiai(x[,i],x[,j])
  rownames(m) <- colnames(x)
  colnames(m) <- colnames(x)
  return(m)
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
percent.occupancy.by.guild <- function(t, gld, taxon, sitedat){
  pr <- data.frame(row.names = colnames(t))
  for(i in 1:length(gld))  {
    temp <- t[rownames(t) %in% rownames(spp[which(spp$guild == gld[i] & spp$life.form == taxon),]), ] 
    temp <- colSums(temp)/nrow(temp)
    pr[,i] <- temp  
  }
  
  colnames(pr) <- gld
  pr <- merge(pr, sitedat[,c("sample.no", "altered.habitat")], by.x = 0, by.y = "sample.no") %>% namerows

  prm <- melt(pr, id.vars = "altered.habitat")
  prm$altered_habitat <- factor(prm$altered.habitat, levels = c("combined", "cropland", "disturbed forest", 
                                                                "fragment", "inhabited area", "pasture", 
                                                                "plantation", "secondary forest", 
                                                                "unaltered")) %>%
    replace_na("unaltered")
  
    return(prm)
}

# Randomization significance test ######

besttest <- function(obsexp, split.var,  ...){
  group.vars <- enquos(...)
  alt <- obsexp %>% split(.[split.var]) %>% 
    purrr::map(~group_by(., !!! group.vars)) %>%
    purrr::map(~group_map(.,~pull(., altered))) 
    
  unalt <- obsexp %>% split(.[split.var]) %>% 
    purrr::map(~group_by(., !!! group.vars)) %>%
    purrr::map(~group_map(.,~pull(., unaltered))) 
  b <- list(unalt = map2(unalt$Different, unalt$Same, 
                         function(x, y) if(all(x) != 0 && all(y) != 0 && length(x[!is.na(x)]) > 1 && length(y[!is.na(y)]) > 1) {
                           return(BESTmcmc(x, y))
                           } else {return(NULL)}), 
            alt   = map2(alt$Different, alt$Same, 
                         function(x, y) if(all(x) != 0 && all(y) != 0 && length(x[!is.na(x)]) > 1 && length(y[!is.na(y)]) > 1) {
                           return(BESTmcmc(x, y))
                           } else {return(NULL)}))
  return(b)
}

bayesPairedTtest <- function(obsexp, split.var,  ...){
  group.vars <- enquos(...)
  
  data <- melt(obsexp) %>% spread(!!split.var, value)
  n <- data %>% group_by(., variable, !!! group.vars) %>% summarise(placeholder = "") %>% unite(name, sep = "_") %>% pull(name)
  b <- data %>% group_by(., variable, !!! group.vars) %>% 
  group_map(~if(length(.$Different[!is.na(.$Different)]) > 1 && 
                 length(.$Same[!is.na(.$Same)]) > 1 &&
                all(.$Different) != 0 && all(.$Same) != 0) {
    bayes.t.test(.$Different, .$Same, paired = TRUE)}
    else{return(NULL)}, keep = TRUE) %>% setNames(n) 
  b <- b[!sapply(b, is_null)]
  return(b)
}



