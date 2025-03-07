library(sp)
library(rlang)

#### TAXONOMY ######

# get a species taxonomy from the open tree of life

get_otl_taxonomy <- function(spp, lifeform = c("bat", "bird"), contxt = c("Mammals", "Birds")){
  library(rotl)
  library(ape)
  
  spp$verbatimScientificName <- rownames(spp)
  
  tax_otl <- map2(lifeform, contxt, ~spp %>% filter(life.form == .x) %>% 
                    pull(species) %>% tnrs_match_names(context = .y) %>%
                    mutate(search_string = str_to_sentence(search_string)))
  
  return(tax_otl)
  
}

# get pairwise phylogenentic distances from open tree of life for matched taxon names

get_pairwise_dist <- function(tax, contxt){
  intree <- replace_na(tax$ott_id, -9999) %>% is_in_tree()
  
  tree <- tax[intree,] %>% 
    filter(!flags %in% c("hidden", "major_rank_conflict")) %>% pull(unique_name) %>% 
    tnrs_match_names(context = contxt) %>% # need to re-match names after filtering since ids are apparently pulled from the larger object, not the table.
    ott_id() %>% tol_induced_subtree(label = "name")
  
  dist <- tree %>% compute.brlen() %>% cophenetic() %>% 
    data.frame(Sp2 = rownames(.)) %>% as_tibble() %>% 
    pivot_longer(cols = contains("_"), names_to = "Sp1", values_to = "dist") %>% 
    distinct() %>% select(Sp1, Sp2, dist) %>% filter(dist > 0) %>% 
    mutate(Sp1 = str_replace(Sp1, "_", " "), 
           Sp2 = str_replace(Sp2, "_", " "))
  
  return(dist)
}

##### DATA MANIPULATION #####
### Change character vectors in df to factors
tofac <- function(df){
  df[,map_lgl(df, is.character)] <- map(df[,map_lgl(df, is.character)], factor)
  return(df)
}

# return non-overlapping values of two vectors USED
outersect <- function(x, y) {
  sort(c(setdiff(x, y),
         setdiff(y, x)))
}
### Abundance to presence-absence USED
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

# First column to rownames USED
namerows <- function(table){
  rownames(table) <- pull(table, 1)
  table <- table[,2:ncol(table)]
  return(table)
}

# remove empty rows and columns in 1 matrix
# or remove rows and columns with too few observations USED
clean.empty <- function(x, mincol = 1, minrow = 1){
  x <- x[which(rowSums(x) > minrow-1),]
  x <- x[,which(colSums(x) > mincol-1)]
  return(x)
}

# triangular distance matrix to long format USED
dist2edgelist <- function(z, sppDat){  #edge list with link types attached
  k3 = as.matrix(z)
  dimnames(k3) <- list(rownames(sppDat), rownames(sppDat)) 
  
  xy <- t(combn(colnames(k3), 2))
  k3 <- data.frame(xy, dist=k3[xy], stringsAsFactors = F)
  
  k3 <- data.frame(k3, id = paste(k3$X1, k3$X2, sep = "-"), stringsAsFactors = F)
  colnames(k3) <- c("Sp1", "Sp2", "Z.Score", "id")
  
  return(k3)
}

# Biogeographic matching algorithm (see supplement) USED
matchbiogeo <- function(coords1, coords2) {
  if(min(nrow(coords1), nrow(coords2))==0) {return(0)
  }else{
    dists <- apply(coords1, 1, function(x) spDistsN1(pts = as.matrix(dplyr::select(coords2, longitude, latitude)), pt = as.numeric(x[c("longitude", "latitude")]), longlat = TRUE))
    mindist <- min(max(apply(dists,1,min)), max(apply(dists, 2, min))) # take the maximum distance of the closest point of other type for each type of point, take the lesser of these.
    
    keep1 <- coords1[which(apply(dists, 2, min)<=mindist), "p.sample"]
    keep2 <- coords2[which(apply(dists, 1, min)<=mindist), "p.sample"]
    
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

# USED
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
  
  #t$presSp1absSp2 <- t$presSp1- t$presBoth
  #t$presSp2absSp1 <- t$presSp2- t$presBoth
  
  #t$absBoth <- t$samples - t$presBoth - t$presSp1absSp2 - t$presSp2absSp1
  
  return(t)
}

# USED
diet_cat <- function(x, g, related = TRUE){
  x$diet.Sp1 <- g[x$Sp1]
  x$diet.Sp2 <- g[x$Sp2]
  x <- na.omit(x)
  
  x$diet.pair <- map2_chr(x$diet.Sp1, x$diet.Sp2, function(x, y) c(x,y) %>% sort() %>% paste(collapse = "-"))
  
  rltd <- unique(x$diet.pair) %>% strsplit("") %>% map(duplicated) %>% map_int(sum) %>% `==`(1)
  rltd <- unique(x$diet.pair)[which(nchar(unique(x$diet.pair)) > 3 & rltd)]
  
  if(!related) x <- x[!which(paste(x$diet.Sp1, x$diet.Sp2, sep = "-") %in% rltd),]
  
  x$diet.match <- as.numeric(x$diet.Sp1 == x$diet.Sp2)
  x$diet.match[x$diet.match == 0] <- "Different"
  x$diet.match[x$diet.match == 1] <- "Same"
  x$diet.match[x$diet.pair %in% rltd] <- "Related"
  x$diet.match[x$diet.Sp1 == "O" | x$diet.Sp2 == "O"] <- "Related"  #all pairs including an omnivore are "related"
  
  return(x)
  
}

# fixed-fixed randomisation function (uses vegan package) USED 
rand_mat <- function(x, method = "curveball", i = 1){
  y <- nullmodel(x, method = method) |> simulate(nsim = 1, burnin = 1e8, seed = i)
  y <- y[,,1]
  attributes(y) <- attributes(x)
  return(y)
}

##### ANALYSES ######

## Single run of FETmP USED
FETmP <- function(Talt, Tunalt, altered, unaltered){
  occurrences <- altered + unaltered
  p <- choose(Talt, 0:altered) * choose(Tunalt, occurrences:unaltered) / choose(Talt+Tunalt, occurrences)
  return(sum(p)-0.5*last(p))
}
## FETmP along vectors USED
FETmP_ <- function(Talt, Tunalt, altered, unaltered){
  out <- numeric()
  for(i in seq_along(Talt)){
    out[i] <- FETmP(Talt[i], Tunalt[i], altered[i], unaltered[i])
  }
  return(out)
}

# Fishers noncentral hypergeometric probability mass function
log_sum_exp <- function(x) log(sum(exp(x - max(x)))) + max(x)
fnchypergeo_lpmf_ <- function(sb, s1, s2, n_site, theta) {
  sb_min <- max(c(0,s1 + s2 - n_site))
  sb_max <- min(c(s1,s2))
  u <- ll <- ans <- numeric()
  for(i in sb_min:sb_max) {
    u[i - sb_min + 1] = i;
    ll[i - sb_min + 1] = lchoose(s1, i) + lchoose(n_site - s1, s2 - i)
  }
  if(length(sb) > 1) {
    for(i in 1:length(sb)) {
      ans[i] = ll[sb[i] - sb_min + 1] + sb[i]*theta - log_sum_exp(ll + u*theta)
    }
  } else {
    for(i in 1:length(theta)) {
      ans[i] = ll[sb - sb_min + 1] + sb*theta[i] - log_sum_exp(ll + u*theta[i])
    }
  }    
  return(ans)
}
fnchypergeo_lpmf <- function(sb, s1, s2, n_site, theta){
  d <- data.frame(sb, s1, s2, n_site, theta)
  Nab_theta <- numeric()
  for(s in 1:nrow(d)){
    Nab_theta[s] <- fnchypergeo_lpmf_(d$sb[s], d$s1[s], d$s2[s], d$n_site[s], d$theta[s])
  }
  return(Nab_theta)
}

#### similarity ####

## CJ1 From John Alroy (2020) USED

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

### Chao1 richness estimator USED
chao1 <- function(n) {
  n <- n[n>0]
  S <- length(n)
  s1 <- length(which(n == 1))
  s2 <- length(which(n == 2))
  if (s2 == 0)
    return(NA)
  return(S + (s1^2 / (2*s2)))
}

# Plotting 
library(rstan)
library(ggridges)

format_stanfit <- function(stanfit, name = "mu"){
  stanfit %>% rstan::extract() %>% `[[`(name) %>% data.frame() %>% 
    pivot_longer(names_to = "group", cols = 1:4, values_to = name) %>% 
    mutate(group = as.factor(group) %>% recode(`X1` = "Intact control", 
                                               `X2` = "Altered control", 
                                               `X3` = "Intact intraguild", 
                                               `X4` = "Altered intraguild")) %>%
    separate(group, into = c("status", "interaction"), remove = FALSE)
}  
plot_stanfit <- function(formattedstanfit){
  formattedstanfit %>%
    ggplot(aes(x = mu, y = status, fill = interaction, col = interaction)) + 
    geom_density_ridges(lwd = 1, alpha = 0.4) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_color_manual(values = c("#045FB4", "#2E2E2E")) +
    scale_fill_manual(values = c("#045FB4", "#2E2E2E")) + 
    coord_cartesian(clip = "off") +
    labs(x = "group θ", y = element_blank()) + theme_ridges()
}

# Easy panel labels in base R plot
# from https://logfc.wordpress.com/2017/03/15/adding-figure-labels-a-b-c-in-the-top-left-corner-of-the-plotting-region/
fig_label <- function(text, region="figure", pos="topleft", cex=NULL, ...) {
  
  region <- match.arg(region, c("figure", "plot", "device"))
  pos <- match.arg(pos, c("topleft", "top", "topright", 
                          "left", "center", "right", 
                          "bottomleft", "bottom", "bottomright"))
  
  if(region %in% c("figure", "device")) {
    ds <- dev.size("in")
    # xy coordinates of device corners in user coordinates
    x <- grconvertX(c(0, ds[1]), from="in", to="user")
    y <- grconvertY(c(0, ds[2]), from="in", to="user")
    
    # fragment of the device we use to plot
    if(region == "figure") {
      # account for the fragment of the device that 
      # the figure is using
      fig <- par("fig")
      dx <- (x[2] - x[1])
      dy <- (y[2] - y[1])
      x <- x[1] + dx * fig[1:2]
      y <- y[1] + dy * fig[3:4]
    } 
  }
  
  # much simpler if in plotting region
  if(region == "plot") {
    u <- par("usr")
    x <- u[1:2]
    y <- u[3:4]
  }
  
  sw <- strwidth(text, cex=cex) * 60/100
  sh <- strheight(text, cex=cex) * 60/100
  
  x1 <- switch(pos,
               topleft     =x[1] + sw, 
               left        =x[1] + sw,
               bottomleft  =x[1] + sw,
               top         =(x[1] + x[2])/2,
               center      =(x[1] + x[2])/2,
               bottom      =(x[1] + x[2])/2,
               topright    =x[2] - sw,
               right       =x[2] - sw,
               bottomright =x[2] - sw)
  
  y1 <- switch(pos,
               topleft     =y[2] - sh,
               top         =y[2] - sh,
               topright    =y[2] - sh,
               left        =(y[1] + y[2])/2,
               center      =(y[1] + y[2])/2,
               right       =(y[1] + y[2])/2,
               bottomleft  =y[1] + sh,
               bottom      =y[1] + sh,
               bottomright =y[1] + sh)
  
  old.par <- par(xpd=NA)
  on.exit(par(old.par))
  
  text(x1, y1, text, cex=cex, ...)
  return(invisible(c(x,y)))
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

