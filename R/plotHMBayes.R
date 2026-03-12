#last updated 11/29/2025 by Quincy Mo at Moffitt Cancer Center

plotHMBayes = function(fit, xList, type = c("gaussian", "binomial", "poisson"),
    sample.order = NULL, feature.order = NULL, dist.method="euclidean", hclust.method="ward.D",
    sparse = NULL, threshold=rep(0.5,length(xList)), feature.scale = rep(F,length(xList)), 
    col.scheme=rep(list(bluered(256)),length(xList)),width=5,chr=NULL,plot.chr=NULL,cap=rep(0,length(xList))) 
{
    m = length(xList)
    if(m > length(type)){
        stop("Error:  data type is missing for some data. \n")        
    }
    
    dttype = c("gaussian","binomial","poisson")
    if(!all(type %in% dttype)){
        cat("Error: ",type[!all(type %in% dttype)],"\n")
        stop("Allowed data types are gaussian, binomial and poisson. \n")
    }

    if (is.null(feature.order)) {
        feature.order = rep(T, m)
    }
    if (is.null(feature.scale)) {
      feature.scale = rep("none", m)
    }
    if (is.null(sparse)) {
        sparse = rep(F, m)
    }
    if (is.null(plot.chr)) {
        plot.chr = rep(F, m)
    }
    
    clusters = fit$clusters
    k = length(unique(clusters))
    if (is.null(sample.order)) {
        sorder = order(clusters)
    }else {
        sorder = sample.order
    }

    n = dim(xList[[1]])[1]
    a = clusters[sorder]
    l = length(a)
    brkpoints = which(a[2:l] != a[1:(l - 1)])
    cluster.start = c(1, brkpoints + 1)
    my.panel.levelplot <- function(...) {
        panel.levelplot(...)
        panel.abline(v = (cluster.start[-1] - 0.5), col = "black", 
            lwd = 1, lty = 1)
        panel.scales = list(draw = FALSE)
    }
    
    image.data = alist()
    feature.hclust = alist()
    for (i in 1:m) {
        pp = fit$beta.pp[[i]]
        upper = threshold[i]
        #cat(i," ", sum(pp > upper),"\n")
        if(sparse[i] == T & sum(pp > upper) > 1){
            image.data[[i]] = xList[[i]][sorder, which(pp > upper)]
        }else{
          warning(paste("No features selected in dataset", i),call. = FALSE)
          image.data[[i]] = xList[[i]][sorder, ]
        }
        
        if(feature.scale[i] == T){
          image.data[[i]]=scale(image.data[[i]],center = TRUE, scale = TRUE)
        }
        
        if(feature.order[i]==T){
          if(dist.method == "correlation"){
            diss=as.dist(1-cor(image.data[[i]],use="na.or.complete"))       
          }else{
            diss = dist(t(image.data[[i]]))
          }
          feature.hclust[[i]] = hclust(diss,method=hclust.method)
          image.data[[i]]=(image.data[[i]])[,(feature.hclust[[i]])$order]
        }
        
        if (plot.chr[i] == T) {
            if (sparse[i]) {
                chr = chr[which(pp > upper)]
            }
            len = length(chr)
            chrom.ends <- rep(NA, length(table(chr)))
            d = 1
            for (r in unique(chr)) {
                chrom.ends[d] <- max(which(chr == r))
                d = d + 1
            }
            chrom.starts <- c(1, chrom.ends[-length(table(chr))] + 
                1)
            chrom.mids <- (chrom.starts + chrom.ends)/2
            my.panel.levelplot.2 <- function(...) {
                panel.levelplot(...)
                panel.abline(v = (cluster.start[-1] - 0.5), col = "black", 
                  lwd = 1, lty = 1)
                panel.abline(h = len - chrom.starts[-1], col = "gray", 
                  lwd = 1)
                panel.scales = list(x = list(), y = list(at = len - 
                  chrom.mids), z = list())
            }
            my.panel = my.panel.levelplot.2
            scales = list(x = list(draw = F), y = list(at = len - 
                chrom.mids, labels = names(table(chr))), z = list(draw = F))
        }else{
           my.panel = my.panel.levelplot
           scales = list(draw = F)
        }

        image.data[[i]] = as.matrix(rev(as.data.frame(image.data[[i]])))
        #reverse the rows so the image will be top-down ordering
        
        if (type[i] == "binomial") {
            colorkey = list(space = "right", height = 0.3, at = c(0, 
                0.5, 1), tick.number = 1)
        }else{
            colorkey = list(space = "right", height = 0.3, tick.number = 5)
        }

        if(cap[i] > 0){
            lattice.options(
                layout.heights=list(bottom.padding=list(x=-0.85), top.padding=list(x=-0.85)),
                layout.widths=list(left.padding=list(x=0), right.padding=list(x=0)))
            cut = quantile(image.data[[i]], prob = cap[i], na.rm = T)
            #cat(cut, " ")
            p = levelplot(image.data[[i]], panel = my.panel, scales = scales, 
                col.regions = col.scheme[[i]], at = c(-Inf, seq(-cut, 
                  cut, length = 256), Inf), xlab = "", ylab = "", 
                colorkey = colorkey)
        }else{
            lattice.options(
                layout.heights=list(bottom.padding=list(x=-0.85), top.padding=list(x=-0.85)),
                layout.widths=list(left.padding=list(x=0), right.padding=list(x=0)))
            p = levelplot(image.data[[i]], panel = my.panel, scales = scales, 
                col.regions = col.scheme[[i]], xlab = "", ylab = "", 
                colorkey = colorkey)
        }
        if (i == m) {
            print(p, split = c(1, i, 1, m), more = F, panel.width = list(width, 
                "inches"))
        }else{
            print(p, split = c(1, i, 1, m), more = T, panel.width = list(width, 
                "inches"))
        }
    }
    list(feature.hclust, image.data)
}
