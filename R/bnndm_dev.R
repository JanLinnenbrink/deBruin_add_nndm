#' NN distance based block CV
#' @description
#' This function calculates the nearest neighbor distance
#' between training locations and prediction locations (Gij), as well as
#' between different cross-validation folds  (Gjstar). The Wasserstein-Test is used
#' to select the CV-folds yielding the lowest distance between Gjstar and Gij.
#' @param tpoints A sf object containing the training points.
#' @param modeldomain A raster or terra object containing the predictor variables.
#' @param ppoints sf or sfc point object. Contains the target prediction points. Optional. Alternative to modeldomain (see Details).
#' @param samplesize Optional. The number of samplepoints generated in the modeldomain to describe the prediction space.
#' If NULL (default), the samplesize equals the number of cells in the predictor raster divided by ten.
#' @param sampling The type in which the samplepoints are generated. The default value is "Fibonacci", which is recommended for global studies.
#' @param n_steps The number of steps until termination.
#' @param k The  number of clusters to be fitted. Either a vector of values, or a single integer.
#' @param rm Boolean. Should training points be allowed to be excluded? By default TRUE.
#' @param prop_out Number of steps with removal of tpoints.
#' @param max_pout The maximum proportion of a cluster that can be removed.
#' @details The nearest neighbor distance distribution between training and prediction locations,
#' as well as between cross-validation folds are calculated. Based on that, the
#' cross-validation split yielding the lowest distance between those distributions is chosen.
#' @details TBC.
#' @note Experimental.
#' @export
#'
#' @examples
#' \dontrun{
#' library(sf)
#' library(raster)
#' library(ggplot2)
#'
#' ### prepare sample data:
#' dat <- get(load(system.file("extdata","Cookfarm.RData",package="CAST")))
#' dat <- aggregate(dat[,c("DEM","TWI", "NDRE.M", "Easting", "Northing","VW")],
#'    by=list(as.character(dat$SOURCEID)),mean)
#' pts <- dat[,-1]
#' pts <- st_as_sf(pts,coords=c("Easting","Northing"))
#' st_crs(pts) <- 26911
#' studyArea <- raster::stack(system.file("extdata","predictors_2012-03-25.grd",package="CAST"))
#'
#' bnndm_folds <- bnndm(pts, modeldomain=studyArea)
#' bnndm_folds$stat # Integral between the two curves
#' plot(bnndm_folds)
#' plot_geodist(pts, studyArea, cvfolds = bnndm_folds$indx_test)
#' ggplot() +
#'   geom_sf(data = pts, col = bnndm_folds$clusters)
#'
#' #use for cross-validation:
#' library(caret)
#' ctrl <- trainControl(method="cv",
#'    index=bnndm_folds$indx_train)
#' model_bnndm <- train(dat[,c("DEM","TWI", "NDRE.M")],
#'    dat$VW,
#'    method="rf",
#'    trControl = ctrl)
#' model_bnndm
#'}


bnndm <- function(tpoints=NULL,
                  modeldomain=NULL,
                  ppoints=NULL,
                  samplesize=1000,
                  sampling = "regular",
                  n_steps = 50,
                  k = c(2,5),
                  rm = TRUE,
                  prop_out = 0.5,
                  max_pout = 0.5) {




  # input formatting --------------------------------------------------------

  if (is.null(ppoints) & !is.null(modeldomain)) {
    message(paste0(samplesize, "prediction points are sampled from the modeldomain"))
    ppoints <- sampleFromArea(modeldomain, samplesize, type = "geo", variables = NULL, sampling)
  }

  if (any(class(tpoints) %in% "sfc")) {
   tpoints <- sf::st_sf(geom = tpoints)
  }
  if (any(class(ppoints) %in% "sfc")) {
    ppoints <- sf::st_sf(geom = ppoints)
  }
  if(sf::st_crs(tpoints) != sf::st_crs(ppoints)){
   tpoints<- sf::st_transform(tpoints, sf::st_crs(ppoints))
    message("tpoints and ppoints must have the same CRS. tpoints have been transformed.")
  }


  # calculate NN distance of samples to prediction --------------------------

  # Gj: NND function for a cluster per point, i.e. LOO CV
  clust <- 1:nrow(tpoints)
  distmat <- sf::st_distance(tpoints)
  units(distmat) <- NULL
  Gj <- distclust(distmat, clust)

  # Gij: prediction to training NN distances
  Gij <- sf::st_distance(ppoints, tpoints)
  Gij <- apply(Gij, 1, min)
  units(Gij) <- NULL

  # Check if there's is clustering in the data in the first place
  testks <- stats::ks.test(Gj, Gij, alternative = "great")
  if(testks$p.value >= 0.05){

    # if no evidence of clustering found: random 10 fold CV
    k = 10
    blocks <- sample(rep(1:k, ceiling(nrow(tpoints)/k)), size = nrow(tpoints), replace=T)
    Gjstar <- distclust(distmat, blocks)
    ws_dist_min <- twosamples::wass_stat(Gjstar, Gij)[[1]]
    cfolds <- CAST::CreateSpacetimeFolds(data.frame(blocks=blocks), spacevar = "blocks", k = k)

    warning("No evidence of clustering has been found, a random CV assignment is returned")

    # Output
    res <- list(blocks = blocks,
                indx_train = cfolds$index, indx_test = cfolds$indexOut,
                Gij = Gij, Gj = Gj, Gjstar = Gjstar, stat = ws_dist_min)
    class(res) <- c("knndm", "list")
    res

  } else {


    # assign initial groups to the tpoints -----------------------------------

    coords <- sf::st_coordinates(tpoints) |> as.matrix()
    if (length(k) == 1) {
      clust <- as.data.frame(stats::kmeans(x=coords, centers=k)$cluster)
      colnames(clust) <- paste0("groups_", k)
      } else {
        clust <- lapply(k, function(y) {
          c <- stats::kmeans(x=coords, centers = y)
          c$cluster})
        clust <- do.call(cbind.data.frame, clust)
        colnames(clust) <- paste0("groups_", k)
      }

    pts_clust <- cbind(tpoints, clust)

    if (n_steps == 1) {

      # CV based on kmeans-clustering
      pts_id_l <- sf::st_drop_geometry(pts_clust) |>
        as.list()

    } else {

      if (rm==TRUE) {
        # add nearest neighbor distances resulting of kmeans clustering to each point
        pts_clust <- pts_clust[,grepl("groups", colnames(pts_clust))]
        pts_df <- sf::st_drop_geometry(pts_clust)
        pts_df$id <- 1:nrow(pts_clust); pts_clust$id <- 1:nrow(pts_clust)

        pts_clust_l <- lapply(k,  function(x) pts_clust[,c(paste0("groups_", x),"id")])
        pts_df_l <- lapply(k,  function(x) pts_df[,c(paste0("groups_", x),"id")])

        for(i in 1:length(unique(k))) {
          pts_df_l[[i]] <- pts_df_l[[i]][order(pts_df_l[[i]][,1]),]
          pts_df_l[[i]][,"dist"] <- distclust(distmat, pts_df_l[[i]][[1]])
          pts_df_l[[i]] <- pts_df_l[[i]][order(pts_df_l[[i]][["dist"]]),]
          pts_df_l[[i]] <- pts_df_l[[i]][order(pts_df_l[[i]][[paste0("groups_", k[[i]])]]),]
          pts_clust_l[[i]] <- pts_clust_l[[i]][match(pts_df_l[[i]][["id"]],pts_clust_l[[i]][["id"]]),]
          pts_clust_l[[i]][,"dist"] <- pts_df_l[[i]][,"dist"]
          pts_clust_l[[i]][,"id"] <- NULL
        }

        # 50% of steps include removal of points
        if(rm==TRUE) {
            n_steps_rm <- n_steps*prop_out
            n_steps <- n_steps-n_steps_rm
        }

        if (rm == FALSE | prop_out < 1) {
        # define probabilities based on n_steps (without removal)
        pmin <- 1/k
        probs <- lapply(pmin, function(pmin) {
          pr <- seq(pmin, 1, length.out=n_steps) |>
            as.data.frame() |>
            t()
        }) }

        if (rm == TRUE & prop_out == 1) {
          # define probabilities based on n_steps (without removal)
          pmin <- 1/k
          probs <- lapply(pmin, function(pmin) {
            pr <- seq(pmin, 1, length.out=n_steps_rm) |>
              as.data.frame() |>
              t()
          }) }

        # sample folds based on those probabilities
        pts_df <- sf::st_drop_geometry(pts_clust)
        pts_df <- as.data.frame(pts_df[,startsWith(names(pts_df), "groups")])
        pts_f <- lapply(1:length(pmin), function(tpoints) tpoints=pts_df)

        for(i in 1:ncol(pts_df)) {
          for(j in 1:nrow(pts_df)) {

            # use the initial groups as possible folds
            all_folds = unique(pts_df[,i])
            current_fold = pts_df[j,i]
            other_folds = all_folds[all_folds != current_fold]

            for (p in 1:length(probs[[i]])) {

              # use the defined probabilities
              current_prob = probs[[i]][p]
              other_probs =  rep((1-current_prob)/(length(other_folds)),
                                 length(other_folds))

              # sample folds
              pts_f[[i]][j,paste0("fold_", p)] <- sample(c(current_fold, other_folds),
                                                         size = 1,
                                                         prob = c(current_prob, other_probs))
            }
          }
          pts_f[[i]] <- pts_f[[i]][,grepl("fold_", colnames(pts_f[[i]]))]
          if (length(k) == 1) {
            colnames(pts_f[[i]]) <- paste0(k,"_fold_p_", 1:p)
          } else {
            colnames(pts_f[[i]]) <- paste0(k[[i]],"_fold_p_", 1:p)
          }
        }

      if(rm==TRUE) {
        # sample with removal
        rm_prop_start <- ifelse(1/n_steps_rm < 0.5, 1/n_steps_rm, 0.1)
        rem_st <- seq(rm_prop_start,max_pout,length.out=n_steps_rm)

        pts_dist <- lapply(pts_clust_l,sf::st_drop_geometry)
        pts_out <- pts_dist

        for (i in seq_along(pts_dist)) {

          for (r in seq_along(rem_st)) {
            coln <- paste0(k[[i]], "_fold_r" ,r)
            pts_out[[i]][,coln] <- pts_out[[i]][,paste0("groups_", k[[i]])]

            for (p in 1:k[[i]]) {
              row_ind_group <- which(pts_dist[[i]][paste0("groups_",k[[i]])]==p)
              nrow_p <- length(row_ind_group)
              nrem <- as.integer(nrow_p*rem_st[[r]])
              pts_out[[i]][row_ind_group,"id"] <- 1:nrow_p
              id_rem <- which(pts_out[[i]][,"id"] < nrem)

              pts_out[[i]][id_rem, coln]  <- 0

              pts_out[[i]][,"id"] <- NULL


            }

          }

          pts_out[[i]] <- pts_out[[i]][,grepl("_fold_r", colnames(pts_out[[i]]))]
          pts_out[[i]][,"ID"] <- as.numeric(rownames(pts_out[[i]]))
          pts_out[[i]] <- pts_out[[i]][order(pts_out[[i]][,"ID"]),]

        }
      }
    }

      # collapse the folds and with colnames: folds _ number of probability step
      if(prop_out < 1 & rm == TRUE) {
        pts_id_w <- do.call(cbind.data.frame, pts_f)
        pts_rm <- do.call(cbind.data.frame, pts_out)
        pts_id_w <- cbind(pts_id_w, pts_rm)
      } else if (rm==TRUE) {
        pts_id_w <- do.call(cbind.data.frame, pts_out)
      } else {
        pts_id_w <- do.call(cbind.data.frame, pts_df)
      }
      pts_id_w <- pts_id_w[,!names(pts_id_w) %in% "ID"]
      pts_id_l <- as.list(pts_id_w)
      indx_removed <- lapply(pts_id_l, function(y) {
        if(0 %in% y)  which(0 == y)
        else NULL
      })
      pts_id_l <- lapply(pts_id_l, function(y) y[y!=0])

    }


    # calculate NN distance of CV-folds (Gjstar) ---------------------------------------

    # loop over folds to calculate CV-distance for each block size
    Gjstar <- lapply(seq_along(pts_id_l), function(i){
      if (n_steps==1) distclust(distmat, pts_id_l[[i]])
      else {
        if (length(indx_removed[[i]] > 0)) {
          distclust(distmat[-indx_removed[[i]],-indx_removed[[i]]], pts_id_l[[i]])}
        else distclust(distmat, pts_id_l[[i]])
      }})

    # calculate the WS.distance of Gjstar to Gij for each CV-fold design
    ws_dist <- lapply(Gjstar, function(y) twosamples::wass_stat(Gij, y)) |>
      unlist()
    ws_dist_min <- min(ws_dist)
    ws_dist_min_indx <- which(ws_dist==min(ws_dist))[[1]]

    # get indices of removed points
    if(rm==TRUE & length(indx_removed[[ws_dist_min_indx]]>1)) {
      indx_removed_sel <- indx_removed[[ws_dist_min_indx]]
    } else indx_removed_sel <- integer(0)

    # select the CV-fold design yielding the lowest WS.distance
    fold_min_dist <- names(pts_id_l[ws_dist_min_indx])[[1]] # the colname of the selected CV design
    if(rm == FALSE & n_steps>1) {
      message(paste(fold_min_dist, "was selected"))
      } else if(rm==TRUE) {
      if(n_steps_rm > 1) message(paste(fold_min_dist, "was selected"))
    } else message(paste(max(unique(blocks)), "folds were selected"))

    Gjstar <- Gjstar[[ws_dist_min_indx]] # the Gjstar of the selected design
    fold_cl <- data.frame(cl=pts_id_l[[ws_dist_min_indx]])
    cfolds <- CAST::CreateSpacetimeFolds(fold_cl, spacevar="cl", k=length(unique(fold_cl$cl)))

    if(n_steps==1) {
      blocks <- pts_id_l[[ws_dist_min_indx]]
      res <- list(blocks = blocks,
                  indx_train = cfolds$index, indx_test = cfolds$indexOut,
                  Gij = Gij, Gj = Gj, Gjstar = Gjstar, stat = ws_dist_min)
      class(res) <- c("knndm", "list")
      res
    }else  {
      blocks <- pts_id_w[,ws_dist_min_indx] # the assignment of tpoints to CV-folds
      # Output
      res <- list(blocks = blocks,
                  indx_train = cfolds$index, indx_test = cfolds$indexOut, indx_rem=indx_removed_sel,
                  Gij = Gij, Gj = Gj, Gjstar = Gjstar, stat = ws_dist_min)
      class(res) <- c("knndm", "list")
      res
    }
  }
}




# Helper function: Compute out-of-fold NN distance
distclust <- function(distm, folds){

  alldist <- c()
  for(f in unique(folds)){
    alldist <- c(alldist, apply(distm[f == folds, f != folds, drop=FALSE], 1, min))
  }
  alldist
}
