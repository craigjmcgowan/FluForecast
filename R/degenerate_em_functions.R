# These functions are taken from the CMU Delphi Epiforecast package
# https://github.com/cmu-delphi/epiforecast-R

cv_apply_helper = function(train_data, test_data, indexer_list, fn, parallel_dim_i=0L, ...) {
  current_dim_i = length(indexer_list)
  stopifnot(dim(train_data)[current_dim_i] == dim(test_data)[current_dim_i])
  if (current_dim_i == 0L) {
    result = fn(train_data, test_data, ...)
    return (result)
  } else {
    current_dim_size = dim(train_data)[current_dim_i]
    current_dim_seq = seq_len(current_dim_size)
    current_dim_inpnames = dimnames(train_data)[[current_dim_i]]
    indexer_name = names(indexer_list)[current_dim_i]
    indexer_val = indexer_list[[current_dim_i]]
    switch(indexer_name,
           each={
             train_inds = current_dim_seq
             test_inds = train_inds
             current_dim_outnames = current_dim_inpnames
           },
           all={
             train_inds = list(TRUE)
             test_inds = train_inds
             current_dim_outnames = "all"
           },
           smear={
             train_inds = lapply(current_dim_seq, function(center_ind) {
               window_inds = center_ind + indexer_val
               window_inds <- window_inds[1L <= window_inds & window_inds <= current_dim_size]
               if (length(window_inds)==0L) {
                 stop ("smear argument led to selection of no data")
               }
               window_inds
             })
             test_inds = current_dim_seq
             current_dim_outnames = current_dim_inpnames
           },
           ablation={
             train_inds = lapply(current_dim_seq, function(left_out_ind) {
               current_dim_seq[-left_out_ind]
             })
             test_inds = train_inds
             current_dim_outnames = current_dim_inpnames
           },
           subsets={
             train_inds = indexer_val
             test_inds = train_inds
             current_dim_outnames = names(indexer_val)
           },
           loo={
             train_inds = lapply(current_dim_seq, function(train_out_ind) {
               current_dim_seq[-train_out_ind]
             })
             test_inds = current_dim_seq
             current_dim_outnames = current_dim_inpnames
           },
           oneahead={
             test_start_ind = match.single.nonna.integer(indexer_val)
             if (test_start_ind <= 1L || current_dim_size < test_start_ind) {
               stop ("oneahead argument outside the index range for the corresponding dimension, or equal to one, resulting in no training data for the first fold")
             }
             test_inds = test_start_ind-1L + seq_len(current_dim_size-test_start_ind+1L)
             train_inds = lapply(test_inds-1L, seq_len)
             current_dim_outnames = current_dim_inpnames[test_inds]
           },
           loo_oneahead={
             oneahead_start_ind = match.single.nonna.integer(indexer_val)
             if (oneahead_start_ind <= 2L || current_dim_size < oneahead_start_ind) {
               stop ("loo_onahead argument outside the index range for the corresponding dimension, or equal to one or two, resulting in no training data for the first fold")
             }
             loo_test_inds = seq_len(oneahead_start_ind-1L)
             loo_train_inds = lapply(loo_test_inds, function(train_out_ind) {
               loo_test_inds[-train_out_ind]
             })
             oneahead_test_inds = oneahead_start_ind-1L + seq_len(current_dim_size-oneahead_start_ind+1L)
             oneahead_train_inds = lapply(oneahead_test_inds-1L, seq_len)
             train_inds = c(loo_train_inds, oneahead_train_inds)
             test_inds = c(loo_test_inds, oneahead_test_inds)
             current_dim_outnames = current_dim_inpnames[test_inds]
           },
           {
             stop("Unrecognized indexer name.")
           }
    )
    stopifnot(length(train_inds) == length(test_inds))
    current_dim_lapply = if (current_dim_i == parallel_dim_i) {
      print("parallel::mclapply")
      parallel::mclapply
    } else {
      lapply
    }
    subresult.list =
      setNames(current_dim_lapply(seq_along(train_inds), function(indset_i) {
        train_data_inds = as.list(rep(TRUE, length(dim(train_data))))
        if (length(train_inds[[indset_i]])==0L) {
          stop ("Some indexing operation resulted in 0 indices.")
        }
        train_data_inds[[current_dim_i]] <- train_inds[[indset_i]]
        inner_train_data = do.call(`[`, c(list(train_data, drop=FALSE), train_data_inds))
        test_data_inds = as.list(rep(TRUE, length(dim(test_data))))
        test_data_inds[[current_dim_i]] <- test_inds[[indset_i]]
        inner_test_data = do.call(`[`, c(list(test_data, drop=FALSE), test_data_inds))
        cv_apply_helper(inner_train_data, inner_test_data, head(indexer_list,-1L), fn, parallel_dim_i=parallel_dim_i, ...)
      }), current_dim_outnames)
    result = simplify2array(subresult.list)
    return (result)
  }
}


cv_apply = function(data, indexer_list, fn, parallel_dim_i=0L, ...) {
  ## If =data= is 1-D, convert it to an array so it will have non-NULL dim:
  if (is.null(dim(data))) {
    data <- as.array(data)
  }
  if (length(dim(data)) != length(indexer_list)) {
    stop ("Need exactly one indexer_list entry per dimension of data (or exactly 1 entry for vector data).")
  }
  if (is.null(names(indexer_list))) {
    stop ("Indexer types must be specified using names in indexer_list; indexer_list has no names.")
  }
  ## Use recursive cv_apply_helper to compute the entries of the result:
  result = cv_apply_helper(data, data, indexer_list, fn, parallel_dim_i=parallel_dim_i, ...)
  ## --- Adjust the class and dimnames: ---
  ## Make sure the result is an array:
  result <- as.array(result)
  ## Make sure the dimnames are a list, not NULL:
  if (is.null(dimnames(result))) {
    dimnames(result) <- rep(list(NULL), length(dim(result)))
  }
  ## Make sure the dimnames names are a character vector, no NULL:
  if (is.null(names(dimnames(result)))) {
    names(dimnames(result)) <- rep("",length(dim(result)))
  }
  ## If the original input =data= had dimnames names, assign them to the
  ## corresponding dimensions in the result:
  if (!is.null(names(dimnames(data)))) {
    names(dimnames(result))[tail(seq_along(dim(result)), length(dim(data)))] <- names(dimnames(data))
  }
  return (result)
}


degenerate_em_weights = function(distr.cond.lkhds,
                                 init.weights=rep(1/dim(distr.cond.lkhds)[2],dim(distr.cond.lkhds)[2]),
                                 stop.eps = sqrt(.Machine[["double.eps"]])) {
  if (any(init.weights < 1e-10)) {
    stop("All init.weight's must be >=1e-10")
  }
  if (!isTRUE(all.equal(1, sum(init.weights)))) {
    stop("Sum of init.weight's must be all.equal to 1.")
  }
  
  ## Set some constants:
  n.obs = dim(distr.cond.lkhds)[1]
  n.distr = dim(distr.cond.lkhds)[2]
  t.distr.cond.lkhds = t(distr.cond.lkhds) # dim: n.distr x n.obs
  
  ## Set initial values of variables adjusted each step:
  weights = init.weights # length: n.distr
  t.lkhds = init.weights*t.distr.cond.lkhds # dim: n.distr x n.obs
  marginals = colSums(t.lkhds) # length: n.obs
  log.lkhd = mean(log(marginals)) # scalar
  ## log.lkhds = list(log.lkhd)
  if (log.lkhd == -Inf) {
    stop ("All methods assigned a probability of 0 to at least one observed event.")
  } else {
    repeat {
      old.log.lkhd = log.lkhd # scalar
      weights <- colMeans(t(t.lkhds)/marginals) # length: n.distr
      t.lkhds <- weights*t.distr.cond.lkhds # dim: n.distr x n.obs
      marginals <- colSums(t.lkhds) # length: n.obs
      log.lkhd <- mean(log(marginals)) # scalar
      ## xxx inefficient
      ## log.lkhds <- c(log.lkhds,list(log.lkhd))
      stopifnot (log.lkhd >= old.log.lkhd)
      if (log.lkhd-old.log.lkhd <= stop.eps || (log.lkhd-old.log.lkhd)/-log.lkhd <= stop.eps) {
        break
      }
    }
  }
  return (weights)
}


generate_indexer_list_weights = function(component.score.array, indexer.list) {
  component.score.array %>>%
    cv_apply(
      indexer.list,
      function(train, test) {
        if (!identical(dimnames(train)[[5L]], "some log score")) {
          stop ('Weighting routine only supports optimizing for "some log score".')
        }
        instance.method.score.mat =
          R.utils::wrap(train, list(1:5,6L))
        na.counts = rowSums(is.na(instance.method.score.mat))
        if (any(! na.counts %in% c(0L, ncol(instance.method.score.mat)))) {
          stop ("NA appeared in probability matrix, but not from nonexistent EW53.")
        }
        instance.method.score.mat <- instance.method.score.mat[na.counts==0L,,drop=FALSE]
        neginf.counts = rowSums(instance.method.score.mat==-Inf)
        if (any(neginf.counts == ncol(instance.method.score.mat))) {
          print(names(neginf.counts)[neginf.counts==ncol(instance.method.score.mat)])
          stop ("All methods assigned a log score of -Inf for some instance.")
        }
        degenerate.em.weights = degenerate_em_weights(exp(instance.method.score.mat))
        return (degenerate.em.weights)
      },
      parallel_dim_i=1L # use parallelism across seasons (only helps for LOSOCV)
    ) %>>%
    ## --- Fix up dim, dimnames: ---
    {
      d = dim(.)
      dn = dimnames(.)
      ## Remove Model="all" dimension:
      stopifnot(identical(tail(dn, 1L), list(Model="all")))
      new.d = head(d, -1L)
      new.dn = head(dn, -1L)
      ## Call the new, unnamed dim 1 "Model":
      stopifnot(identical(
        head(dn, 1L),
        stats::setNames(dimnames(component.score.array)["Model"], "")))
      names(new.dn)[[1L]] <- "Model"
      dim(.) <- new.d
      dimnames(.) <- new.dn
      .
    } %>>%
    ## --- Convert to data.frame with desired format: ---
    {
      ## Rename dimensions:
      stopifnot(identical(
        names(dimnames(.)),
        c("Model", "Season", "Model Week", "Location", "Target", "Metric")
      ))
      names(dimnames(.)) <-
        c("component_model_id","season","Model Week","location","target","Metric")
      ## Drop ="all" dimensions:
      d = dim(.)
      dn = dimnames(.)
      keep.dim.flags = !sapply(dn, identical, "all")
      dim(.) <- d[keep.dim.flags]
      dimnames(.) <- dn[keep.dim.flags]
      .
    } %>>%
    ## Melt into tibble:
    reshape2::melt(value.name="weight") %>>% tibble::as_tibble() %>>%
    dplyr::mutate_if(is.factor, as.character) %>>%
    {.}
}