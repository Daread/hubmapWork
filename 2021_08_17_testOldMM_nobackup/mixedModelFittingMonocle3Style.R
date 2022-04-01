
# Helper function to send a vector of fixed and random effects. Returns a string that can be converted 
#    to a formula that's understood by lme4
getMixedModelFormula <- function(inputFixedEffects, inputRandomEffects){
  # Formatted Random Effects
  fixedFormulaPart = paste(inputFixedEffects, collapse = " + ")

  # For now, assume there are random effects (otherwise we'd just use the fit_models code in monocle3...)
  randomFormulaPart = paste(inputRandomEffects, collapse = ") + (1|")
  # Add the front and back of the random intercepts formula
  randomFormulaPart = paste0("(1|", randomFormulaPart, ")")

  fullFormula = paste0("~", fixedFormulaPart, " + ", randomFormulaPart)
  return(fullFormula)
}


# Test whether a matrix is one of our supported sparse matrices
is_sparse_matrix_MM <- function(x){
  class(x) %in% c("dgCMatrix", "dgTMatrix", "lgCMatrix")
}


sparse_apply_MM <- function(Sp_X, MARGIN, FUN, convert_to_dense, ...){
  if (convert_to_dense){
    if (MARGIN == 1){
      Sp_X <- Matrix::t(Sp_X)
      res <- lapply(colnames(Sp_X), function(i, FUN, ...) {
        FUN(as.matrix(Sp_X[,i]), ...)
      }, FUN, ...)
    }else{
      res <- lapply(colnames(Sp_X), function(i, FUN, ...) {
        FUN(as.matrix(Sp_X[,i]), ...)
      }, FUN, ...)
    }
  }else{
    if (MARGIN == 1){
      Sp_X <- Matrix::t(Sp_X)
      res <- lapply(colnames(Sp_X), function(i, FUN, ...) {
        FUN(Sp_X[,i], ...)
      }, FUN, ...)
    }else{
      res <- lapply(colnames(Sp_X), function(i, FUN, ...) {
        FUN(Sp_X[,i], ...)
      }, FUN, ...)
    }
  }

  return(res)

}

smart_es_apply_MM <- function(cds, MARGIN, FUN, convert_to_dense,
                           reduction_method="UMAP", ...) {
  parent <- environment(FUN)
  if (is.null(parent))
    parent <- emptyenv()
  e1 <- new.env(parent=parent)
  coldata_df = as.data.frame(colData(cds))
  tryCatch({
    coldata_df$cluster = clusters(cds, reduction_method)[colnames(cds)]
    coldata_df$partition = partitions(cds, reduction_method)[colnames(cds)]
    coldata_df$pseudotime = pseudotime(cds)
  }, error = function(e) {} )
  Biobase::multiassign(names(as.data.frame(coldata_df)),
                       as.data.frame(coldata_df), envir=e1)
  environment(FUN) <- e1

  if (is_sparse_matrix_MM(SingleCellExperiment::counts(cds))){
    res <- sparse_apply_MM(SingleCellExperiment::counts(cds), MARGIN, FUN, convert_to_dense, ...)
  } else {
    res <- pbapply::pbapply(SingleCellExperiment::counts(cds), MARGIN, FUN, ...)
  }

  if (MARGIN == 1)
  {
    names(res) <- row.names(cds)
  }else{
    names(res) <- colnames(cds)
  }

  res
}




#' @title Helper function for model fitting
#' @description test
#' @param x test
#' @param model_formula_str a formula string specifying the model to fit for
#'   the genes.
#' @param expression_family specifies the family function used for expression
#'   responses
#' @param disp_func test
#' @param verbose Whether to show VGAM errors and warnings. Only valid for
#'   cores = 1.
#' @param ... test
#' @name fit_model_helper_MM
#' @keywords internal
fit_model_helper_MM <- function(x,
                             model_formula_str,
                             expression_family,
                             disp_func = NULL,
                             clean_model = TRUE,
                             verbose = FALSE,
                             ...) {
  print("Running fit_model_helper_MM now")

  model_formula_str <- paste("f_expression", model_formula_str,
                           sep = "")
  orig_x <- x
  # FIXME: should we be using this here?
  # x <- x + pseudocount
  if (expression_family %in% c("negbinomial", "poisson", "zinegbinomial",
                               "zipoisson", "quasipoisson")) {
    x <- x / Size_Factor
    f_expression <- round(x)
  }
  else if (expression_family %in% c("binomial", "gaussian")) {
    f_expression <- x
  }
  else {
    # FIXME: maybe emit a warning or error here instead.
    f_expression <- log10(x)
  }
  f_expression = as.numeric(f_expression)
  model_formula = stats::as.formula(model_formula_str)

  # browser()

  tryCatch({
    if (verbose) messageWrapper = function(expr) { expr }
    else messageWrapper = suppressWarnings


    FM_fit = messageWrapper(switch(expression_family,
                    "negbinomial" = glmer.nb(model_formula,    #varsAsForm, data=inputDF, 
                                              nAGQ=0,
                                            control=glmerControl(optimizer = "nloptwrap"), # Faster optimizer. Not used
                                              # for laplacian or GHQ fitting
                                              verbose=FALSE)
                    # "negbinomial" = MASS::glm.nb(model_formula, epsilon=1e-3,
                    #                              model=FALSE, y=FALSE, ...),
                    # "poisson" = speedglm::speedglm(model_formula,
                    #                                family = stats::poisson(),
                    #                                acc=1e-3, model=FALSE,
                    #                                y=FALSE, ...),
                    # "quasipoisson" = speedglm::speedglm(model_formula,
                    #                                     family = stats::quasipoisson(),
                    #                                     acc=1e-3, model=FALSE,
                    #                                     y=FALSE, ...),
                    # "binomial" = speedglm::speedglm(model_formula,
                    #                                 family = stats::binomial(),
                    #                                 acc=1e-3, model=FALSE,
                    #                                 y=FALSE, ...),
                    # "gaussian" = speedglm::speedglm(model_formula,
                    #                                 family = stats::gaussian(),
                    #                                 acc=1e-3, model=FALSE,
                    #                                 y=FALSE, ...),
                    # "zipoisson" = pscl::zeroinfl(model_formula,
                    #                              dist="poisson", ...),
                    # "zinegbinomial" = pscl::zeroinfl(model_formula,
                    #                                  dist="negbin", ...)
                    ))
    
    FM_summary = summary(FM_fit)

    if (clean_model)
    # browser()
      # FM_fit = clean_model_object(FM_fit) # May want to put this back in later to clean up the fitting output of GLMM
    df = list(model=FM_fit, model_summary=FM_summary)
    df
  }, error = function(e) {
    if (verbose) { print (e) }
    list(model=NA, model_summary=NA)
  })
}

#' Fits a model for each gene in a cell_data_set object.
#'
#' This function fits a generalized linear model for each gene in a
#' cell_data_set. Formulae can be provided to account for additional covariates
#' (e.g. day collected, genotype of cells, media conditions, etc).
#'
#' @param cds The cell_data_set upon which to perform this operation.
#' @param model_formula_str A formula string specifying the model to fit for
#'   the genes.
#' @param expression_family Specifies the family function used for expression
#'   responses. Can be one of "quasipoisson", "negbinomial", "poisson",
#'   "binomial", "gaussian", "zipoisson", or "zinegbinomial". Default is
#'   "quasipoisson".
#' @param reduction_method Which method to use with clusters() and
#'   partitions(). Default is "UMAP".
#' @param cores The number of processor cores to use during fitting.
#' @param clean_model Logical indicating whether to clean the model. Default is
#'   TRUE.
#' @param verbose Logical indicating whether to emit progress messages.
#' @param ... Additional arguments passed to model fitting functions.
#'
#' @return a tibble where the rows are genes and columns are
#'   * id character vector from `rowData(cds)$id`
#'   * gene_short_names character vector from `rowData(cds)$gene_short_names`
#'   * num_cells_expressed int vector from `rowData(cds)$num_cells_expressed`
#'   * gene_id character vector from row.names(rowData(cds))`
#'   * model GLM model list returned by speedglm
#'   * model_summary model summary list returned by `summary(model)`
#'   * status character vector of model fitting status: OK when model converged, otherwise FAIL
#'
#' @export
fit_mixed_models <- function(cds,
                     model_formula_str,
                     fixedAndRandomVars,
                     expression_family = "negbinomial",
                     reduction_method="UMAP",
                     cores = 1,
                     clean_model = TRUE,
                     verbose = FALSE,
                     ...) {

  print("Fitting MM Now")
  model_form <- stats::as.formula(model_formula_str)
  if (!"num_cells_expressed" %in% names(rowData(cds))) {
    cds <- detect_genes(cds)
  }
  coldata_df = colData(cds)
  tryCatch({
    coldata_df$cluster = clusters(cds, reduction_method)[colnames(cds)]
    coldata_df$partition = partitions(cds, reduction_method)[colnames(cds)]
    coldata_df$pseudotime = pseudotime(cds, reduction_method)
  }, error = function(e) {} )

  # Test model formula validity.
  # Notes:
  #  o  allow for formulas in model formula: e.g.,  ~ splines::ns(pseudotime, df=3)
  #  o  watch for NA, NaN, Inf in terms
  #  o  is.na counts NA and NaN
  #  o  splines::ns(pseudotime, df=3) fails if there are Inf values, at least,
  #     which causes model.frame( model_formula, ...) to fail
  #  o  model.frame catches mis-spelled functions
  err_msg <- NULL
  mf_terms <- all.vars(model_form)

  
  for( mf_term in mf_terms )
  {
    if(!( mf_term %in% names(coldata_df)))
    {
      err_msg <- paste0(err_msg,'  \'', mf_term, '\': not in cds\n')
      next
    }
  }


  if(length(err_msg) > 0)
    stop( '\n-- bad fit_mixed_models terms --\n', err_msg )

  # browser()
  print("Getting df subset for fitting")

  tryCatch({
    modelFixedForm = stats::as.formula(paste0("~", paste(fixedAndRandomVars, collapse = " + ")))
    stats::model.frame(modelFixedForm, data=coldata_df)
    # browser()
    # stats::model.frame(model_form, data=coldata_df) # Original from fit_models
  }, error = function( cnd ) {
       info_msg <- ''
       for( mf_term in mf_terms )
       {
         mf_length  <- length(coldata_df[[mf_term]])
         mf_num_inf <- sum(is.infinite(coldata_df[[mf_term]]))
         mf_num_nan <- sum(is.nan(coldata_df[[mf_term]]))
         mf_num_na  <- sum(is.na(coldata_df[[mf_term]]))
         if( mf_num_inf > 0 )
           info_msg <- paste0(info_msg, '  \'', mf_term, '\': ' , mf_num_inf, ' of ', mf_length, ' values are Inf\n')
         if( mf_num_nan > 0 )
           info_msg <- paste0(info_msg, '  \'', mf_term, '\': ' , mf_num_nan, ' of ', mf_length, ' values are NaN\n')
         if( mf_num_na - mf_num_nan > 0 )
           info_msg <- paste0(info_msg, '  \'', mf_term, '\': ' , mf_num_na - mf_num_nan, ' of ', mf_length, ' values are NA\n')
       }
       rm( mf_term, mf_length, mf_num_inf, mf_num_nan, mf_num_na )
       stop (paste0( 'Error in model formula: ', conditionMessage( cnd ), '\n', info_msg ) )
  })

  disp_func <- NULL

  print("Fitting Now")
  if (cores > 1) {
    fits <-
      mc_es_apply(
        cds,
        1,
        fit_model_helper_MM,
        required_packages = c("BiocGenerics", "Biobase", "MASS", "purrr",
                              "pscl", "speedglm", "dplyr", "Matrix"),
        cores = cores,
        reduction_method = reduction_method,
        model_formula_str = model_formula_str,
        expression_family = expression_family ,
        disp_func = disp_func,
        clean_model = clean_model,
        verbose = verbose,
        ...
      )
    fits
  } else{
    fits <- smart_es_apply_MM(
      cds,
      1,
      fit_model_helper_MM,
      convert_to_dense = TRUE,
      model_formula_str = model_formula_str,
      expression_family = expression_family,
      reduction_method = reduction_method,
      disp_func = disp_func,
      clean_model = clean_model,
      verbose = verbose,
      ...
    )
    fits
  }

  # browser()

  rowData(cds)$gene_id <- row.names(rowData(cds))
  fits <- tibble::as_tibble(purrr::transpose(fits))
  M_f <- tibble::as_tibble(rowData(cds))
  M_f <- dplyr::bind_cols(M_f, fits)
  M_f <- M_f %>%
    dplyr::mutate(status = purrr::map(.f = purrr::possibly(
      extract_model_status_helper, NA_real_), .x = model)) %>%
    tidyr::unnest(status)



  return(M_f)
}

extract_model_status_helper <- function(model){
  if (class(model)[1] == "speedglm") {
    status_str <- ifelse(model$convergence, "OK", "FAIL")
    return (status_str)

  } else if (class(model)[1] == "negbin"){
    status_str <- ifelse(model$converged, "OK", "FAIL")
    return (status_str)
  } else if (class(model) == "zeroinfl"){
    status_str <- ifelse(model$converged, "OK", "FAIL")
    return (status_str)
  }else {
    return("FAIL")
  }
}

extract_coefficient_helper = function(model, model_summary,
                                      pseudo_count = 0.01) {
  if (class(model)[1] == "speedglm") {
    coef_mat <- model_summary$coefficients # first row is intercept
    # We need this because some summary methods "format" the coefficients into
    # a factor...
    coef_mat <- apply(coef_mat, 2, function(x) {as.numeric(as.character(x)) })
    row.names(coef_mat) = row.names(model_summary$coefficients)
    colnames(coef_mat) = c('estimate',
                           'std_err',
                           'test_val',
                           'p_value')
    log_eff_over_int = log2((model$family$linkinv(coef_mat[, 1] +
                                                    coef_mat[1, 1]) +
                               pseudo_count) /
                              rep(model$family$linkinv(coef_mat[1, 1]) +
                                    pseudo_count, times = nrow(coef_mat)))
    log_eff_over_int[1] = 0
    coef_mat = tibble::as_tibble(coef_mat, rownames = "term")
    coef_mat$normalized_effect = log_eff_over_int
    coef_mat$model_component = "count"
    return (coef_mat)

  } else if (class(model)[1] == "negbin"){
    coef_mat = model_summary$coefficients # first row is intercept
    # We need this because some summary methods "format" the coefficients into
    # a factor...
    coef_mat = apply(coef_mat, 2, function(x) {as.numeric(as.character(x)) })
    row.names(coef_mat) = row.names(model_summary$coefficients)
    colnames(coef_mat) = c('estimate',
                           'std_err',
                           'test_val',
                           'p_value')
    log_eff_over_int = log2((model$family$linkinv(coef_mat[, 1] +
                                                    coef_mat[1, 1]) +
                               pseudo_count) /
                              rep(model$family$linkinv(coef_mat[1, 1]) +
                                    pseudo_count, times = nrow(coef_mat)))
    log_eff_over_int[1] = 0
    coef_mat = tibble::as_tibble(coef_mat, rownames = "term")
    coef_mat$normalized_effect = log_eff_over_int

    coef_mat$model_component = "count"
    return (coef_mat)
  } else if (class(model) == "zeroinfl"){
    count_coef_mat = model_summary$coefficients$count # first row is intercept
    colnames(count_coef_mat) = c('estimate',
                           'std_err',
                           'test_val',
                           'p_value')
    log_eff_over_int = log2((model$linkinv(count_coef_mat[, 1] +
                                             count_coef_mat[1, 1]) +
                               pseudo_count) /
                              rep(model$linkinv(count_coef_mat[1, 1]) +
                                    pseudo_count,
                                  times = nrow(count_coef_mat)))
    log_eff_over_int[1] = 0
    count_coef_mat = tibble::as_tibble(count_coef_mat, rownames = "term")
    count_coef_mat$normalized_effect = log_eff_over_int
    count_coef_mat$model_component = "count"

    zero_coef_mat = model_summary$coefficients$zero # first row is intercept
    colnames(zero_coef_mat) = c('estimate',
                                 'std_err',
                                 'test_val',
                                 'p_value')
    zero_coef_mat = tibble::as_tibble(zero_coef_mat, rownames = "term")
    zero_coef_mat$normalized_effect = NA
    zero_coef_mat$model_component = "zero"
    coef_mat = dplyr::bind_rows(count_coef_mat, zero_coef_mat)
    return (coef_mat)
  } else {
    coef_mat = matrix(NA_real_, nrow = 1, ncol = 5)
    colnames(coef_mat) = c('estimate',
                           'std_err',
                           'test_val',
                           'p_value',
                           'normalized_effect')
    coef_mat = tibble::as_tibble(coef_mat)
    coef_mat$term = NA_character_
    coef_mat$model_component = NA_character_
    return(coef_mat)
  }
}

#' Extracts a table of coefficients from a tibble containing model objects
#'
#' @param model_tbl A tibble of model objects, generally the output of
#'   \code{\link{fit_mixed_models}}.
#' @importFrom dplyr %>%
#' @export
coefficient_table <- function(model_tbl) {
  M_f = model_tbl %>%
    dplyr::mutate(terms = purrr::map2(.f = purrr::possibly(
      extract_coefficient_helper, NA_real_), .x = model,
      .y = model_summary)) %>%
    tidyr::unnest(terms)
  M_f = M_f %>% dplyr::group_by(model_component, term) %>%
    dplyr::mutate(q_value = stats::p.adjust(p_value)) %>% dplyr::ungroup()
  return(M_f)
}