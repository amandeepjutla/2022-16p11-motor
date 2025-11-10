# _functions.r 
# created 20220826
# last updated 20221002

# functions
clean_df <- function(df_to_clean) {
    cleaned_df <- 
        unique(df_to_clean) |>
        group_by(individual) |>
        mutate(duplicate = n()) |>
        filter(duplicate == 1) |>
        select(-duplicate) |>
        ungroup()
    return(cleaned_df)
}

load_data <- function(path, file_type = 'csv') {
    extension <- paste0('*.', file_type)    
    file_names <- list.files(path, pattern=extension)
    
    df <- file_names |>
        map(function(x) str_c(path, '/', x)) |>
        map(fread) |>
        map(as_tibble)
    
    names(df) <- file_names |>
        gsub(pattern = extension, replacement= '') |>
        gsub(pattern = '-', replacement = '_')

    return(df)
}
            
strict_left_join <- function(x, y, by = NULL, ...) {
    by <- common_by(by, x, y)
    if(any(duplicated(y[by$y]))) {
        stop("Duplicate values in foreign key")
    } else left_join(x, y, by = by, ...)
}
            
summarize_continuous_variable <- function(data, variable) {
    summary <- data |> summarize(
        count = n(),
        mean = round(mean({{variable}}, na.rm = TRUE), 2),
        sd = round(sd({{variable}}, na.rm = TRUE), 2),
        variable = deparse(quote({{variable}})) # probably a better way to do this
    )
    if(!is.grouped_df(data)) {
        summary <- add_column(summary, "group" = "all", .before = 1)
    }
    return(summary)
}
       
summarize_discrete_variable <- function(data, variable) {
    my_variable <- enquo(variable)
    summary <- data |> count(!!my_variable) |> mutate(percent = round(n/sum(n)*100, 2))
    
    if(!is.grouped_df(data)) {
        summary <- add_column(summary, "group" = "all", .before = 1)
    }
    return(summary)
}
            
# This function is adapted from Dustin Fife's "fifer" package
# at https://github.com/dustinfife/fifer
chisq.post.hoc <- function(
    tbl,
    test=c("fisher.test"), 
    popsInRows=TRUE,
    control=c("fdr","BH","BY","bonferroni","holm","hochberg","hommel"),
    digits=4, ...) {
        control <- match.arg(control) # extract correction method
        test = match.fun(test) # extract which test (fisher or chi square)

        if (!popsInRows) tbl <- t(tbl)  # test rows or columns
        popsNames <- rownames(tbl)
        
        prs <- combn(1:nrow(tbl),2) # come up with all possible comparisons
        
        tests <- ncol(prs) # preallocate
        pvals <- numeric(tests)
        lbls <- character(tests)
        
        for (i in 1:tests) {
            pvals[i] <- test(tbl[prs[,i],], ...)$p.value
            lbls[i] <- paste(popsNames[prs[,i]],collapse=" vs. ")
        }
    
        adj.pvals <- p.adjust(pvals,method=control)
        cat("Adjusted p-values used the",control,"method.\n\n")
        data.frame(
            comparison=lbls,
            raw.p=round(pvals,digits),
            adj.p=round(adj.pvals,digits))
}

# This function is adapted from Hause Lin's "hausekeep" package 
# at https://github.com/hauselin/hausekeep

#' @title Identify outliers using robust median absolute deviation approach
#' @name outliersMAD
#'
#' @description identify outliers in vectors using Leys et al.'s (2013) median absolute deviation approach.
#'
#' @param x a vector of numbers
#' @param MADCutOff value to use as cutoff (Leys e tal. recommend 2.5 or 3.0 as default)
#' @param replaceOutliersWith if value is an outlier, what to replace it with? NA by default
#' @param showMADValues if TRUE, will show deviation score of each value
#' @param outlierIndices return index/position of outlier
#' @param bConstant b = 1/qnorm(0.75) (1.4826 if data are normally distributed)
#' @param digits how many digits/decimals to round output to
#'
#' @return A vector with outliers identified (default converts outliers to NA)
#'
#' @references \itemize{
#' \item Leys, C., Ley, C., Klein, O., Bernard, P., & Licata, L. (2013). Detecting outliers: Do not use standard deviation around the mean, use absolute deviation around the median. Journal of Experimental Social Psychology, 49(4), 764-766. doi:10.1016/j.jesp.2013.03.013 (\url{https://www.sciencedirect.com/science/article/pii/S0022103113000668})}
#'
#' @author Hause Lin
#'
#' @export
#'
#' @usage
#' outliersMAD(x, MADCutOff = 3.0, replaceOutliersWith = NA,
#' showMADValues = FALSE, outlierIndices = FALSE, bConstant = 1.4826, digits = 2)
#'
#' @examples
#' example <- c(1, 3, 3, 6, 8, 10, 10, 1000, -1000) # 1000 is an outlier
#' outliersMAD(example)
#' outliersMAD(example, MADCutOff = 3.0)
#' outliersMAD(example, MADCutOff = 2.5, replaceOutliersWith = -999)
#' outliersMAD(example, MADCutOff = 1.5, outlierIndices = TRUE)
#' outliersMAD(example, MADCutOff = 1.5, showMADValues = TRUE)
#' outliersMAD(example, MADCutOff = 1.5, showMADValues = TRUE, replaceOutliersWith = -88)
outliersMAD <- function(x, MADCutOff = 3.0, replaceOutliersWith = NA, showMADValues = FALSE, outlierIndices = FALSE, bConstant = 1.4826, digits = 2) {
  # bConstant: usually, b = 1.4826, a constant linked to the assumption of normality of the data, disregarding the abnormality induced by out- liers (Rousseeuw & Croux, 1993).

  # compute number of absolute MADs away for each value: formula: abs( ( x - median(x) ) )/ mad(x)
  MADAway <- (x - stats::median(x, na.rm = T)) / stats::mad(x, constant = bConstant, na.rm = T)
  absMADAway <- abs(MADAway)
  # subset data that has absMADAway greater than the MADCutOff and replace them with replace
  x[absMADAway > MADCutOff] <- replaceOutliersWith
  outliers <- length(x[absMADAway > MADCutOff])
  if (showMADValues) { # if values == TRUE, return number of mads for each value
    message("Showing MAD from median for each value.")
    message(paste0(outliers, " outliers detected."))
    return(round(MADAway, digits))
  } else if (outlierIndices) {
    message("Showing indices of outliers.")
    if (is.na(replaceOutliersWith)) {
      return(which(is.na(x)))
    } else {
      return(x[x == replaceOutliersWith])
    }

  } else {
    message(paste0(outliers, " outliers detected."))
    message(paste0("Outliers replaced with ", replaceOutliersWith))
    return(round(x, digits)) # otherwise, return original with outliers replaced
  }
}

#' Detect outliers using the modified Z score method of Iglewicz and Hoaglin
#' 
#' Identify outliers within a distribution of numeric values using
#' the modified Z score method. Iglewicz and Hoaglin recommend an absolute
#' Z score threshold of 3.5 to identify potential outliers.
#' 
#' Full details are provided in: 
#' 
#' Boris Iglewicz and David Hoaglin (1993), 
#' "Volume 16: How to Detect and Handle Outliers", 
#' The ASQC Basic References in Quality Control: Statistical Techniques, 
#' Edward F. Mykytka, ed.
#' 
#' @param x distribution to find outliers in
#' @param threshold absolute value of the modified z score threshold above which
#'   to consider a value an outlier; defaults to \code{3.5} on the 
#'   recommendation of Iglewicz and Hoaglin
#' @param return_scores optionally, return the modified z score of each 
#'   observation instead of a masked version of the input vector
#' 
#' @return a vector of the same length as the input, with outliers masked 
#'   (or, if \code{return_scores} is true, the modified z scores of each 
#'   observation)
#' 
#' @export
iglewicz_hoaglin = function(x, threshold = 3.5, return_scores = F) {
  # check input
  if (!is.numeric(x)) 
    stop("could not identify outliers: input is not numeric")  
  # calculate modified Z scores
  med = median(x, na.rm = T)
  MAD = median(abs(x - med), na.rm = T)
  Mi = 0.6745 * (x - med) / MAD
  if (return_scores) {
    return(Mi)
  } else {
    # mask vector
    x[abs(Mi) > threshold] = NA
    return(x)
  }
}            

message("Loaded functions")