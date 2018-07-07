
#' @noRd
check.dimensionality <- function(dim_vec, should_be)
{
  if (all(dim_vec != should_be))
  {
    stop(paste(
      sprintf("Prior should have dimensionalty: (%d x %d).",
              should_be[1], should_be[2]),
      sprintf("You provided a prior of dimension:(%d x %d)",
              dim_vec[1], dim_vec[2]),
      sep="\n"), call.=FALSE)
  }
}