#' Fish stable isotope dataset.
#'
#' A dataset containing values for three stable isotopes measured in the muscle tissue of four species of arctic fish. For use in examples.
#'
#' @details This dataset contains \eqn{\delta^{15}N}, \eqn{\delta^{13}C}, and \eqn{delta^{34}S} values for the following fish species:
#' \itemize{
#'   \item ARCS - Arctic Cisco (*Coregonus autumnalis*), \eqn{n = 69}.
#'   \item BDWF - Broad Whitefish (*Coregonus nasus*), \eqn{n = 71}.
#'   \item LKWF - Lake Whitefish (*Coregonus clupeaformis*), \eqn{n = 67}.
#'   \item LSCS - Least Cisco (*Coregonus sardinella*), \eqn{n = 70}
#' }
#'
#' @source Fish were collected between 2007 and 2008 from an estuarine area of the Beaufort Sea, North and West of the Mackenzie Delta at Phillips Bay, Yukon Territory, Canada (69.28 N, 138.49 W).
#'
#' @format A data frame with 278 rows (observations) and 4 columns (species, \eqn{\delta^{15}N}, \eqn{\delta^{13}C}, and \eqn{\delta^{34}S}).
#' @examples
#' head(fish)
#' aggregate(fish[2:4], fish[1], mean)
"fish"
