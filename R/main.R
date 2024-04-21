#' @title Calculate Coefficient of Variation
#' @description Calculates the coefficient of variation (CV) between two numerical vectors.
#' @param a The first numerical vector.
#' @param b The second numerical vector.
#' @return A numeric value representing the coefficient of variation.
compute_covariance <- function(a, b) {
  cov(a, b, use = "complete.obs") / (mean(a, na.rm = T) * mean(b, na.rm = T))
}

#' @title Calculate Validation Statistics with Bootstrapping
#' @description Calculates bootstrapped validation statistics (coefficients of variation)
#'   using user-provided growth and death functions.
#' @param data A data frame containing the variables needed for the calculations.
#' @param species_variable Species variable in the data frame (provided as a string).
#' @param growth_function Expression of the growth function (provided as a string).
#' @param death_function Expression of the death function (provided as a string).
#' @param Nbootstraps The number of bootstrap iterations to perform (default is 500).
#' @param paired A logical flag (TRUE/FALSE) indicating whether the growth and death
#'   covariances should be treated as paired for analysis (default is FALSE).
#' @return  A data frame containing z-scores of the bootstrapped coefficients of variation (CoV)
#'   for growth and death.
#' @examples
#' set.seed(1010)
#' covariance_statistics(
#'   example_data,
#'   species_variable = "prey",
#'   growth_function = "prey",
#'   death_function = "prey*predator + prey^2"
#' )
#' @export
covariance_statistics <- function(data, species_variable,
                                  growth_function,
                                  death_function,
                                  Nbootstraps = 500,
                                  paired = F) {
  z_score <- function(cov_growth, cov_death, paired) {
    if(paired){
      diff <- cov_growth - cov_death
      z_score <- abs(mean(diff))/sd(diff)
    } else{
      z_score <- abs((mean(cov_growth) - mean(cov_death))) /
        sqrt(var(cov_growth) + var(cov_death))
    }
    z_score
  }

  data %>%
    rsample::bootstraps(times = Nbootstraps, apparent = TRUE) %>%
    dplyr::mutate(covariance = purrr::map(splits, function(splits) {
      ts <- rsample::analysis(splits) %>%
        na.omit(.)
      with(
        ts,
        c(
          cov_growth = compute_covariance(
            eval(parse(text = growth_function)),
            eval(parse(text = species_variable))
          ),
          cov_death = compute_covariance(
            eval(parse(text = death_function)),
            eval(parse(text = species_variable))
          )
        )
      )
    })) %>%
    dplyr::select(-splits, -id) %>%
    tidyr::unnest_wider(covariance) %>%
    with(z_score(cov_growth, cov_death, paired))
}

#' @title Calculate Covariance Estimates for Growth and Death
#' @description Estimates covariances between growth/death variables and a specified species variable.
#' @param data A data frame containing the necessary variables.
#' @param species_variable The name of the species variable within the data frame (provided as a string).
#' @param growth_function The name of the growth function (provided as a string).
#' @param death_function The name of the death function (provided as a string).
#' @return A named list containing the estimated covariances:
#'   * cov_growth: Covariance between the growth variable and the species variable.
#'   * cov_death: Covariance between the death variable and the species variable.
#' @examples
#' covariance_estimate(
#'   example_data,
#'   species_variable = "prey",
#'   growth_function = "prey",
#'   death_function = "prey*predator + prey^2"
#' )
#' @export
covariance_estimate <- function(data,
                                species_variable,
                                growth_function,
                                death_function) {
  with(
    data,
    c(
      cov_growth = compute_covariance(
        eval(parse(text = growth_function)),
        eval(parse(text = species_variable))
      ),
      cov_death = compute_covariance(
        eval(parse(text = death_function)),
        eval(parse(text = species_variable))
      )
    )
  )
}

#' @title Visualize Convergence of Simulated Covariance Ratios
#' @description Generates a plot demonstrating the convergence of covariance ratios under
#'   simulated process noise, aiding in analysis of model stability.
#' @param data A data frame containing variables for modeling, including a 'time' column for simulation steps.
#' @param species_variable Species variable in the data frame (provided as a string).
#' @param growth_function Expression of the growth function (provided as a string).
#' @param death_function Expression of the death function (provided as a string).
#' @param Ntrials The number of simulation trials to perform (default is 10).
#' @param process_noise The level of process noise to add during simulations (default is 0.05).
#' @return A ggplot object visualizing the convergence of covariance ratios across different simulation lengths.
#' @import ggplot2
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import purrr
#' @examples
#' plot <- convergence_visualization(
#'   example_data,
#'   species_variable = "prey",
#'   growth_function = "prey",
#'   death_function = "prey*predator + prey^2"
#' )
#' print(plot)
#' @export
convergence_visualization <- function(data, species_variable,
                                      growth_function,
                                      death_function,
                                      Ntrials = 10,
                                      process_noise = 0.05) {
  tibble(end_step = floor(nrow(data) / 10):nrow(data)) %>%
    mutate(ts = map(end_step, ~ data[1:.x, ])) %>%
    expand_grid(label = 1:Ntrials) %>%
    mutate(covariance = map(ts, function(ts) {
      ts <- ts %>%
        select(-time) %>%
        mutate_all(
          ~ . * rnorm(n(), 1, process_noise)
        )
      with(
        ts,
        c(
          cov_growth = compute_covariance(
            eval(parse(text = growth_function)),
            eval(parse(text = species_variable))
          ),
          cov_death = compute_covariance(
            eval(parse(text = death_function)),
            eval(parse(text = species_variable))
          )
        )
      )
    })) %>%
    dplyr::select(-ts) %>%
    tidyr::unnest_wider(covariance) %>%
    dplyr::group_by(end_step) %>%
    dplyr::summarise(
      mean = mean(cov_growth / cov_death),
      sd = sd(cov_growth / cov_death),
      .groups = "drop"
    ) %>%
    ggplot(aes(x = end_step, y = mean)) +
    geom_line() +
    geom_ribbon(aes(ymin = mean - 2 * sd, ymax = mean + 2 * sd), alpha = 0.2) +
    geom_hline(color = "firebrick2", yintercept = 1, linetype = "dashed") +
    theme_bw()
}

#' An example dataset of ecological time series
#' @format ## `example_data`
#' A data frame with 198 row and 3 columns:
#' \describe{
#'   \item{time}{Measurement time}
#'   \item{prey}{Prey abundance}
#'   \item{predator}{Predator abundance}
#' }
#' @source <https://www.nature.com/articles/s41586-019-1857-0>
"example_data"
