# ----------------------------
# Import packages
# ----------------------------

load_packages <- function(packages) {
    for (pkg in packages) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
            install.packages(pkg, repos = "https://cloud.r-project.org")
        }
        library(pkg, character.only = TRUE)
    }
}

# load all required packages
load_packages(c("metafor", "ggplot2", "dplyr", "ggstance"))


# ----------------------------
# Data handling helpers
# ----------------------------
#' Generic sensitivity analysis framework
#' @param dat Data frame
#' @param mods_formula Model formula
#' @param random_structure Random effects structure
#' @param filters List of filter functions to apply
#' @param grouping_var Optional grouping variable
#' @param plot Whether to plot results
#' @param plot_func Plotting function to use
sensitivity_analysis_framework <- function(
    dat,
    mods_formula,
    random_structure,
    filters,
    grouping_var = NULL,
    plot = TRUE,
    plot_func = plot_extreme_filtering) {
    results_list <- list()
    model_info_list <- list()

    # determine if we're doing grouped or ungrouped analysis
    if (is.null(grouping_var)) {
        groups <- list("All" = dat)
    } else {
        groups <- split(dat, dat[[grouping_var]])
    }
    # order groups by size (largest first)
    groups <- groups[order(sapply(groups, nrow), decreasing = TRUE)]

    # process each group
    for (group_name in names(groups)) {
        group_data <- groups[[group_name]]

        # skip groups that are too small for analysis
        if (nrow(group_data) < 5) {
            message(sprintf("Skipping group %s: too few observations (%d)", group_name, nrow(group_data)))
            next
        }

        message(sprintf("Processing group: %s", group_name))

        # fit baseline model
        baseline_model <- safe_fit_model(
            group_data, mods_formula, random_structure,
            sprintf("%s baseline", group_name)
        )
        if (!is.null(baseline_model)) {
            group_label <- if (is.null(grouping_var)) NA else group_name
            results_list[[length(results_list) + 1]] <- extract_results(
                baseline_model, "All data", "None (all data)", NA, group_label, group_data
            )
            model_info_list[[paste0(group_name, "_baseline")]] <- list(
                model = baseline_model, data = group_data
            )
        }

        # apply each filter
        for (filter_name in names(filters)) {
            filter_func <- filters[[filter_name]]
            filter_results <- filter_func(
                group_data, group_name, baseline_model,
                mods_formula, random_structure
            )

            # add results from this filter
            if (!is.null(filter_results$results)) {
                results_list <- c(results_list, filter_results$results)
            }
            if (!is.null(filter_results$models)) {
                model_info_list <- c(model_info_list, filter_results$models)
            }
        }
    }

    # combine and process results
    if (length(results_list) > 0) {
        # ensure ExclusionValue is character for consistency
        results_list <- lapply(results_list, function(df) {
            if ("ExclusionValue" %in% names(df)) {
                df$ExclusionValue <- as.character(df$ExclusionValue)
            }
            df
        })

        sens_results <- bind_rows(results_list)

        if (plot && nrow(sens_results) > 0) {
            plot_func(sens_results)
        }

        return(list(
            results = sens_results,
            models = model_info_list
        ))
    } else {
        warning("No successful model fits")
        return(list(results = data.frame(), models = list()))
    }
}



#' Standard result extraction from fitted models
#' @param model Fitted rma.mv model
#' @param scenario Scenario description
#' @param exclusion_type Type of exclusion applied
#' @param exclusion_value Value/level of exclusion
#' @param group Optional group identifier
#' @param dat Original data used to fit the model (needed for scaling info)
extract_results <- function(model, scenario, exclusion_type = NA, exclusion_value = NA, group = NA, dat = NULL) {
    coefs <- coef(summary(model))

    # create base result dataframe
    return(data.frame(
        Scenario = scenario,
        ExclusionType = exclusion_type,
        ExclusionValue = exclusion_value,
        Group = group,
        Term = rownames(coefs),
        est = coefs[, "estimate"],
        ci.lb = coefs[, "ci.lb"],
        ci.ub = coefs[, "ci.ub"],
        QE = model$QE,
        QEp = model$QEp,
        QM = model$QM,
        QMp = model$QMp,
        k = model$k,
        stringsAsFactors = FALSE
    ))
}

#' Safe model fitting with error handling
#' @param dat Data frame
#' @param mods_formula Model formula
#' @param random_structure Random effects structure
#' @param scenario_name Name for logging
safe_fit_model <- function(dat, mods_formula, random_structure, scenario_name = "model") {
    if (nrow(dat) < 5) {
        warning(sprintf("Too few rows (%d) for %s, skipping.", nrow(dat), scenario_name))
        return(NULL)
    }

    tryCatch(
        rma.mv(yi, vi,
            mods = mods_formula,
            random = random_structure,
            data = dat, method = "REML"
        ),
        error = function(e) {
            message(sprintf("Model failed for %s: %s", scenario_name, e$message))
            return(NULL)
        }
    )
}


#' Compute Cook's Distance for metafor models
#' @param fit A fitted rma.mv object
#' @param parallel Whether to use parallel processing
#' @param ncpus Number of CPU cores to use
#' @return Vector of Cook's distances
compute_cooks_distance <- function(fit, parallel = "multicore", ncpus = 64) {
    cd <- cooks.distance(fit,
        progbar = TRUE,
        reestimate = FALSE,
        parallel = parallel,
        ncpus = ncpus,
        cl = NULL
    )
    # return as datafram
    return(data.frame(cooks.distance = cd))
}


# ----------------------------
# Filtering
# ----------------------------
#' Extreme value percentile filter
#' @param dat Data frame
#' @param mods_formula Model formula
#' @param random_structure Random effects structure
#' @param plot Whether to plot results
#' @param percentiles Vector of percentiles for extreme value trimming
#' @param extreme_limits Vector of length 2 specifying the lower and upper limits for the effect size
extreme_filtering <- function(
    dat,
    mods_formula = NULL,
    random_structure,
    plot = TRUE,
    percentiles = 1:10,
    extreme_limits = NULL) {
    # define filters to apply
    filters <- list()

    # add percentile filters with closures to capture parameters
    filters[["extreme_percentile"]] <- function(group_data, group_name, baseline_model, mods_formula, random_structure) {
        extreme_percentile_filter(group_data, group_name, baseline_model, mods_formula, random_structure, percentiles)
    }

    filters[["variance_percentile"]] <- function(group_data, group_name, baseline_model, mods_formula, random_structure) {
        variance_percentile_filter(group_data, group_name, baseline_model, mods_formula, random_structure, percentiles)
    }

    # add user limits filter if provided
    if (!is.null(extreme_limits)) {
        filters[["user_limits"]] <- function(group_data, group_name, baseline_model, mods_formula, random_structure) {
            user_limits_filter(group_data, group_name, baseline_model, mods_formula, random_structure, extreme_limits)
        }
    }

    # create plotting function
    plot_func <- if (plot) {
        function(sens_results) {
            plot_extreme_filtering(sens_results)
        }
    } else {
        NULL
    }

    # use the framework
    return(sensitivity_analysis_framework(
        dat = dat,
        mods_formula = mods_formula,
        random_structure = random_structure,
        filters = filters,
        grouping_var = NULL,
        plot = plot,
        plot_func = plot_func
    ))
}


#' Extreme value percentile filter
#' @param group_data Data frame containing meta-analysis data for a specific group
#' @param group_name Name of the group
#' @param baseline_model Fitted rma.mv model for the baseline model
#' @param mods_formula Model formula for fixed effects
#' @param random_structure Random effects structure
#' @param percentiles Vector of percentiles for extreme value trimming
#' @return List containing the results and model information
extreme_percentile_filter <- function(group_data, group_name, baseline_model, mods_formula, random_structure, percentiles = 1:10) {
    results_list <- list()
    model_info_list <- list()

    for (p in percentiles) {
        q_low <- quantile(group_data$yi, p / 100, na.rm = TRUE)
        q_high <- quantile(group_data$yi, 1 - p / 100, na.rm = TRUE)
        dat_trim <- subset(group_data, yi >= q_low & yi <= q_high)
        scenario <- sprintf("Trim extremes (%.0f%%)", p)

        print(sprintf("Fitting model with extremes trimmed at %.0f%% (%.2f, %.2f)...", p, q_low, q_high))

        res_trim <- safe_fit_model(dat_trim, mods_formula, random_structure, scenario)
        if (!is.null(res_trim)) {
            results_list[[length(results_list) + 1]] <- extract_results(
                res_trim, scenario, "Effect size percentile", p, group_name, group_data
            )
            model_info_list[[paste0(group_name, "_trim_", p)]] <- list(
                model = res_trim, data = dat_trim, q_low = q_low, q_high = q_high
            )
        }
    }

    return(list(results = results_list, models = model_info_list))
}

#' Variance percentile filter
#' @param group_data Data frame containing meta-analysis data for a specific group
#' @param group_name Name of the group
#' @param baseline_model Fitted rma.mv model for the baseline model
#' @param mods_formula Model formula for fixed effects
#' @param random_structure Random effects structure
#' @param percentiles Vector of percentiles for variance trimming
#' @return List containing the results and model information
variance_percentile_filter <- function(group_data, group_name, baseline_model, mods_formula, random_structure, percentiles = 1:10) {
    results_list <- list()
    model_info_list <- list()

    for (p in percentiles) {
        v_thresh <- quantile(group_data$vi, 1 - p / 100, na.rm = TRUE)
        dat_var <- subset(group_data, vi <= v_thresh)
        scenario <- sprintf("High variance (top %.0f%%) excluded", p)

        print(sprintf("Fitting model with high-variance points excluded at %.0f%% (v_thresh=%.2f)...", p, v_thresh))

        res_var <- safe_fit_model(dat_var, mods_formula, random_structure, scenario)
        if (!is.null(res_var)) {
            results_list[[length(results_list) + 1]] <- extract_results(
                res_var, scenario, "Variance percentile", p, group_name, group_data
            )
            model_info_list[[paste0(group_name, "_var_", p)]] <- list(
                model = res_var, data = dat_var, v_thresh = v_thresh
            )
        }
    }

    return(list(results = results_list, models = model_info_list))
}

#' User-defined extreme limits filter
#' @param group_data Data frame containing meta-analysis data for a specific group
#' @param group_name Name of the group
#' @param baseline_model Fitted rma.mv model for the baseline model
#' @param mods_formula Model formula for fixed effects
#' @param random_structure Random effects structure
#' @param extreme_limits Vector of length 2 specifying the lower and upper limits for the effect size
#' @return List containing the results and model information
user_limits_filter <- function(group_data, group_name, baseline_model, mods_formula, random_structure, extreme_limits = NULL) {
    results_list <- list()
    model_info_list <- list()

    if (!is.null(extreme_limits)) {
        # accept either a list of limits or a single vector
        if (is.list(extreme_limits)) {
            limits_list <- extreme_limits
        } else if (is.numeric(extreme_limits) && length(extreme_limits) == 2) {
            limits_list <- list(extreme_limits)
        } else {
            stop("extreme_limits must be a list of length-2 numeric vectors or a single length-2 numeric vector.")
        }

        for (lims in limits_list) {
            lim_low <- lims[1]
            lim_high <- lims[2]
            dat_lim <- subset(group_data, yi >= lim_low & yi <= lim_high)
            scenario <- sprintf("User limits (%.2f, %.2f)", lim_low, lim_high)

            print(sprintf("Fitting model with user-specified limits: yi in [%.2f, %.2f]...", lim_low, lim_high))

            res_lim <- safe_fit_model(dat_lim, mods_formula, random_structure, scenario)
            if (!is.null(res_lim)) {
                exclusion_value_str <- sprintf("%.2f:%.2f", lim_low, lim_high)
                userlimits_str <- sprintf("User limits (%.1f, %.1f)", lim_low, lim_high)

                results_list[[length(results_list) + 1]] <- extract_results(
                    res_lim, scenario, userlimits_str, exclusion_value_str, group_name, group_data
                )
                model_info_list[[paste0(group_name, "_limits_", length(results_list))]] <- list(
                    model = res_lim, data = dat_lim, lim_low = lim_low, lim_high = lim_high
                )
            }
        }
    }

    return(list(results = results_list, models = model_info_list))
}


#' Filter data based on dependent variables and  climate bounds
#' @param group_data Data frame containing meta-analysis data for a specific group
#' @param group_name Name of the group
#' @param baseline_model Fitted rma.mv model for the baseline model
#' @param mods_formula Model formula for fixed effects
#' @param random_structure Random effects structure
#' @param dt_bounds Vector of length 2 specifying the lower and upper bounds for the temperature moderator
#' @param dph_bounds Vector of length 2 specifying the lower and upper bounds for the pH moderator
climate_filter <- function(group_data, group_name, baseline_model, mods_formula, random_structure,
                           dt_bounds = c(0, 4), dph_bounds = c(-0.4, 0)) {
    results_list <- list()
    model_info_list <- list()

    print("Fitting climate filtered model...")
    dat_clim <- subset(
        group_data,
        dt >= dt_bounds[1] & dt <= dt_bounds[2] &
            dph >= dph_bounds[1] & dph <= dph_bounds[2]
    )

    res_clim <- safe_fit_model(dat_clim, mods_formula, random_structure, "Climate filtered")
    if (!is.null(res_clim)) {
        results_list[[length(results_list) + 1]] <- extract_results(
            res_clim, "Climate filtered", "Climate bounds",
            sprintf("dt[%.1f,%.1f],dph[%.1f,%.1f]", dt_bounds[1], dt_bounds[2], dph_bounds[1], dph_bounds[2]),
            group_name, group_data
        )
        model_info_list[[paste0(group_name, "_clim_filtered")]] <- list(
            model = res_clim, data = dat_clim
        )
    }

    return(list(results = results_list, models = model_info_list))
}


#' Climate bounds filter
#' @param dat Data frame containing meta-analysis data
#' @param mods_formula Model formula for fixed effects
#' @param random_structure Random effects structure
#' @param grouping_var Grouping variable
#' @param dt_bounds Vector of length 2 specifying the lower and upper bounds for the temperature moderator
#' @param dph_bounds Vector of length 2 specifying the lower and upper bounds for the pH moderator
#' @param plot Logical. If TRUE, plots the results
clim_filtering <- function(dat,
                           mods_formula = NULL,
                           random_structure,
                           grouping_var = "core_grouping",
                           dt_bounds = c(-2, 4),
                           dph_bounds = c(-0.4, 0),
                           plot = TRUE) {
    filters <- list()
    filters[["climate"]] <- function(group_data, group_name, baseline_model, mods_formula, random_structure) {
        climate_filter(group_data, group_name, baseline_model, mods_formula, random_structure, dt_bounds, dph_bounds)
    }

    plot_func <- if (plot) {
        function(sens_results) {
            plot_coefficient_comparison(sens_results)
        }
    } else {
        NULL
    }

    return(sensitivity_analysis_framework(
        dat = dat,
        mods_formula = mods_formula,
        random_structure = random_structure,
        filters = filters,
        grouping_var = grouping_var,
        plot = plot,
        plot_func = plot_func
    ))
}


#' Cook's distance influence filter
#' @param group_data Data frame containing meta-analysis data for a specific group
#' @param group_name Name of the group
#' @param baseline_model Fitted rma.mv model for the baseline model
#' @param mods_formula Model formula for fixed effects
#' @param random_structure Random effects structure
#' @return List containing the results and model information
cooks_distance_filter <- function(group_data, group_name, baseline_model, mods_formula, random_structure) {
    results_list <- list()
    model_info_list <- list()

    if (is.null(baseline_model)) {
        return(list(results = list(), models = list()))
    }

    n_obs <- nrow(group_data)
    n_coef <- length(coef(baseline_model))

    # Calculate thresholds
    thresholds <- list(
        "conservative" = 4 / (n_obs - n_coef),
        "liberal" = 2 * sqrt(n_coef / (n_obs - n_coef - 1))
    )

    cooks_scenario_labels <- c(
        "conservative" = "Conservative threshold",
        "liberal" = "Liberal threshold"
    )

    cooks <- compute_cooks_distance(baseline_model)

    for (thr_name in names(thresholds)) {
        thr <- thresholds[[thr_name]]
        print(paste("Applying", thr_name, "threshold:", round(thr, 4)))

        to_exclude <- which(cooks > thr)
        scenario_label <- cooks_scenario_labels[[thr_name]]

        if (length(to_exclude) > 0 && length(to_exclude) < nrow(group_data)) {
            dat_cook <- group_data[-to_exclude, ]
            res_cook <- safe_fit_model(dat_cook, mods_formula, random_structure, scenario_label)

            if (!is.null(res_cook)) {
                results_list[[length(results_list) + 1]] <- extract_results(
                    res_cook, scenario_label, "Cook's distance", thr_name, group_name, group_data
                )
                model_info_list[[paste0(group_name, "_", scenario_label)]] <- list(
                    model = res_cook, data = dat_cook, excluded = to_exclude, cooks_threshold = thr
                )
            }
        } else if (length(to_exclude) == 0) {
            message(sprintf("No influential points found for %s threshold in group %s.", thr_name, group_name))
            results_list[[length(results_list) + 1]] <- extract_results(
                baseline_model, scenario_label, "Cook's distance", thr_name, group_name, group_data
            )
        } else {
            message(sprintf("All points would be excluded for %s threshold in group %s; skipping.", thr_name, group_name))
        }
    }

    return(list(results = results_list, models = model_info_list))
}


# ----------------------------
# Plotting
# ----------------------------

#' Plot coefficient comparison across sensitivity analysis scenarios
#' @param sens_results Data frame containing sensitivity analysis results
#' @param title Plot title
#' @param width Plot width
#' @param height Plot height
plot_coefficient_comparison <- function(sens_results, title = "", width = 8, height = 6, return_plot = FALSE) {
    # reduce vertical spacing between points/bars
    scenario_spacing <- 0.15

    # convert Term to numeric for positioning with tighter spacing
    unique_terms <- unique(sens_results$Term)
    n_terms <- length(unique_terms)

    # use compressed spacing instead of integer spacing
    term_spacing <- 0.8 # Reduce from default 1.0 spacing
    term_positions <- seq(1, by = term_spacing, length.out = n_terms)
    names(term_positions) <- unique_terms

    sens_results$Term_num <- term_positions[sens_results$Term]

    # calculate vertical offsets for each coefficient term within each group
    sens_results$y_offset <- sens_results$Term_num

    for (group in unique(sens_results$Group)) {
        group_data <- sens_results[sens_results$Group == group, ]

        for (term in unique(group_data$Term)) {
            term_data <- group_data[group_data$Term == term, ]
            term_num <- unique(term_data$Term_num)

            # get unique scenarios for this term within this group
            scenarios <- unique(term_data$Scenario)
            n_scenarios <- length(scenarios)

            if (n_scenarios > 1) {
                # center the scenarios around the base term position
                scenario_positions <- seq(-scenario_spacing * (n_scenarios - 1) / 2,
                    scenario_spacing * (n_scenarios - 1) / 2,
                    length.out = n_scenarios
                )
                names(scenario_positions) <- scenarios

                # apply offsets
                for (scenario in scenarios) {
                    row_idx <- which(sens_results$Group == group &
                        sens_results$Term == term &
                        sens_results$Scenario == scenario)

                    sens_results$y_offset[row_idx] <- term_num + scenario_positions[scenario]
                }
            } else {
                # single scenario, just use term position
                row_idx <- which(sens_results$Group == group & sens_results$Term == term)
                sens_results$y_offset[row_idx] <- term_num
            }
        }
    }

    # determine filtering type and create appropriate color scheme and legend
    exclusion_types <- unique(sens_results$ExclusionType)
    scenarios <- unique(sens_results$Scenario)

    # create color palette and legend title based on filtering type
    if ("Cook's distance" %in% exclusion_types) {
        # Cook's distance filtering
        color_values <- c(
            "Conservative threshold" = "#E31A1C", # red
            "Liberal threshold" = "#1F78B4", # blue
            "All data" = "#2CA02C" # green
        )
        legend_title <- "Cook's Distance\nThreshold"
        sens_results$Group_label <- NA_character_
        for (grp in unique(sens_results$Group)) {
            n_full <- unique(sens_results[sens_results$Group == grp & sens_results$Scenario == "All data", ]$k)
            n_liberal <- unique(sens_results[sens_results$Group == grp & sens_results$Scenario == "Liberal threshold", ]$k)
            n_conservative <- unique(sens_results[sens_results$Group == grp & sens_results$Scenario == "Conservative threshold", ]$k)
            label <- sprintf(
                "%s (Total n = %.0f | Liberal filtered n = %.0f (%.1f%%) | Conservative filtered n = %.0f (%.1f%%))",
                grp,
                n_full,
                n_liberal,
                100 * (n_full - n_liberal) / n_full,
                n_conservative,
                100 * (n_full - n_conservative) / n_full
            )
            idx <- which(sens_results$Group == grp)
            sens_results$Group_label[idx] <- label
        }
    } else if (any(c("Effect size percentile", "Variance percentile") %in% exclusion_types)) {
        n_scenarios <- length(scenarios)
        # create discrete colour map using viridis_pal from scales package
        color_values <- scales::viridis_pal(option = "plasma")(n_scenarios)
        names(color_values) <- scenarios
        legend_title <- "Filtering\nLevel"
    } else if ("Climate bounds" %in% exclusion_types) {
        color_values <- c(
            "Climate filtered" = "#FF7F00", # orange
            "All data" = "#2CA02C" # green
        )
        legend_title <- "Climate\nBounds"
        sens_results$Group_label <- NA
        for (grp in unique(sens_results$Group)) {
            n_full <- unique(sens_results[sens_results$Group == grp & sens_results$Scenario == "All data", ]$k)
            n_clim_filtered <- unique(sens_results[sens_results$Group == grp & sens_results$Scenario == "Climate filtered", ]$k)
            idx <- which(sens_results$Group == grp)
            sens_results$Group_label[idx] <- sprintf(
                "%s (Total n = %.0f | Post climate filtering n = %.0f | Removed %.1f%%)",
                grp,
                n_full,
                n_clim_filtered,
                100 * (n_full - n_clim_filtered) / n_full
            )
        }
    } else if ("Scaling comparison" %in% exclusion_types) {
        color_values <- c(
            "Raw model coefficients" = "#1F78B4", # blue
            "Scaled model coefficients" = "#E31A1C" # red
        )
        legend_title <- "Model\nType"
        sens_results$Group_label <- NA_character_
        for (grp in unique(sens_results$Group)) {
            n_obs <- unique(sens_results[sens_results$Group == grp, ]$k)[1]
            idx <- which(sens_results$Group == grp)
            sens_results$Group_label[idx] <- sprintf(
                "%s (n = %.0f)",
                grp,
                n_obs
            )
        }
    } else {
        n_scenarios <- length(scenarios)
        color_values <- scales::brewer_pal(type = "qual", palette = "Set1")(n_scenarios)
        names(color_values) <- scenarios
        legend_title <- "Analysis\nType"
        for (grp in unique(sens_results$Group)) {
            n_full <- unique(sens_results[sens_results$Group == grp & sens_results$Scenario == "All data", ]$k)
            idx <- which(sens_results$Group == grp)
            sens_results$Group_label[idx] <- sprintf(
                "%s (Total n = %.0f)",
                grp,
                n_full
            )
        }
    }

    # ensure all scenarios in data have colors
    missing_scenarios <- setdiff(scenarios, names(color_values))
    if (length(missing_scenarios) > 0) {
        additional_colors <- rainbow(length(missing_scenarios))
        names(additional_colors) <- missing_scenarios
        color_values <- c(color_values, additional_colors)
    }

    # ensure Group_label is a factor with levels in the order of appearance in the data
    sens_results$Group_label <- factor(sens_results$Group_label, levels = unique(sens_results$Group_label))

    p <- ggplot(sens_results, aes(x = est, y = y_offset)) +
        geom_point(aes(color = Scenario), size = 2) +
        geom_errorbar(
            aes(xmin = ci.lb, xmax = ci.ub, color = Scenario),
            width = 0.02
        ) +
        geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
        facet_wrap(~Group_label, ncol = 1, scales = "free_y") +
        scale_y_continuous(
            breaks = unique(sens_results$Term_num),
            labels = unique(sens_results$Term),
            name = "Model Coefficient",
            limits = c(min(sens_results$y_offset) - 0.3, max(sens_results$y_offset) + 0.3),
            expand = expansion(mult = 0.02)
        ) +
        scale_color_manual(
            values = color_values,
            name = legend_title
        ) +
        labs(
            x = "Coefficient Estimate (±95% CI)",
            title = title
        ) +
        theme_minimal(base_size = 12) +
        theme(
            panel.grid.minor.y = element_blank(),
            strip.background = element_rect(fill = "grey90", colour = NA),
            strip.text = element_text(face = "bold"),
            legend.position = "bottom"
        )

    if (return_plot) {
        return(p)
    } else {
        print(p)
    }
}


#' Plot extreme value filtering
#' @param sens_results Data frame with sensitivity results
#' @param title Plot title
plot_extreme_filtering <- function(sens_results, title = "") {
    # create manual vertical offsets for each ExclusionType and ExclusionValue combination

    type_spacing <- 0.25 # spacing between different ExclusionTypes
    value_spacing <- 0.02 # spacing between different ExclusionValues within the same type

    # drop effect size percentile exclusion type
    sens_results <- sens_results[sens_results$ExclusionType != "Effect size percentile", ]

    # convert Term to numeric for positioning
    sens_results$Term_num <- as.numeric(as.factor(sens_results$Term))

    # create offset groups: first by ExclusionType, then by ExclusionValue within type
    sens_results$.group_key <- paste(sens_results$ExclusionType, sens_results$ExclusionValue, sep = "_")

    # calculate vertical offsets
    sens_results$y_offset <- sens_results$Term_num

    # create mapping for y-axis labels (needed for scale_y_continuous)
    term_mapping <- unique(sens_results[, c("Term", "Term_num")])
    term_mapping <- term_mapping[order(term_mapping$Term_num), ]

    for (term in unique(sens_results$Term)) {
        term_data <- sens_results[sens_results$Term == term, ]
        term_num <- unique(term_data$Term_num)

        # get unique ExclusionTypes for this term
        exclusion_types <- unique(term_data$ExclusionType)
        n_types <- length(exclusion_types)

        # center the ExclusionTypes around the base position
        type_positions <- seq(-type_spacing * (n_types - 1) / 2,
            type_spacing * (n_types - 1) / 2,
            length.out = n_types
        )
        names(type_positions) <- exclusion_types

        # for each ExclusionType, offset the ExclusionValues
        for (etype in exclusion_types) {
            type_data <- term_data[term_data$ExclusionType == etype, ]
            exclusion_values <- unique(type_data$ExclusionValue)
            n_values <- length(exclusion_values)

            if (n_values > 1) {
                # create offsets for different exclusion values
                value_positions <- seq(-value_spacing * (n_values - 1) / 2,
                    value_spacing * (n_values - 1) / 2,
                    length.out = n_values
                )
                names(value_positions) <- as.character(exclusion_values)

                # apply offsets
                for (i in 1:nrow(type_data)) {
                    row_idx <- which(sens_results$Term == term &
                        sens_results$ExclusionType == etype &
                        sens_results$ExclusionValue == type_data$ExclusionValue[i])

                    sens_results$y_offset[row_idx] <- term_num +
                        type_positions[etype] +
                        value_positions[as.character(type_data$ExclusionValue[i])]
                }
            } else {
                # single value, just use type position
                row_idx <- which(sens_results$Term == term & sens_results$ExclusionType == etype)
                sens_results$y_offset[row_idx] <- term_num + type_positions[etype]
            }
        }
    }

    p <- ggplot(
        sens_results,
        aes(
            x = est,
            y = y_offset,
            shape = ExclusionType,
            color = as.numeric(ExclusionValue)
        )
    ) +
        geom_point(size = 2) +
        geom_errorbar(
            aes(xmin = ci.lb, xmax = ci.ub),
            width = 0.02
        ) +
        geom_vline(xintercept = 0, linetype = "dashed") +
        scale_color_viridis_c(
            option = "plasma",
            name = "Exclusion Level (top n %)",
            na.value = "grey50"
        ) +
        scale_y_continuous(
            breaks = term_mapping$Term_num,
            labels = as.character(term_mapping$Term),
            name = "Coefficient",
            limits = c(min(sens_results$y_offset) - 0.3, max(sens_results$y_offset) + 0.3),
            expand = expansion(mult = 0.02)
        ) +
        labs(
            x = "Coefficient estimate (±95% CI)",
            title = title,
            shape = "Exclusion Type"
        ) +
        theme_minimal(base_size = 14) +
        theme(
            panel.grid.minor.y = element_blank()
        )
    print(p)
}
