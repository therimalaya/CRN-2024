## ---- Source script ----
source(here::here("Scripts/Functions.R"))

## ---- Load packages ----
pkgs <- c("ggplot2", "data.table", "purrr", "dplyr", "stringr")
for (pkg in pkgs) require(pkg, character.only = TRUE)

## ---- Local functions ----
get_risktable <- function(fit, vars = c("n.risk", "cum.event", "cum.censor"), combine_outcome = FALSE, y = NULL, ... ) {
  out <- tidy_survfit(fit, ...) %>%
    as.data.table() %>% 
    .[, c("outcome", "time", vars), with = FALSE]
  return(out)
}

## ---- Set defaults ----
ggplot2::theme_set(
  ggthemes::theme_few(base_size = 14) +
    theme(
      panel.grid = element_line(color = "#f0f0f0")
    )
)

## ---- Set color theme ----
pal <- paletteer::palettes_d
col_pal_l <- c(pal$ggsci$light_uchicago[c(4, 5, 7, 9)], "#aaaaaa")
col_pal_d <- c(pal$ggsci$default_uchicago[c(4, 5, 7, 9)], "#aaaaaa")

## ---- Plot commons ----
plot_commons <- function(x) {
  col_pal_l <- pal$ggsci$light_uchicago[c(5, 4, 2)]
  col_pal_d <- pal$ggsci$dark_uchicago[c(5, 4, 2)]
  list(
    ggthemes::theme_few(base_size = 14),
    ggplot2::scale_linetype_manual(
      values = 1:3,
      breaks = c("Melanoma", "Other", "Alive")
    ),
    ggplot2::scale_fill_manual(
      values = col_pal_l,
      breaks = c("Melanoma", "Other", "Alive")
    ),
    ggplot2::scale_color_manual(
      values = col_pal_d,
      breaks = c("Melanoma", "Other", "Alive")
    ),
    scale_x_continuous(
      breaks = seq(0, 40 * 365.241, 5 * 365.241),
      labels = \(x) x / 365.241,
      minor_breaks = seq(0, 40 * 365.241, 2.5 * 365.241)
    ),
    scale_y_continuous(
      breaks = seq(0, 1, 0.2),
      labels = scales::label_percent(suffix = ""),
      limits = c(0, 1)
    ),
    ggplot2::theme(
      legend.position = c(0.005, 0.995),
      legend.justification = c(0, 1),
      legend.direction = "horizontal",
      panel.grid = element_line(color = "#f0f0f0"),
      panel.spacing.y = unit(12, "pt"),
      strip.text = element_text(size = unit(14, "pt")),
      axis.title = element_text(size = unit(16, "pt")),
      axis.text.x = element_text(size = unit(12, "pt")),
      axis.text.y = element_text(size = unit(10, "pt"))
    ),
    labs(
      x = "Follow-up time, years",
      y = "Cumulative incidence of death (%)",
      color = NULL,
      fill = NULL
    )
  )
}


## ---- Count data ----
counts <- readRDS(here::here("Data/Counts.Rds"))

## ---- Plots: BarNS: By Cstage and Tstage ----
plot_data <- readRDS(here::here("Data/bar-cstage.Rds"))
count_data <- counts[["by-sex-tstage-year_cat-clinical_stage"]] %>% 
    copy() %>%
    .[sex != "Overall" & tstage != "Overall" & year_cat == "1983-2019"] %>% 
    .[!clinical_stage %in% c("M", "Overall")] %>% 
    .[, sex := fifelse(sex == "female", "Women", "Men")] %>% 
    .[, N := if_else(N <= 5, "≤5", as.character(N))]

plot_by_cstage <- function(dta, count, fup = NULL, layout = "horizontal") {
  if (!is.null(fup)) {
    plot_data <- dta %>% dplyr::filter(end == fup)
  } else {
    plot_data <- dta
  }
  
  plot_out <- plot_data %>% 
    ggplot(aes(clinical_stage, cns_pp, ymin = lo_cns_pp, ymax = hi_cns_pp)) +
    geom_col(
      aes(
        group = tstage, 
        fill = tstage,
        color = tstage
      ),
      position = "dodge"
    ) +
    geom_pointrange(
      aes(color = tstage, group = tstage),
      position = position_dodge(width = 0.9),
      linewidth = 0.75,
      fill = "whitesmoke",
      shape = 21
    ) +
    facet_grid(
      rows = vars(sex),
      labeller = labeller(
        tstage = \(x) stringr::str_replace(get_tstage_label(x), "M", "Unspecified"),
        end = \(x) paste(x, "year net survival")
      )
    ) +
    scale_y_continuous(
      breaks = scales::breaks_width(0.2),
      expand = expansion(add = c(0, 0.1))
    ) +
    scale_x_discrete(
      expand = expansion(0.1, 0.1),
      labels = \(x) gsub(" ", "\n", x)
    ) +
    scale_color_manual(
      values = col_pal_d,
      labels = \(x) stringr::str_replace(get_tstage_label(x), "M", "Unspecified"),
    ) +
    scale_fill_manual(
      values = alpha(col_pal_l, 0.8),
      labels = \(x) stringr::str_replace(get_tstage_label(x), "M", "Unspecified"),
    ) +
    expand_limits(y = c(0, 1)) +
    theme(
      legend.position = if (layout == "vertical") "right" else "bottom",
      legend.justification = if (layout == "vertical") "bottom" else "left",
      panel.spacing.x = unit(10, "pt"),
      panel.spacing.y = unit(10, "pt"),
      axis.text.x = element_text(size = unit(16, "pt")),
      axis.text.y = element_text(size = unit(14, "pt")),
      axis.title.y = element_text(size = unit(18, "pt")),
      strip.text = element_text(size = unit(18, "pt"))
    ) +
    labs(
      x = NULL,
      y = "Net survival (%)",
      color = "T category",
      fill = "T category"
    )
  
  
  # if (!is.null(fup)) {
  #   plot_out <- plot_out +
  #     annotate(
  #       geom = "rect", xmin = -Inf, xmax = Inf, ymin = -0.3, ymax = 0,
  #       fill = "whitesmoke"
  #     ) +
  #     geom_hline(yintercept = 0, color = "#0f0f0f") +
  #     geom_text(
  #       data = count,
  #       aes(
  #         x = clinical_stage, y = -0.05, label = N,
  #         color = tstage, group = tstage
  #       ),
  #       size = rel(3.5), show.legend = FALSE,
  #       position = position_dodge(width = 0.9),
  #       angle = 90, hjust = 1,
  #       inherit.aes = FALSE
  #     ) +
  #     coord_cartesian(ylim = c(-0.3, min(max(dta[["cns_pp"]]), 1.5)))
  # } else {
  #   plot_out <- plot_out +
  #     coord_cartesian(ylim = c(0, min(max(dta[["cns_pp"]]), 1.5)))
  # }
  plot_out <- plot_out +
    coord_cartesian(ylim = c(0, min(max(dta[["cns_pp"]]), 1.5)))
  
  return(plot_out)
}

plt <- plot_by_cstage(plot_data, count_data, fup = 5, layout = "vertical")
plt_ann <- plt + annotate(
    geom = "rect", xmin = 2.5, ymin = 0, xmax = Inf, ymax = Inf,
    fill = "maroon", alpha = 0.05, color = "maroon"
  )
ggsave(plt, filename = "plots/Table-1a-Plot-5v.svg", width = 7, height = 6, scale = 1.4)
ggsave(plt_ann, filename = "plots/Table-1a-Plot-5v-ann.svg", width = 7, height = 6, scale = 1.4)


## ---- Plots: BarNS: Local: By Sex and Tstage ----
plot_data <- readRDS(here::here("Data/bar-local.Rds"))
count_data <- counts[["by-sex-tstage-year_cat-clinical_stage"]] %>% 
  copy() %>%
  .[sex != "Overall" & tstage != "Overall" & year_cat != "1983-2019"] %>% 
  .[, sex := fifelse(sex == "female", "Women", "Men")] %>% 
  .[, N := if_else(N <= 5, "≤5", as.character(N))] %>% 
  .[grepl("Local", clinical_stage)]

plot_by_tstage <- function(dta, count, fup = NULL) {
  if (!is.null(fup)) {
    plot_data <- plot_data %>% dplyr::filter(end == fup)
  }  

  plot_out <- plot_data %>% 
    ggplot(aes(year_cat, cns_pp, ymin = lo_cns_pp, ymax = hi_cns_pp)) +
    geom_col(aes(fill = sex, color = sex), position = "dodge") +
    geom_pointrange(
      aes(color = sex),
      position = position_dodge(width = 0.9),
      linewidth = 0.75,
      fill = "whitesmoke",
      shape = 21
    ) +
    ggh4x::facet_grid2(
      cols = vars(tstage),
      rows = if (is.null(fup)) vars(end) else NULL,
      strip = ggh4x::strip_nested(),
      labeller = labeller(
        tstage = \(x) stringr::str_replace(get_tstage_label(x), "M", "Unspecified"),
        end = \(x) paste(x, "year")
      )
    ) +
    scale_y_continuous(
      breaks = scales::breaks_width(0.2),
      expand = expansion(add = c(0, 0.1))
    ) +
    scale_x_discrete(
      expand = expansion(0.1, 0.1)
    ) +
    scale_color_manual(values = c("steelblue3", "tomato3")) +
    scale_fill_manual(values = alpha(c("steelblue3", "tomato3"), 0.8)) +
    expand_limits(y = c(0, 1)) +
    theme(
      legend.position = "bottom",
      legend.justification = "left",
      panel.spacing.x = unit(10, "pt"),
      panel.spacing.y = unit(10, "pt")
    ) +
    labs(
      x = "Period of diagnosis",
      y = "Net survival (%)",
      color = "Sex",
      fill = "Sex",
      caption = glue::glue("Localized Melanoma, Follow up: {fup} years")
    )
  
  if (!is.null(fup)) {
    plot_out <- plot_out +
      annotate(
        geom = "rect", xmin = -Inf, xmax = Inf, ymin = -0.3, ymax = 0,
        fill = "whitesmoke"
      ) +
      geom_hline(yintercept = 0, color = "#0f0f0f") +
      geom_text(
        data = count,
        aes(
          x = year_cat, y = -0.05, label = N,
          color = sex, group = sex
        ),
        size = rel(3), show.legend = FALSE,
        position = position_dodge(width = 0.9),
        angle = 90, hjust = 1,
        inherit.aes = FALSE
      ) +
      coord_cartesian(ylim = c(-0.3, min(max(dta[["cns_pp"]]), 1.5)))
  } else {
    plot_out <- plot_out +
      coord_cartesian(ylim = c(0, min(max(dta[["cns_pp"]]), 1.5)))
  }
  return(plot_out)
}

plt <- plot_by_tstage(dta = plot_data, count = count_data, fup = 5)
ggsave(plt, filename = "plots/Table-1-Plot-5.svg", width = 10, height = 4, scale = 1.2)

## ---- Plot: CumInc: By tstage per cstage ----
data <- readRDS(here::here("Data/cuminc.Rds"))
rt <- readRDS(here::here("Data/cuminc-rt.Rds"))
get_plt <- function(data, rate_table, pattern = "Local") {
  dta <- data[!is.na(tstage)][grepl(pattern, Stage)]
  rt <- na.omit(rate_table)[grepl(pattern, Stage)]
  caption <- fcase(
    grepl("[lL]ocal", pattern), "Localized cases",
    grepl("[aA]ll", pattern), "All cases",
    grepl("[rR]egion", pattern), "Regional metastasis",
    grepl("[dD]istant", pattern), "Distant metastasis",
    default = NULL
  )
  rt[, tstage := forcats::fct_relabel(
    tstage,
    stringr::str_replace,
    "M", "Unspecified"
  )]
  dta[, tstage := forcats::fct_relabel(
    tstage,
    stringr::str_replace,
    "M", "Unspecified"
  )]

  rt[, hjust := fifelse(time %in% range(time), "inward", "0.5")]
  rt <- rt[dta[, max(time), by = .(Stage, sex, tstage, outcome)], 
    on = .(Stage, sex, tstage, outcome)][time <= V1 + 2 * 365.241][ , V1 := NULL]
  
  plt_tstage_sex <- ggplot(
    data = dta,
    aes(
      x = time,
      y = estimate,
      fill = outcome,
      color = outcome,
      group = interaction(sex, outcome, Stage),
      linetype = Stage
    )
  ) +
    geom_ribbon(
      aes(ymin = conf.low, ymax = conf.high),
      color = NA, alpha = 0.35,
    ) +
    geom_step(
      na.rm = TRUE
    ) +
    facet_grid(
      rows = vars(sex),
      cols = vars(tstage),
      labeller = ggplot2::labeller(
        sex = snakecase::to_title_case,
        tstage = get_tstage_label,
        year_cat = \(x) gsub("--", "–", x)
      )
    ) +
    plot_commons() +
    scale_linetype_discrete(guide = "none") +
    labs(
      caption = caption
    ) +
    guides(
      color = guide_legend(title.position = "top"),
      fill = guide_legend(title.position = "top")
    )

  plt_tstage_sex +
    annotate(
      xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0,
      geom = "rect", fill = "whitesmoke"
    ) +
    geom_hline(yintercept = 0) +
    expand_limits(y = c(-0.4, 1)) +
    coord_cartesian(y = c(-0.4, 1)) +
    geom_text(
      data = rt,
      aes(label = value, y = y, x = time, hjust = hjust),
      color = "grey20",
      size = rel(3.2),
      show.legend = FALSE,
      na.rm = TRUE,
      vjust = "inward"
    ) +
    scale_y_continuous(
      breaks = c(seq(-0.4, 0, 0.1), seq(0, 1, 0.2)),
      expand = ggplot2::expansion(add = c(0.05, 0)),
      labels = function(x) {
        x <- scales::percent(x, accuracy = 1, suffix = "")
        x[4:1] <- c("At risk", "Melanoma", "Other", "Censored")
        return(x)
      }
    )
}

plt_local <- get_plt(data, rt, "Local")
ggsave(plt_local, filename = "plots/cuminc-local.svg", width = 10, height = 5, scale = 1.4)

plt_regional <- get_plt(data, rt, "Regional")
ggsave(plt_regional, filename = "plots/cuminc-regional.svg", width = 10, height = 5, scale = 1.4)

plt_distant <- get_plt(data, rt, "Distant")
ggsave(plt_distant, filename = "plots/cuminc-distant.svg", width = 10, height = 5, scale = 1.4)

## ---- Plot: CumInc: Treatment ----
data <- readRDS(here::here("Data/cuminc-trt.Rds"))
rt <- readRDS(here::here("Data/cuminc-trt-rt.Rds"))

get_trt_plt <- function(data, rt) {
  ggplot(
    data = data,
    aes(x = time, y = estimate, color = outcome)
  ) +
    ggsurvfit::stat_stepribbon(
      aes(ymin = conf.low, ymax = conf.high, fill = outcome),
      alpha = 0.25, color = NA,
      na.rm = TRUE
    ) +
    geom_step(na.rm = TRUE, linewidth = 0.75) +
    ggh4x::facet_grid2(
      rows = vars(sex),
      cols = vars(period),
      scales = "free_x",
      independent = "x",
      labeller = ggplot2::labeller(
        sex = snakecase::to_title_case,
        period = \(x) paste0("Diag. period: ", gsub("-", "–", x))
      ),
      strip = ggh4x::strip_nested()
    ) +
    plot_commons() +
    scale_x_continuous(
      breaks = scales::breaks_width(365.241),
      labels = \(x) x / 365.241,
    ) +
    annotate(
      xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0,
      geom = "rect", fill = "whitesmoke"
    ) +
    geom_hline(yintercept = 0) +
    geom_text(
      data = rt,
      aes(label = value, y = y, x = time),
      color = "grey20",
      size = rel(4),
      family = "mono",
      show.legend = FALSE,
      na.rm = TRUE,
      vjust = "inward"
    ) +
    expand_limits(y = c(-0.3, 1)) +
    scale_y_continuous(
      breaks = c(seq(-0.3, 0, 0.1), seq(0, 1, 0.2)),
      expand = ggplot2::expansion(add = c(0.05, 0)),
      labels = function(x) {
        x <- scales::percent(x, accuracy = 1, suffix = "")
        x[3:1] <- c("At risk", "Melanoma", "Other")
        return(x)
      }
    )
}

plt <- get_trt_plt(data, rt)
ggsave(plt, filename = "plots/cuminc-trt.svg", width = 7, height = 6, scale = 1.2)

## ---- Plot: CondSurv: CRN ----
dta <- readRDS("Data/crn-cond-surv-5.Rds")
plot_df <- dta[[1]]
axis_df <- dta[[2]]
get_cond_plot <- function(plot_df, axis_df, show_cond = TRUE) {
  plot_out <- plot_df %>%
    ggplot(aes(end, cns_pp, color = tstage)) +
    geom_step(
      data = ~ subset(.x, cond == 0 & end <= 15 & type == "surv"),
      aes(linetype = "Net survival")
    )
  if (!show_cond) {
    plot_out <- plot_out +
      ggsurvfit::stat_stepribbon(
        data = ~ subset(.x, cond == 0 & end <= 15 & type == "surv"),
        aes(ymin = lo_cns_pp, ymax = hi_cns_pp, fill= tstage),
        color = NA, alpha = 0.15
      )
  }
  if (show_cond) {
    plot_out <- plot_out + 
      geom_line(
        data = ~ subset(.x, type == "cond" & end <= 20),
        aes(group = interaction(id, tstage), linetype = "Conditional net survival")
      )
  }
  plot_out <- plot_out +
    facet_grid(
      cols = vars(sex),
      labeller = labeller(
        sex = snakecase::to_title_case,
      )
    ) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.justification = "left",
      legend.box.just = "left",
      plot.caption = element_text(hjust = 0)
    ) +
    scale_fill_manual(
      name = "T category",
      values = col_pal_d,
      labels = get_tstage_label,
      guide = guide_legend(title.position = "top", nrow = 2, order = 1)
    ) +
    scale_color_manual(
      name = "T category",
      values = col_pal_d,
      labels = get_tstage_label,
      guide = guide_legend(title.position = "top", nrow = 2, order = 1)
    ) +
    scale_linetype_manual(
      values = c(1, 2),
      breaks = c("Net survival", "Conditional net survival"),
      guide = if(show_cond) {
        guide_legend(nrow = 2, order = 2, title.position = "top") 
      } else {
        "none"
      }
    ) +
    scale_y_continuous(
      breaks = scales::breaks_width(0.1),
      minor_breaks = scales::breaks_width(0.05),
      labels = scales::label_percent(suffix = ""),
      limits = c(0, 1.05)
    ) +
    scale_x_continuous(
      breaks = scales::breaks_width(2)
    ) +
    labs(
      x = "Follow-up time, years",
      y = "Net survival (%)",
      linetype = "5-year",
    ) +
    expand_limits(y = 0)
  
  if (show_cond) {
  plot_out <- plot_out +
    geom_segment(
      data = axis_df[id == 5],
      aes(y = 0, yend = 0.02, x = end, xend = end),
      color = "grey20", linewidth = 0.25
    ) +
    geom_segment(
      data = unique(axis_df[id == 5, -c(7, 8), with = FALSE]),
      aes(y = 0, yend = 0, x = xmin, xend = xmax),
      color = "grey20", linewidth = 0.25
    ) +
    geom_text(
      data = axis_df[id == 5 & (end %% 2 == 0 & id == 10) | (end %% 2 == 1 & id == 5)],
      aes(y = 0.06, x = end, label = label),
      size = rel(4),
      color = "grey10"
    ) +
    annotate(
      geom = "text",
      label = "Conditional year of survival",
      x = 10, y = 0,
      hjust = 0.5, vjust = 1.5,
      color = "grey10",
      size = rel(3)
    )
  }
  return(plot_out)
}

plt_ns <- get_cond_plot(plot_df, axis_df, show_cond = FALSE)
plt_cns <- get_cond_plot(plot_df, axis_df, show_cond = TRUE)

ggsave(plt_ns, filename = "plots/net-surv-5.svg", width = 7.8, height = 6.5, scale = 1)
ggsave(plt_cns, filename = "plots/condsurv-crn5.svg", width = 7.8, height = 6.5, scale = 1)



