# required libraries ----
library(shiny)            # creates shiny app
library(tidyverse)        # group of tools for data manipulation and creating output graphics

# server general section ----

## colorblind friendly palette ----
cb_palette <- c("#000000", "#CC79A7", "#D55E00", "#56B4E9", "#E69F00", "#009E73", "#999999", "#0072B2")

## plot theme ----
mytheme <-  theme(
  panel.background = element_rect(fill = "white", colour = NA),
  panel.grid.major = element_line(colour = "grey70", linewidth =  0.2),
  panel.grid.minor = element_line(colour = "grey85", linewidth =  0.5),
  legend.background = element_blank(),
  legend.key = element_blank(),
  legend.key.width = unit(4, "line"),
  legend.title = element_text(face = "bold", size = 18),
  legend.text = element_text(face = "bold", size = 16),
  legend.position = "right",
  legend.direction = "vertical",
  strip.background = element_blank(),
  strip.text = element_text(face = "bold", size = rel (1.5)),
  axis.ticks = element_line(colour = "black", linewidth =  1),
  axis.line = element_line(colour = "black", linewidth =  1, lineend = "square"),
  axis.text.x = element_text(colour = "black", size = 16),
  axis.text.y = element_text(colour = "black", size = 16),
  axis.title.x = element_text(size = 18),
  axis.title.y = element_text(size = 18),
  plot.title = element_text(face = "bold", size = 20, hjust = 0.5))

## molar change in enthalpy of reaction delta H (kJ/mol) ----

### free bromine acid-base ----
# Benjamin (2015) Table A.1, Page 848, under Bromide (Br)
#   Hypobromous acid, HOBr = -113.0 kJ/mol
#   Hypobromite ion, OBr- = -94.1 kJ/mol
#   Hydrogen ion, H+ = 0 kJ/mol
#   HOBr <--> OBr- + H+
#   Sum of Products - Sum of Reactants = (-94.1 + 0) - (-113.0) = 18.9
delta_H_br <- 18.9

### orthophosphate acid-base ----
# Goldberg et al. (2002), Page 321
delta_H_p_1 <- -8.0
delta_H_p_2 <- 3.6
delta_H_p_3 <- 16.0

## molar change in heat capacity of reaction delta Cp (kJ/K/mol) ----

### orthophosphate acid-base ----
# Goldberg et al. (2002), Page 321
delta_Cp_p_1 <- -0.141
delta_Cp_p_2 <- -0.230
delta_Cp_p_3 <- -0.242

## constants for equilibrium constant equation types 1 or 2 ----
# Type 1 is a 5 constant equation
# Type 2 is a 3 constant equation (parameters 4 and 5 are set to zero)

### carbonate acid-base ----
# Plummer & Busenberg (1982)
# Type 1
A_c_1 <- c(-356.3094, -0.06091964, 21834.37, 126.8339, -1684915)
A_c_2 <- c(-107.8871, -0.03252849, 5151.79, 38.92561, -563713.9)

### free ammonia acid-base ----
# Bates & Pinching (1949)
# Type 2
A_n <- c(0.6322, -0.001225, -2835.76, 0, 0)

### free chlorine acid-base ----
# Morris (1966)
# Type 2
A_cl <- c(10.0686, -0.0253, -3000, 0, 0)

### orthosilicate acid-base ----
# Adjusted from Busey & Mesmer (1977) to as shown in Nordstrom et al. (1990)
# Type 1
A_s_1 <- c(-302.3724, -0.05069842, 15669.69, 108.18466, -1119669)
# Volosov et al. (1972)
# Type 2
A_s_2 <- c(8.354, -0.021962, -4465.2, 0, 0)

### water dissociation ----
# Nordstrom et al. (1990)
# Type 1
A_w <- c(-283.971, -0.05069842, 13323, 102.24447, -1119669)

## molecular weights ----
CaCO3_equiv_MW <- 50.0
carbon_MW <- 12.011
chlorine_MW <- 70.906
nitrogen_MW <- 14.0067
orthophosphate_MW <- 94.97136
silica_MW <- 60.0848

## 25 Celsius value in Kelvin ----
T_K_25C <- 25 + 273.15

## universal gas constant R (kJ/K/mol) ----
R <- 0.00831451

## function to calculate equilibrium constants by equation types 1 or 2 ----
equilibrium_1_2 <- function(A, T_K) {
  10^(A[1] + A[2] * T_K + A[3] / T_K + A[4] * log10(T_K) + A[5] / T_K^2)
}

## function to calculate equilibrium constants by equation types 3 or 4 ----
equilibrium_3_4 <- function(A, T_K) {
  10^(A[1] / (2.303 * R) * (1 / T_K_25C - 1 / T_K) +
        A[2] / (2.303 * R) * (T_K_25C / T_K - 1 - 2.303 * log10(T_K_25C / T_K)) + log10(A[3]))
}

## equilibrium constants at 25 C & 0 M ionic strength ----

### carbonate acid-base ----
K_c_1_25_0 <- equilibrium_1_2(A_c_1, T_K_25C)
K_c_2_25_0 <- equilibrium_1_2(A_c_2, T_K_25C)

### free ammonia acid-base ----
K_n_25_0 <- equilibrium_1_2(A_n, T_K_25C)

### free bromine acid-base ----
# Benjamin (2015) Table A.2, Page 856
K_br_25_0 <- 10^-8.63

### free chlorine acid-base ----
K_cl_25_0 <- equilibrium_1_2(A_cl, T_K_25C)

### orthophosphate acid-base ----
# Goldberg et al. (2002), Page 321
K_p_1_25_0 <- 10^-2.148
K_p_2_25_0 <- 10^-7.198
K_p_3_25_0 <- 10^-12.35

### orthosilicate acid-base ----
K_s_1_25_0 <- equilibrium_1_2(A_s_1, T_K_25C)
K_s_2_25_0 <- equilibrium_1_2(A_s_2, T_K_25C)

### water dissociation ----
K_w_25_0 <- equilibrium_1_2(A_w, T_K_25C)

## constants for equilibrium constant equation types 3 or 4 ----
# Type 3 includes, change in enthalpy, change in heat capacity, and 25C equilibrium constant
# Type 4 includes only change in enthalpy and 25C equilibrium constant (change in heat capacity is set to zero)

### free bromine acid-base ----
# Benjamin (2015)
# Type 3
A_br <- c(delta_H_br, 0, K_br_25_0)

### orthophosphate acid-base ----
# Goldberg et al. (2002), Page 321
# Type 4
A_p_1 <- c(delta_H_p_1, delta_Cp_p_1, K_p_1_25_0)
A_p_2 <- c(delta_H_p_2, delta_Cp_p_2, K_p_2_25_0)
A_p_3 <- c(delta_H_p_3, delta_Cp_p_3, K_p_3_25_0)


## function to calculate activity coefficient (gamma) ----
davies_eq <- function(charge, is_M_gamma, T_C_gamma, T_K_gamma) {

  # calculate dielectric constant for given temperature ----
  # Malmberg & Maryott (1956)
  epsilon <- 87.74 - 0.4008 * T_C_gamma + T_C_gamma^2 * 0.0009398 - T_C_gamma^3 * 0.00000141

  # calculate constant A for Davies equation ----
  # Benjamin (2015)
  A <- 1824830 * (epsilon * T_K_gamma)^-1.5

  # calculate activity coefficients per Davies equation ----
  # Benjamin (2015)
  gamma <- 10^(-A * charge^2 * (sqrt(is_M_gamma) / (1 + sqrt(is_M_gamma)) - 0.30 * is_M_gamma))

  # return desired activity coefficient ----
  return (gamma)

}

## function to calculate alpha values and buffer intensity term ----
alpha_buffer_term <- function(acid_type, alpha_or_buffer, H, K_1, K_2, K_3) {

  # monoprotic acids ----
  if (acid_type == "mono") {

    ## calculate alpha values ----
    alpha_0 <- 1 / (1       + K_1 / H)
    alpha_1 <- 1 / (H / K_1 + 1)

    ## return alpha values ----
    if (alpha_or_buffer == "alpha") {
      return(list(alpha_0, alpha_1))
    }

    ## calculate buffer term ----
    buffer_term <- 2.303 * alpha_0 * alpha_1

    ## return buffer term ----
    if (alpha_or_buffer == "buffer") {
      return(buffer_term)
    }
  }

  # diprotic acids ----
  if (acid_type == "di") {

    ## calculate alpha values ----
    alpha_0 <- 1 / (1                 + K_1 / H + K_1 * K_2 / H^2)
    alpha_1 <- 1 / (H / K_1           + 1       + K_2 / H)
    alpha_2 <- 1 / (H^2 / (K_1 * K_2) + H / K_2 + 1)

    ## return alpha values ----
    if (alpha_or_buffer == "alpha") {
      return(list(alpha_0, alpha_1, alpha_2))
    }

    ## calculate buffer term ----
    buffer_term <- 2.303 * (alpha_0 * alpha_1 +
                              4 * alpha_0 * alpha_2 +
                              alpha_1 * alpha_2)

    ## return buffer term ----
    if (alpha_or_buffer == "buffer") {
      return(buffer_term)
    }
  }

  # triprotic acids ----
  if (acid_type == "tri") {

    ## calculate alpha values ----
    alpha_0 <- 1 / (1                       + K_1 / H           + K_1 * K_2 / H^2 + K_1 * K_2 * K_3 / H^3)
    alpha_1 <- 1 / (H / K_1                 + 1                 + K_2 / H         + K_2 * K_3 / H^2)
    alpha_2 <- 1 / (H^2 / (K_1 * K_2)       + H / K_2           + 1               + K_3 / H)
    alpha_3 <- 1 / (H^3 / (K_1 * K_2 * K_3) + H^2 / (K_2 * K_3) + H / K_3         + 1)

    ## return alpha values ----
    if (alpha_or_buffer == "alpha") {
      return(list(alpha_0, alpha_1, alpha_2, alpha_3))
    }

    ## calculate buffer term ----
    buffer_term <- 2.303 * (alpha_0 * alpha_1 +
                              4 * alpha_0 * alpha_2 +
                              9 * alpha_0 * alpha_3 +
                              alpha_1 * alpha_2 +
                              4 * alpha_1 * alpha_3 +
                              alpha_2 * alpha_3)

    ## return buffer term ----
    if (alpha_or_buffer == "buffer") {
      return(buffer_term)
    }
  }
}

## function to calculate buffer intensity ----
buffer_solver <- function(pH_range, pH_delta, temp, is_input, is_mM, TDS_mg_L, EC_micros_cm, DIC_input,
                          temp_alk, pH_alk, alk_mg_L, DIC_mg_C_L, TOTPO4_mg_PO4_L,
                          TOTNH3_mg_N_L, TOTBr_mg_Cl_L, TOTCl_mg_Cl_L, TOTSiO2_mg_SiO2_L) {

  # set pH range ----
  pH <- seq(pH_range[1], pH_range[2], by = 0.01)

  # calculate temperature in Kelvin ----
  T_K <- temp + 273.15

  # calculate ionic strength in molar concentration ----

  ## known ionic strength input ----
  if (is_input == "known_is") {
    is_M <- is_mM / 1000
  }

  ## estimate ionic strength from total dissolved solids ----
  # Tchobanoglous et al. (2003), page 40
  if (is_input == "est_is_TDS") {
    is_M <- 2.5e-5 * TDS_mg_L
  }

  ## estimate ionic strength from electrical conductivity ----
  # Tchobanoglous et al. (2003), page 56
  if (is_input == "est_is_EC") {
    is_M <- 1.6e-5 * EC_micros_cm
  }

  # calculate activity corrections ----
  gamma_1 <- davies_eq(1, is_M, temp, T_K)
  gamma_2 <- davies_eq(2, is_M, temp, T_K)
  gamma_3 <- davies_eq(3, is_M, temp, T_K)

  # calculate ionic strength corrected equilibrium constants at specified temperature ----

  ## carbonate acid-base ----
  K_c_1 <- equilibrium_1_2(A_c_1, T_K) / gamma_1
  K_c_2 <- gamma_1 * equilibrium_1_2(A_c_2, T_K) / gamma_2

  ## free ammonia acid-base ----
  K_n <- gamma_1 * equilibrium_1_2(A_n, T_K)

  ## free bromine acid-base ----
  K_br <- equilibrium_3_4(A_br, T_K) / gamma_1

  ## free chlorine acid-base ----
  K_cl <- equilibrium_1_2(A_cl, T_K) / gamma_1

  ## orthophosphate acid-base ----
  K_p_1 <- equilibrium_3_4(A_p_1, T_K) / gamma_1
  K_p_2 <- gamma_1 * equilibrium_3_4(A_p_2, T_K) / gamma_2
  K_p_3 <- gamma_2 * equilibrium_3_4(A_p_3, T_K) / gamma_3

  ## orthosilicate acid-base ----
  K_s_1 <- equilibrium_1_2(A_s_1, T_K) / gamma_1
  K_s_2 <- gamma_1 * equilibrium_1_2(A_s_2, T_K) / gamma_2

  ## water dissociation ----
  K_w <- equilibrium_1_2(A_w, T_K)

  # initial concentrations converted to molar ----

  ## total carbonate (M) ----

  ### initialize DIC estimate flag ----
  est_TOTCO3_flag <- 0

  ### known DIC concentration input ----
  if (DIC_input == "known_DIC") {
    TOTCO3 <- DIC_mg_C_L / carbon_MW / 1000
  }

  ### estimate DIC from total alkalinity water sample ----
  # Benjamin (2015)
  if (DIC_input == "est_alk") {

    #### total alkalinity water sample temperature in Kelvin ----
    T_K_alk <- temp_alk + 273.15

    #### calculate activity corrections for total alkalinity water sample ----

    ##### calculate dielectric constant for given total alkalinity water sample temperature ----
    gamma_1_alk <- davies_eq(1, is_M, temp_alk, T_K_alk)
    gamma_2_alk <- davies_eq(2, is_M, temp_alk, T_K_alk)

    #### total alkalinity water sample water dissociation ----
    K_w_alk <- equilibrium_1_2(A_w, T_K_alk)

    #### total alkalinity water sample carbonate acid-base ----
    K_c_1_alk <- equilibrium_1_2(A_c_1, T_K_alk) / gamma_1_alk
    K_c_2_alk <- gamma_1_alk * equilibrium_1_2(A_c_2, T_K_alk) / gamma_2_alk

    #### total alkalinity water sample water species activities ----
    H_plus_alk <- 10^-pH_alk
    OH_minus_alk <- K_w_alk / H_plus_alk

    #### total alkalinity water sample carbonate alpha values ----
    alpha_c_alk <- alpha_buffer_term("di", "alpha", H_plus_alk, K_c_1_alk, K_c_2_alk, NA)
    alpha_1_c_alk <- alpha_c_alk[[2]]
    alpha_2_c_alk <- alpha_c_alk[[3]]

    #### estimated DIC from total alkalinity water sample ----
    # Benjamin (2015)
    TOTCO3 <- (alk_mg_L / CaCO3_equiv_MW / 1000 + gamma_1_alk * H_plus_alk - gamma_1_alk * OH_minus_alk) / (alpha_1_c_alk + 2 * alpha_2_c_alk)

    #### set DIC estimate to zero if DIC estimate is negative and set flag noting impossible ----
    if(TOTCO3 < 0) {
      TOTCO3 <- 0
      est_TOTCO3_flag <- -1
    }
  }

  ## total orthophosphate (M) ----
  TOTPO4 <- TOTPO4_mg_PO4_L / orthophosphate_MW / 1000

  ## total ammonia (M) ----
  TOTNH3 <- TOTNH3_mg_N_L / nitrogen_MW / 1000

  ## total bromine (M) ----
  TOTBr <- TOTBr_mg_Cl_L / chlorine_MW / 1000

  ## total chlorine (M) ----
  TOTCl <- TOTCl_mg_Cl_L / chlorine_MW / 1000

  ## total orthosilicate (M) ----
  TOTSiO2 <- TOTSiO2_mg_SiO2_L / silica_MW / 1000

  # total carbonate (mg carbon/L) ----
  TOTCO3_mg_C_L <- TOTCO3 * carbon_MW * 1000

  # water species activities ----
  H_plus <- 10^-pH
  OH_minus <- K_w / H_plus

  # buffer intensity calculations ----

  ## water buffer intensity ----
  buffer_water <- 2.303 * gamma_1 * (H_plus + OH_minus)

  ## carbonate buffer intensity ----
  buffer_carbonate <- TOTCO3 * alpha_buffer_term("di", "buffer", H_plus, K_c_1, K_c_2, NA)

  ## orthophosphate buffer intensity ----
  buffer_orthophosphate <- TOTPO4 * alpha_buffer_term("tri", "buffer", H_plus, K_p_1, K_p_2, K_p_3)

  ## free ammonia buffer intensity ----
  buffer_ammonia <- TOTNH3 * alpha_buffer_term("mono", "buffer", H_plus, K_n, NA, NA)

  ## free bromine buffer intensity ----
  buffer_bromine <- TOTBr * alpha_buffer_term("mono", "buffer", H_plus, K_br, NA, NA)

  ## free chlorine buffer intensity ----
  buffer_chlorine <- TOTCl * alpha_buffer_term("mono", "buffer", H_plus, K_cl, NA, NA)

  ## orthosilicate buffer intensity ----
  buffer_orthosilicate <- TOTSiO2 * alpha_buffer_term("di", "buffer", H_plus, K_s_1, K_s_2, NA)

  ## total buffer intensity ----
  buffer_total <- buffer_water + buffer_carbonate + buffer_orthophosphate +
    buffer_ammonia + buffer_bromine + buffer_chlorine + buffer_orthosilicate

  # calculate required strong acid or strong base added for desired pH change ----
  if (!(is.na(pH_delta[1]))) {

    acid_base_added <- data.frame(pH = pH,
                                  total_meq = buffer_total * 1000) %>%

      # filter pH between desired pH change
      filter(between(pH, pH_delta[1], pH_delta[2])) %>%

      # create column of interpolated areas
      mutate(delta_acid_base_meq = (lead(total_meq) + total_meq) / 2 * (lead(pH) - pH))

    # create output line for strong acid or strong base required for desired pH change ----
    pH_delta_acid_base <- paste0("Strong Acid (pH ",
                                 format(round(pH_delta[2], 2), nsmall = 2),
                                 " to ",
                                 format(round(pH_delta[1], 2), nsmall = 2),
                                 ") or Strong Base (pH ",
                                 format(round(pH_delta[1], 2), nsmall = 2),
                                 " to ",
                                 format(round(pH_delta[2], 2), nsmall = 2),
                                 ") Addition for pH Change (milliequivalent/L) = ",
                                 format(round(sum(acid_base_added$delta_acid_base_meq, na.rm = TRUE), 3), nsmall = 3))

    # return data frame from function ----
    return(pH_delta_acid_base)
  }

  if (is.na(pH_delta)) {
    # create data frame for buffer intensity plot ----
    buffer_plot_data <- data.frame(pH = pH,
                                   Total = buffer_total,
                                   Ammonia = buffer_ammonia,
                                   Bromine = buffer_bromine,
                                   Carbonate = buffer_carbonate,
                                   Chlorine = buffer_chlorine,
                                   Orthophosphate = buffer_orthophosphate,
                                   Orthosilicate = buffer_orthosilicate,
                                   Water = buffer_water) %>%

      # make data frame long
      pivot_longer(cols = Total:Water, names_to = "buffer", values_to = "buffer_intensity_M") %>%

      # get buffer intensity in mM
      mutate(buffer_intensity_mM = buffer_intensity_M * 1000,

             # make buffer a factor
             buffer = factor(buffer, levels = c("Total", "Ammonia", "Bromine", "Carbonate",
                                                "Chlorine", "Orthophosphate", "Orthosilicate", "Water")))

    # combine simulated condition names and values in a table ----
    Parameter <- c("pH range lower value",
                   "pH range upper value",
                   "Temperature (Celsius)",
                   "Ionic Strength (mM)",
                   "Dissolved Inorganic Carbon (mg carbon/L)",
                   "Total Dissolved Orthophosphate (mg orthophosphate/L)",
                   "Free Ammonia (mg nitrogen/L)",
                   "Free Bromine (mg chlorine/L)",
                   "Free Chlorine (mg chlorine/L)",
                   "Total Dissolved Orthosilicate (mg silica/L)")

    Value <- c(format(round(pH_range[1], 2), nsmall = 2),
               format(round(pH_range[2], 2), nsmall = 2),
               format(round(temp, 1), nsmall = 1),
               format(round(is_M * 1000, 1), nsmall = 1),
               format(round(TOTCO3_mg_C_L, 1), nsmall = 1),
               format(round(TOTPO4_mg_PO4_L, 2), nsmall = 2),
               format(round(TOTNH3_mg_N_L, 2), nsmall = 2),
               format(round(TOTBr_mg_Cl_L, 2), nsmall = 2),
               format(round(TOTCl_mg_Cl_L, 2), nsmall = 2),
               format(round(TOTSiO2_mg_SiO2_L, 2), nsmall = 2))

    ic_table <- cbind(Parameter, Value)

    # get equilibrium constant values in a table ----
    Dissociation_Equilibrium <- c("Water",
                                  "First Carbonate",
                                  "Second Carbonate",
                                  "First Orthophosphate",
                                  "Second Orthophosphate",
                                  "Third Orthophosphate",
                                  "Free Ammonia",
                                  "Free Bromine",
                                  "Free Chlorine",
                                  "First Orthosilicate",
                                  "Second Orthosilicate")

    Simulated_pK <- c(format(round(-log10(c(K_w, K_c_1, K_c_2, K_p_1, K_p_2)), 3), nsmall = 3),
                      format(round(-log10(K_p_3), 2), nsmall = 2),
                      format(round(-log10(K_n), 3), nsmall = 3),
                      format(round(-log10(c(K_br, K_cl, K_s_1, K_s_2)), 2), nsmall = 2))

    Standard_pK <- c(format(round(-log10(c(K_w_25_0, K_c_1_25_0, K_c_2_25_0, K_p_1_25_0, K_p_2_25_0)), 3), nsmall = 3),
                     format(round(-log10(K_p_3_25_0), 2), nsmall = 2),
                     format(round(-log10(K_n_25_0), 3), nsmall = 3),
                     format(round(-log10(c(K_br_25_0, K_cl_25_0, K_s_1_25_0, K_s_2_25_0)), 2), nsmall = 2))

    pk_table <- cbind(Dissociation_Equilibrium, Standard_pK, Simulated_pK)

    # create output lines for dissolved inorganic carbon and ionic strength ----

    ## create valid DIC estimate output line ----
    if (est_TOTCO3_flag == 0) {
      est_DIC_mg_L <- paste0("Estimated Dissolved Inorganic Carbon (mg carbon/L) = ", format(round(TOTCO3_mg_C_L, 1), nsmall = 1))
    }

    ## create impossible DIC estimate output line ----
    if (est_TOTCO3_flag == -1) {
      est_DIC_mg_L <- paste0("Chemically Impossible Condition - Estimated Dissolved Inorganic Carbon (mg carbon/L) Set to ", format(round(TOTCO3_mg_C_L, 1), nsmall = 1))
    }

    ## create ionic strength output line ----
    est_IS_mM <- paste0("Estimated Ionic Strength (mM) = ", format(round(is_M * 1000, 1), nsmall = 1))

    # create output lines for simulated pK values ----
    pK_sim <- paste0("pK for Standard (25 \u00B0C & 0 M) & Simulated (", format(round(temp, 1), nsmall = 1),
                     " \u00B0C & ", format(round(is_M * 1000, 1), nsmall = 1), " mM) Conditions")

    # create buffer intensity plot ----

    ## create ggplot ----
    plot_gg <- ggplot(data = buffer_plot_data,
                      aes(x = pH, y = buffer_intensity_mM, color = buffer, linetype = buffer)) +
      geom_line(linewidth = 1) +
      scale_color_manual(values = cb_palette) +
      scale_linetype_manual(values = c(1, 4, 4, 3, 4, 3, 3, 2)) +
      xlab("pH") +
      ylab("Buffer Intensity\n(milliequivalent/L)/pH unit") +
      guides(color = guide_legend(ncol = 1)) +
      labs(color = "Component", linetype = "Component") +
      scale_x_continuous(breaks = scales::pretty_breaks()) +
      scale_y_continuous(breaks = scales::pretty_breaks()) +
      ggtitle("Plot of Buffer Intensity versus pH") +
      mytheme

    ## make a plotly plot ----
    plot <- plotly::ggplotly(plot_gg) %>%

      ## modify hover text ----
    plotly::style(hovertemplate = paste0("pH = %{x:.2f}<br>",
                                         "Buffer Intensity = %{y:.4f} (mequiv/L)/pH unit")) %>%
      ## remove plotly logo and make mode bar always visible ----
    plotly::config(displaylogo = FALSE, displayModeBar = TRUE)

    # return items from function ----
    return(list(plot_gg, plot, buffer_plot_data, ic_table, pk_table, est_DIC_mg_L, est_IS_mM, pK_sim))
  }
}

# server definition section ----

# define server logic required to run simulations and produce output
server <- function(input, output, session) {

  # create input for selecting starting and ending pH values for pH change calculation ----
  output$pH_delta <- renderUI({

    # require that input$pH_range exists before running to set limits on pH change slider
    req(input$pH_range)

    # create the slider input for the pH change starting and ending values
    sliderInput("pH_delta",
                label = p("Select Endpoint pH Values for pH Change Calculation",
                          style = "font-size: 20px"),
                min = input$pH_range[1],
                max = input$pH_range[2],
                value = c(input$pH_range[1], input$pH_range[2]),
                step = 0.05)
  }
  )

  # generate simulation of buffer intensity ----
  sim <- reactive({

    # call function to get outputs
    buffer_solver(pH_range = input$pH_range,
                  pH_delta = NA,
                  temp = input$temp,
                  is_input = input$is_input,
                  is_mM = input$is_mM,
                  TDS_mg_L = input$TDS_mg_L,
                  EC_micros_cm = input$EC_micros_cm,
                  DIC_input = input$DIC_input,
                  temp_alk = input$temp_alk,
                  pH_alk = input$pH_alk,
                  alk_mg_L = input$alk_mg_L,
                  DIC_mg_C_L = input$DIC_mg_C_L,
                  TOTPO4_mg_PO4_L = input$TOTPO4_mg_PO4_L,
                  TOTNH3_mg_N_L = input$TOTNH3_mg_N_L,
                  TOTBr_mg_Cl_L = input$TOTBr_mg_Cl_L,
                  TOTCl_mg_Cl_L = input$TOTCl_mg_Cl_L,
                  TOTSiO2_mg_SiO2_L = input$TOTSiO2_mg_SiO2_L)
    })

  # generate simulation for pH change calculation ----
  sim_pH_delta <- reactive({

    # require that input$pH_delta exists before running
    req(input$pH_delta)

    # call function to get outputs
    buffer_solver(pH_range = input$pH_range,
                  pH_delta = input$pH_delta,
                  temp = input$temp,
                  is_input = input$is_input,
                  is_mM = input$is_mM,
                  TDS_mg_L = input$TDS_mg_L,
                  EC_micros_cm = input$EC_micros_cm,
                  DIC_input = input$DIC_input,
                  temp_alk = input$temp_alk,
                  pH_alk = input$pH_alk,
                  alk_mg_L = input$alk_mg_L,
                  DIC_mg_C_L = input$DIC_mg_C_L,
                  TOTPO4_mg_PO4_L = input$TOTPO4_mg_PO4_L,
                  TOTNH3_mg_N_L = input$TOTNH3_mg_N_L,
                  TOTBr_mg_Cl_L = input$TOTBr_mg_Cl_L,
                  TOTCl_mg_Cl_L = input$TOTCl_mg_Cl_L,
                  TOTSiO2_mg_SiO2_L = input$TOTSiO2_mg_SiO2_L)
  })

  # generate buffer intensity plot as output ----
  output$plot <- plotly::renderPlotly(sim()[[2]])

  # generate simulated condition table ----
  output$ic <- DT::renderDataTable(DT::datatable((sim()[[4]]),
                                                 caption = tags$caption(
                                                   style = "caption-side: bottom; text-align: left; margin: 8px 0; color:black;",
                                                   HTML("")),
                                                 options = list(
                                                   dom = "t",
                                                   pageLength = 15,
                                                   columnDefs = list(
                                                     list(className = "dt-center", targets = "_all"),
                                                     list(orderable = FALSE, targets = c(0,1)),
                                                     list(searchable = FALSE, targets = c(0,1))
                                                   )
                                                 ),
                                                 selection = "none"
                                                 )
                                   )

  # generate pK table ----
  output$pk <- DT::renderDataTable(DT::datatable((sim()[[5]]),
                                                 colnames = c("Dissociation Equilibrium", "Standard pK", HTML("Simulated<sup>1</sup> pK")),
                                                 caption = tags$caption(
                                                   style = "caption-side: bottom; text-align: left; margin: 8px 0; color:black;",
                                                   HTML("<b><sup>1</sup>Chemical equilibrium written with chemical concentrations, except for H<sup>+</sup> and OH<sup>-</sup> which used activities.</b>")),
                                                 options = list(
                                                   dom = "t",
                                                   pageLength = 15,
                                                   columnDefs = list(
                                                     list(className = "dt-center", targets = "_all"),
                                                     list(orderable = FALSE, targets = "_all"),
                                                     list(searchable = FALSE, targets = "_all")
                                                   )
                                                 ),
                                                 escape = FALSE,
                                                 selection = "none")
                                   )

  # generate dissolved inorganic carbon concentration as output ----
  output$sim_DIC_mg_L <- renderText(sim()[[6]])

  # generate ionic strength concentration as output ----
  output$sim_is_mM <- renderText(sim()[[7]])

  # generate ionic strength concentration and temperature as output for pK ----
  output$pK_sim <- renderText(sim()[[8]])

  # generate strong acid or strong base concentration required for desired pH change as output ----
  output$pH_delta_acid_base <- renderText(sim_pH_delta())

  # download buffer intensity plot as output ----
  output$downloadPlot <- downloadHandler(
    filename = function() {
      paste('Buffer_Intensity_Plot_', substr(as.character(Sys.time()), 1, 10),
            '_', substr(as.character(Sys.time()), 12, 13),
            '_', substr(as.character(Sys.time()), 15, 16), '.png', sep = '')
      },
    content = function(file) {
      ggsave(file, plot = sim()[[1]], height = 6, width = 10, device = "png")
      }
    )

  # download buffer intensity data as output ----
  output$downloadData <- downloadHandler(
    filename = function() {
      paste('Buffer_Intensity_Data_', substr(as.character(Sys.time()), 1, 10),
            '_', substr(as.character(Sys.time()), 12, 13),
            '_', substr(as.character(Sys.time()), 15, 16), '.csv', sep = '')
      },
    content = function(file) {write.csv(sim()[[3]], file, row.names = TRUE)
      }
    )

  # download simulated condition table as output ----
  output$downloadIC <- downloadHandler(
    filename = function() {
      paste('Buffer_Intensity_Simulation_Conditions_', substr(as.character(Sys.time()), 1, 10),
            '_', substr(as.character(Sys.time()), 12, 13),
            '_', substr(as.character(Sys.time()), 15, 16), '.csv', sep = '')
    },
    content = function(file) {
      write.csv(sim()[[4]], file, row.names = FALSE)
    }
  )

  # download pK table as output ----
  output$downloadpK <- downloadHandler(
    filename = function() {
      paste('Buffer_Intensity_pK_', substr(as.character(Sys.time()), 1, 10),
            '_', substr(as.character(Sys.time()), 12, 13),
            '_', substr(as.character(Sys.time()), 15, 16), '.csv', sep = '')
    },
    content = function(file) {
      write.csv(sim()[[5]], file, row.names = FALSE)
    }
  )

  }

# ui layout ----
ui <- shinyUI(fluidPage(

## set language attribute for accessibility ----
lang = "en",

## title block ----
h1("Buffer Intensity Simulator (BIS) Source Code, Version 1.0.0, Last Updated October 25, 2024"),

h3("Contact: David G. Wahman (wahman.david@epa.gov), United States Environmental Protection Agency"),

h4(HTML("The code simulates the buffer intensity curve for the selected water chemistry and provides an estimate of the strong acid or strong base
   required to change the pH between two user-selected pH values. Please see the following publication for details on code development
   and illustrations of example uses: Wahman, D. G., Schock, M. R., & Lytle, D. A. (2024). Drinking Water Buffer Intensity Simulator (BIS):
     Development and Practical Simulations. <i>AWWA Water Science</i> 6(6), e70006, "),
        a(target = "_blank", href = "https://doi.org/10.1002/aws2.70006", "https://doi.org/10.1002/aws2.70006",
          .noWS = "outside"), '.', .noWS = c("after-begin", "before-end")),

h5(HTML("<b><u>Disclaimer:</u></b> The United States Environmental Protection Agency (EPA) GitHub project code is provided on an 'as is' basis and the user
    assumes responsibility for its use. EPA has relinquished control of the information and no longer has responsibility to
    protect the integrity, confidentiality, or availability of the information. Any reference to specific commercial products,
    processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement,
    recommendation or favoring by EPA. The EPA seal and logo shall not be used in any manner to imply endorsement of any commercial
    product or activity by EPA or the United States Government. No warranty expressed or implied is made regarding the accuracy
    or utility of the system, nor shall the act of distribution constitute any such warranty. This code has been reviewed in
    accordance with EPA policy and has been approved for external and free use. The views expressed in this code do not necessarily
    represent the views or policies of the EPA. Although a reasonable effort has been made to assure that the results obtained are correct,
    this code is experimental. Therefore, the author and the EPA are not responsible and assume no liability whatsoever
    for any results or any use made of the results obtained from this code, nor for any damages or litigation that result
    from the use of the code for any purpose.")),

## sidebar layout ----
sidebarLayout(

  ### sidebar panel ----
  sidebarPanel(width = 4,

               #### set color of input sliders for accessibility ----
               shinyWidgets::chooseSliderSkin("Shiny", color = "#007cc0"),

               h3(HTML("<b>Select Conditions for Buffer Intensity Simulation</b>")),

               br(),

               #### input for pH range ----
               sliderInput("pH_range",
                           label = p("pH Range",
                                     style = "font-size: 16px"),
                           min = 0.00,
                           max = 14.00,
                           value = c(6.00, 10.00),
                           step = 0.05),

               br(),

               #### input for temperature ----
               sliderInput("temp",
                           label = p("Temperature (Celsius)",
                                     style = "font-size: 16px"),
                           min = 5.0,
                           max = 35.0,
                           value = 25.0,
                           step = 0.1),

               br(),

               #### input for ionic strength ----
               wellPanel(

                 ##### selection of method of input -----
                 radioButtons("is_input",
                              label = p("Select Method for Ionic Strength Input", style = "font-size: 16px"),
                              c("Enter Known Ionic Strength" = "known_is",
                                "Estimate Ionic Strength from Total Dissolved Solids" = "est_is_TDS",
                                "Estimate Ionic Strength from Electrical Conductivity" = "est_is_EC"),
                              selected = "known_is"),

                 ##### known ionic strength input ----
                 conditionalPanel(
                   condition = "input.is_input == 'known_is'",

                   sliderInput("is_mM",
                               label = p("Known Ionic Strength (mM)",
                                         style = "font-size: 14px"),
                               min = 0.0,
                               max = 100.0,
                               value = 5.0,
                               step = 0.5)
                   ),

                 ##### estimate from total dissolved solids ----
                 conditionalPanel(
                   condition = "input.is_input == 'est_is_TDS'",


                   sliderInput("TDS_mg_L",
                               label = p("Measured Water Sample Total Dissolved Solids (mg/L)",
                                         style = "font-size: 14px"),
                               min = 0,
                               max = 1000,
                               value = 200,
                               step = 5)
                   ),

                 ##### estimate from electrical conductivity ----
                 conditionalPanel(
                   condition = "input.is_input == 'est_is_EC'",

                   sliderInput("EC_micros_cm",
                               label = p(HTML("Measured Water Sample Electrical Conductivity (&#181;S/cm or &#181;mho/cm)"),
                                         style = "font-size: 14px"),
                               min = 0,
                               max = 1000,
                               value = 310,
                               step = 5)
                   ),

                 ##### display ionic strength if estimated ----
                 conditionalPanel(
                   condition = "input.is_input == 'est_is_TDS' || input.is_input == 'est_is_EC'",
                   h4(HTML(paste0("<b>", textOutput("sim_is_mM"), "</b>")))
                 )
                 ),

               #### input for dissolved inorganic carbon concentration ----
               wellPanel(

                 ##### selection of method of input -----
                 radioButtons("DIC_input",
                              label = p("Select Method for Dissolved Inorganic Carbon Input", style = "font-size: 16px"),
                              c("Enter Known Dissolved Inorganic Carbon Concentration" = "known_DIC",
                                "Estimate Dissolved Inorganic Carbon from a Total Alkalinity Water Sample" = "est_alk"),
                              selected = "known_DIC"),

                 ##### known DIC concentration input ----
                 conditionalPanel(
                   condition = "input.DIC_input == 'known_DIC'",

                   sliderInput("DIC_mg_C_L",
                               label = p(HTML("Known Dissolved Inorganic Carbon Concentration (mg carbon/L)"),
                                         style = "font-size: 14px"),
                               min = 0.0,
                               max = 100.0,
                               value = 5.0,
                               step = 0.5)
                   ),

                 ##### estimate from total alkalinity water sample inputs ----
                 conditionalPanel(
                   condition = "input.DIC_input == 'est_alk'",

                   ###### input for total alkalinity water sample temperature ----
                   sliderInput("temp_alk",
                               label = p("Total Alkalinity Water Sample Measured Temperature (Celsius)",
                                         style = "font-size: 14px"),
                               min = 5.0,
                               max = 35.0,
                               value = 25.0,
                               step = 0.1),

                   br(),

                   ###### input for total alkalinity water sample pH ----
                   sliderInput("pH_alk",
                               label = p("Total Alkalinity Water Sample Measured Initial pH",
                                         style = "font-size: 14px"),
                               min = 6.00,
                               max = 11.00,
                               value = 8.00,
                               step = 0.05),

                   br(),

                   ###### input for total alkalinity water sample measurement ----
                   sliderInput("alk_mg_L",
                               label = p(HTML("Total Alkalinity Water Sample Measurement (mg/L as CaCO<sub>3</sub>)"),
                                         style = "font-size: 14px"),
                               min = 0,
                               max = 500,
                               value = 100,
                               step = 2),

                   ##### display estimated dissolved inorganic carbon if estimated ----
                   h4(HTML(paste0("<b>", textOutput("sim_DIC_mg_L"), "</b>"))),

                   ##### note that ionic strength used is per above entry ----
                   p(HTML("<b><i>Note: Buffer Itensity Simulation Ionic Strength used in Dissolved Inorganic Carbon Estimation from Total Alkalinity Water Sample</i></b>"),
                     style = "font-size: 12px")
                 )
               ),

               #### input for total dissolved orthophosphate concentration ----
               sliderInput("TOTPO4_mg_PO4_L",
                           label = p("Total Dissolved Orthophosphate Concentration (mg orthophosphate/L)",
                                     style = "font-size: 16px"),
                           min = 0.00,
                           max = 10.00,
                           value = 0.00,
                           step = 0.05),

               br(),

               #### input for total free ammonia concentration ----
               sliderInput("TOTNH3_mg_N_L",
                           label = p(HTML("Free Ammonia Concentration (mg nitrogen/L)"),
                                     style = "font-size: 16px"),
                           min = 0.00,
                           max = 10.00,
                           value = 0.00,
                           step = 0.05),

               br(),

               #### input for total free bromine concentration ----
               sliderInput("TOTBr_mg_Cl_L",
                           label = p(HTML("Free Bromine Concentration (mg chlorine/L)"),
                                     style = "font-size: 16px"),
                           min = 0.00,
                           max = 10.00,
                           value = 0.00,
                           step = 0.05),

               br(),

               #### input for total free chlorine concentration ----
               sliderInput("TOTCl_mg_Cl_L",
                           label = p(HTML("Free Chlorine Concentration (mg chlorine/L)"),
                                     style = "font-size: 16px"),
                           min = 0.00,
                           max = 10.00,
                           value = 0.00,
                           step = 0.05),

               br(),

               #### input for total orthosilicate concentration ----
               sliderInput("TOTSiO2_mg_SiO2_L",
                           label = p(HTML("Total Dissolved Orthosilicate Concentration (mg silica/L)"),
                                     style = "font-size: 16px"),
                           min = 0.00,
                           max = 50.00,
                           value = 0.00,
                           step = 0.25),

               br()
               ),

  ### main panel ----
  mainPanel(width = 8,

            wellPanel(

              #### simulated plot title ----
              h3(HTML("<b>Buffer Intensity Simulation Plot</b>")),

              fluidRow(
                ##### button to download simulation data ----
                downloadButton("downloadData",
                               label = "Simulation Data Download (.csv file)"),


                ##### button to download plot ----
                downloadButton("downloadPlot",
                               label = "Simulation Plot Download (.png file)")
              ),

              br(),

              ##### produce buffer intensity plot ----
              plotly::plotlyOutput("plot", height = "600px"),

              br(),

              ##### generate slider for pH change selection ----

              # render slider to select desired pH change
              uiOutput("pH_delta"),

              # create output text for required strong acid or strong base added for pH change
              h4(HTML(paste0("<b>", textOutput("pH_delta_acid_base"), "</b>")))
              ),

            wellPanel(

              h3(HTML("<b>Buffer Intensity Simulation Summary Tables</b>")),

              fluidRow(

                column(6,

                       #### simulated condition table title ----
                       h4(HTML("<b>Simulated Condition Summary</b>"))
                       ),

                column(6,
                       ##### pK summary table title ----
                       h4(HTML(paste0("<b>", textOutput("pK_sim"), "</b>")))
                       )
                ),

              fluidRow(

                column(6,
                       ##### button to download simulated condition table ----
                       downloadButton("downloadIC",
                                      label = "Simulated Condition Table Download (.cvs file)")
                ),

                column(6,
                       ##### button to download pK summary table ----
                       downloadButton("downloadpK",
                                      label = HTML("pK Table Download (.cvs file)"))
                       )
                ),

              fluidRow(

                column(6,
                       ##### generate simulated condition table ----
                       DT::dataTableOutput("ic")
                       ),

                column(6,
                       ##### generate pK summary table ----
                       DT::dataTableOutput("pk")
                       )
                )
              )
            )
  ),
))

# code function call ----
shinyApp(ui = ui, server = server)