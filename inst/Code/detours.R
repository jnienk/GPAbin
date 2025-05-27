library(tidyverse)
library(tourr)
library(mulgar)
library(Rtsne)
library(uwot)
library(crosstalk)
library(plotly)
library(detourr)

#a <- missmi(implist)|> DRT() |> GPAbin() |> biplFig()

#GPAbin

CLP.col <- vector("character",length(a$levels[[1]])) #create a vector with the same length as the number of CLPs
CLP.col <- c(rep("#61223b",length(a$levels[[1]])))

samps.dat <- tibble(samps.1 = a$Z.GPAbin[,1],
                    samps.2 = a$Z.GPAbin[,2])
#sample colour will be specified according to cluster later in the add_trace()

#changes
CLPs.dat <- tibble(name = a$levels[[1]],
                   CLP.1 = a$CLP.GPAbin[,1],
                   CLP.2 = a$CLP.GPAbin[,2],CLP.pch = CLP.pch,
                   CLP.col = CLP.col)
##MCA and data setup##

#gpabin_df <- data.frame(GPAbinX = a$Z.GPAbin[, 1],
#                      GPAbinY = a$Z.GPAbin[, 2])
#comb_df <- bind_cols(apply(comp.sim,2,as.numeric), gpabin_df)

tib_comp <- as.tibble(apply(comp.sim,2,as.numeric))

detour_plot <- detour(tib_comp, tour_aes(
  projection = starts_with("V"))) |>
  tour_path(grand_tour(2),
            max_bases=50, fps = 60) |>
  show_scatter(alpha = 0.7, axes = FALSE,
               width = "100%", height = "450px")

GPAbin_biply <- plot_ly() %>% 
  layout(                   
    yaxis = list(scaleanchor = "x",
                 scaleratio = 1), showlegend = FALSE
  ) %>%
  add_trace(
    type = "scatter",
    mode = "markers",
    x = samps.dat$samps.1, 
    y = samps.dat$samps.2,marker=list(symbol="circle", color="#b79962"),
    text = c(""),
    textposition = "right",
    hoverinfo = "text"
  ) %>%
  add_trace(type = "scatter",
            mode = "markers",
            x = CLPs.dat$CLP.1, 
            y = CLPs.dat$CLP.2, marker=list(color=CLPs.dat$CLP.col, symbol=CLPs.dat$CLP.pch,
                                            size=12),
            text = CLPs.dat$name,
            textposition = "right",
            hoverinfo = "text"
  )  %>%
  layout(hovermode="x")

bscols(
  detour_plot, GPAbin_biply,
  widths = c(5, 6)
)




##interactive plotting##

plot_ly() %>% 
  layout(                   
    yaxis = list(scaleanchor = "x",
                 scaleratio = 1), showlegend = FALSE
  ) %>%
  add_trace(
    type = "scatter",
    mode = "markers",
    x = samps.dat$samps.1, 
    y = samps.dat$samps.2,marker=list(symbol="circle", color="#b79962"),
    text = c(""),
    textposition = "right",
    hoverinfo = "text"
  ) %>%
  add_trace(type = "scatter",
            mode = "markers",
            x = CLPs.dat$CLP.1, 
            y = CLPs.dat$CLP.2, marker=list(color=CLPs.dat$CLP.col, symbol=CLPs.dat$CLP.pch,
                                            size=12),
            text = CLPs.dat$name,
            textposition = "right",
            hoverinfo = "text"
  )  %>%
  layout(hovermode="x")

##interactive plotting##