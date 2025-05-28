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

CLP.col <- vector("character",length(a$lvls[[1]])) #create a vector with the same length as the number of CLPs
CLP.col <- c(rep("#61223b",length(a$lvls[[1]])))

CLP.pch <- vector("character",length(a$lvls[[1]])) #create a vector with the same length as the number of CLPs
CLP.pch <- c(rep('x',length(a$lvls[[1]])))

samps.dat <- tibble(samps.1 = a$Z.GPAbin[,1],
                    samps.2 = a$Z.GPAbin[,2])
#sample colour will be specified according to cluster later in the add_trace()


#changes
CLPs.dat <- tibble(name = a$lvls[[1]],
                   CLP.1 = a$CLP.GPAbin[,1],
                   CLP.2 = a$CLP.GPAbin[,2],
                   CLP.col = CLP.col, CLP.pch = CLP.pch)
##MCA and data setup##

#gpabin_df <- data.frame(GPAbinX = a$Z.GPAbin[, 1],
#                      GPAbinY = a$Z.GPAbin[, 2])
#comb_df <- bind_cols(apply(comp.sim,2,as.numeric), gpabin_df)

tib_comp <- as.tibble(apply(compdat,2,as.numeric))
tib_comp <- as.tibble(apply(tib_comp,2, jitter, factor=0.5))

x <- bind_cols(samps.dat, tib_comp)
shared_dat <- SharedData$new(x)

detour_plot <- detour(shared_dat, tour_aes(
  projection = starts_with("V"))) |>
  tour_path(grand_tour(2),
            max_bases=50, fps = 60) |>
  show_scatter(alpha = 0.7, axes = FALSE,
               width = "100%", height = "450px", palette="#b79962")

GPAbin_biply <- plot_ly(shared_dat, x = ~samps.1, 
                        y = ~samps.2) %>% 
  layout(                   
    yaxis = list(scaleanchor = "x",
                 scaleratio = 1), showlegend = FALSE
  ) %>%
  add_trace(
    type = "scatter",
    mode = "markers",
    marker=list(symbol="circle"),
    text = c(""),
    textposition = "right",
    hoverinfo = "text"
  ) |> 
  highlight(on="plotly_selected",
            off="plotly_doubleclick")
  # add_trace(type = "scatter",
  #           mode = "markers",
  #           x = CLPs.dat$CLP.1, 
  #           y = CLPs.dat$CLP.2, marker=list(color=CLPs.dat$CLP.col, symbol=CLPs.dat$CLP.pch,
  #                                           size=12),
  #           text = CLPs.dat$name,
  #           textposition = "right",
  #           hoverinfo = "text"
  # )  %>%
  #layout(hovermode="x")

  
bscols(
  detour_plot, GPAbin_biply,
  widths = c(5, 6)
)
