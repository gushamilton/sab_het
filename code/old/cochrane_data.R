library(tidyverse)


d  |> 
filter(RCT == "Yes") 

load("data/CDSR.Rdata")

colnames(data)
str(data)
data  |> 
filter(total1 > 500)  |> 
filter(RCT == "yes")  |> 
select(effect.es)  |> 
filter(abs(effect.es) < 10)  |> 
ggplot(aes(x = effect.es)) +
geom_histogram() 

data |>
filter(total1 > 1000) |>
pull(effect.es) |>
abs() |>
quantile(probs = c(0.025,0.5, 0.975))
ex