# TarDis Map - Exploring Target Disease Associations

## Precomputing data

Run `code.R` to regenerate the precomputed dataset (`work.Rda`) and copy it into the `tardis` folder. You may need to change database settings within the code for a new version of TCRD

## Deploying the code

The R code to generate and explore TarDis heatmaps can be deployed as a Shiny
application by executing, from within R,
```
library(shiny)
runApp('tardis')
```

