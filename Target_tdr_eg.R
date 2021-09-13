# Target diagram examples
# From tdr package

## `data` must include these columns: nrmse, nmbe, sdm, sdo

targetDiagram(data, class = "",
              xlab = expression("RMSEc" %.% "sign(" * sigma^"*" * ")"),
              ylab = "MBE",
              auto.key = list(space = "right"),
              default.scales = list(cex = 0.6),
              scales = list(),
              type = "quantile",
              cuts = seq(0.25, 1, 0.25),
              par.settings = tdTheme(),
              ...)

target_diagram(data,
               xlab = expression("RMSEc"%.%"sign("*sigma^"*"*")"),
               ylab = 'MBE', 
               type = 'quantile', cuts = seq(0.25, 1, .25),
               ...)

tdTheme(pch = 21:25,
        col.points = brewer.pal(n= 8, 'Dark2'),
        fill = brewer.pal(n= 8, 'Pastel2'),
        cex = 1.1,
        ...)


data(modelEx)

## Compute statistics
errModel <- applyStats(pvModels, pvObs)

## Display results
## Default settings use type = 'quantile'
targetDiagram(errModel, groups = model)

target_diagram(errModel, fill = 'model')

## whose breaks can be defined with 'cuts'
targetDiagram(errModel, groups = model,
              type = 'quantile',
              cuts = seq(0, 1, .1))

target_diagram(errModel, fill = 'model',
               type = 'quantile',
               cuts = seq(0, 1, .1))

## Alternatively, with type = 'at'
## one can define manually the breaks
targetDiagram(errModel, groups = model,
              type = 'at',
              cuts = seq(0, .1, .02))

target_diagram(errModel, fill = 'model',
               type = 'at',
               cuts = seq(0, .1, .02))