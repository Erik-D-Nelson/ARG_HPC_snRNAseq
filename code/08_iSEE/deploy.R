library("rsconnect")

source("token.R")

options(repos = BiocManager::repositories())
rsconnect::deployApp(
    appFiles = c("app.R", "sce_iSEE.rda", "initial.R"),
    appName = "2022_HPC_ARG",
    account = "libd",
    server = "shinyapps.io"
)

