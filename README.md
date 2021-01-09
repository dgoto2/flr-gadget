# FLR-Gadget
A management strategy evaluation (MSE) framework using FLR and Gadget

## Description 
FLR-Gadget is a management strategy evaluation (MSE) framework using FLR (The Fisheries Library in R) mse (https://github.com/flr/mse) with an R package of customized Gadget (Globally applicable Area Disaggregated General Ecosystem Toolbox, https://github.com/Hafro/gadget2), GadgetR (https://github.com/REDUS-IMR/gadget), as an operating model (OM). This framework is designed to run single and multi- species MSEs. The OM can be age- or length- based. The framework can perform short-cut and full MSEs. A4A (Assessment for All, https://github.com/flr/FLa4a) statistical catch-at-age model and SAM (State-space Assessment Model, https://github.com/flr/FLSAM) are implemented as an assessment model. 

## Prerequisites
Install the following packages:
```
  devtools::install_github("REDUS-IMR/gadget", ref="gadgetr")
  install.packages(c("mse", "FLa4a", "FLash", "FLAssess", "ggplotFL", "FLBRP", "FLCore"), repos="http://flr-project.org/R")
  # The framework requires Version 2.0 or later of mse.
  install.packages(c("dplyr", "MASS", "filelock", "copula","triangle", "coda"))  
```

## Acknowledgements
This work has been conducted as part of Institute of Marine Research/Havforskningsinstituttet's (https://www.hi.no/en) REDUS (Reducing uncertainty in stock assessment) Project (http://redus.no/).
