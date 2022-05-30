# TAVC.seg
Robust multiscale time-average variance estimation for change point detection.

Software accompanying 
> E. T. McGonigle and H. Cho (2022) ["Robust multiscale estimation of time-average variance for time
series segmentation".](https://arxiv.org/abs/2205.11496)

- The main routines are contained in main.R. 

To perform robust TAVC estimation, do the following:

- Source main.R into `R`.
- Read the description for `robust.tavc.est` within main.R.

To perform mean change point detection with the robust TAVC estimation procedure

- Using the multiscale bottom-up MOSUM procedure: install the MOSUM R package and read the description in `mosum.tavc'.
- Using the wild binary segmentation 2 algorithm: install the breakfast R package and read the description in `WBS2.tavc'.

For example,

```{r}

cpt.sig = c(rep(0,200),rep(2,300),rep(4,200),rep(2,300))

set.seed(123)

x = cpt.sig + arima.sim(model = list(ar = 0.5), sd = sqrt(1-0.5^2), n = 1000)
x.m.c = mosum.tavc(x,G = c(30,60,90,150), alpha = 0.05)

x.m.c$cpts

x.w.c = wbs2.tavc(x, min.int.len = 60)

x.w.c$cpts


```

If you have any questions, please contact euan.mcgonigle@bristol.ac.uk
