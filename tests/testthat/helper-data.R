quiet <- function(expr){
  utils::capture.output(value <- suppressMessages(expr))
  value
}

sim_clf_data <- function(seed=42){
  set.seed(seed)
  df <- data.frame(rbind(t(rmultinom(60, 300, c(.21,.6,.12,.38,.099))),
                         t(rmultinom(60, 300, c(.001,.6,.42,.58,.299))),
                         t(rmultinom(60, 300, c(.21,.6,.12,.38,.099)))))
  md <- data.frame(study=factor(rep(c("S1", "S2", "S3"), each=60)),
                   disease=factor(rep(rep(c("Case", "Control"), each=30), 3)),
                   subject=factor(rep(1:90, each=2)),
                   age=rep(1:60, 3))
  df[md$disease=="Case", 1] <- df[md$disease=="Case", 1] * 3
  list(df=df, md=md)
}
