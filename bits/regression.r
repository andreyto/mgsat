library(arm)
## Expected values vs average residual plot for logistic regression
## non-mixed model
binnedplot(predict(m,type="response"),resid(m,type="response"))
## without random effects
binnedplot(predict(m,type="response",re.form=NA),resid(m,type="response"))
## with random effects
binnedplot(predict(m,type="response",re.form=NULL),resid(m,type="response"))
