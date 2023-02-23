
## start with CPI, make it stationary otherwise who knows what the results are telling us!
cpi = read.csv('CPIAUCNS.csv')

plot(cpi$CPIAUCNS,type='l',ylab='CPI',xlab='Months since Jan 1913',main='CPI Plot')
acf(cpi$CPIAUCNS,main='Monthly CPI ACF')
## obviously raw data looks non-stationary,
## so try the first difference of cpi to get stationarity
diff1_cpi = diff(cpi$CPIAUCNS)
plot(diff1_cpi,type='l',ylab='CPI 1st Difference',xlab='Months since Jan 1913',main='CPI 1st Difference Plot')
acf(ts(diff1_cpi),main='CPI 1st Difference ACF')
## still looks pretty non-stationary, let's use an adf test here, and actually
## let's go back and test the original data, too, just to be sure on our
## path of logic here
library(tseries)
adf.test(ts(cpi$CPIAUCNS))
adf.test(ts(diff1_cpi))
## soooo looks like we reject the null that the first difference is non-stationary
## all the way up to the 1% level, so we're sticking with that! Nice
## now just create a date vector that corresponds to the number of values we have here,
## since eventually that'll come in handy (just shave off the first value of original date vector)
cpi_dates = cpi$DATE[-1]

## now move on to Fed Funds rate
fedfunds = read.csv('FEDFUNDS.csv')

plot(fedfunds$FEDFUNDS,type='l',ylab='Federal Funds Rate',xlab='Months since Jul 1954',main='Federal Funds Rate Plot')
acf(fedfunds$FEDFUNDS,main='Federal Funds Rate ACF')
adf.test(ts(fedfunds$FEDFUNDS))
## so plot of ACF seems very non-stationary, but ADF says stationary at the 10% level,
## kind of makes sense looking at the raw data. variance isn't crazy different at any point,
## and obviously it is pretty well bounded below 20(% mind you) and above (for the most part) 0(%).
## however, perhaps we can achieve more certainty with a 1st difference? let's check:
diff1_fedfunds = diff(fedfunds$FEDFUNDS)
plot(diff1_fedfunds,type='l',ylab='Federal Funds Rate 1st Difference',xlab='Months since Jul 1954',main='Federal Funds Rate 1st Difference Plot')
acf(ts(diff1_fedfunds),main='Federal Funds Rate 1st Difference ACF')
adf.test(diff1_fedfunds)
## raw plot actually looks less stationary in terms of variance, but obviously
## way, way more stationary in terms of mean. ACF is more convincing and so is ADF test
## so we'll claim this is our stationary time series for now.
## now make the fedfunds date vector.
fedfunds_dates = fedfunds$DATE[-1]


## now i gotta match up the length of the cpi data to the fedfunds data using...
## just realized I don't need the date vectors for this, I can just shave off
## the first X many elements of the longer vector where X = length of the shorter one - length
## of the longer one...

#diff1_cpi_abb = diff1_cpi[length(diff1_cpi) - length(diff1_fedfunds):length(diff1_cpi)]
## jk let's go with what I originally planned
start_index = which(cpi_dates == fedfunds_dates[1])-1
diff1_cpi_abb = diff1_cpi[start_index:length(diff1_cpi)]

## (the abb is short for abbreviated)
## now I can just do the ADL model can't I?
## Maybe first let's just take a look at their correlation
plot(diff1_cpi_abb,diff1_fedfunds)
## well that wasn't very informative...


## so structure of this analysis is different than I thought when I first set this stuff up,
## but yeah I've set up the data well so far, just need to go through a few more
## steps to find out what the effects of innovations in the interest rate are on inflation
## (at least from an ARDL model --  **not** an ADL btw, should be called ARDL for consistency's sake)

library(forecast)
## steps 1 and 2: find optimal lag level for 1st diff of inflation and 1st diff for fedfunds for AR model
## (only test up to 15 lags, didn't think it was super necessary to test ruther than this)
num_tests = 1:15
cpi_aic = c()
fedfunds_aic = c()
for (lag in num_tests) {
  cpi_test_model = arima(diff1_cpi_abb, order = c(lag, 0, 0))
  cpi_aic[lag] = cpi_test_model[['aic']]
  
  fedfunds_test_model = arima(diff1_fedfunds, order = c(lag, 0, 0))
  fedfunds_aic[lag] = fedfunds_test_model[['aic']]
}
cat('OPTIMAL AR VALUE FOR 1ST DIFF OF CPI:', which(cpi_aic==min(cpi_aic)))
cat('OPTIMAL AR VALUE FOR 1ST DIFF OF FEDFUNDS:', which(fedfunds_aic==min(fedfunds_aic)))

## step 3: obtain residuals from optimal AR model on fedfunds, make sure they're white noise
fedfunds_model = arima(diff1_fedfunds, order = c(13, 0, 0))
fedfunds_res = fedfunds_model[['residuals']]
print(Box.test(fedfunds_res,type='Ljung'))
## comes back with 95.4% p-value, can't reject H0 that they're white noise (all good to go!)
## step 4: "filter" inflation by multiplying it by the lag polynomial estimated for fedfunds
fedfunds_coef = fedfunds_model[['coef']]
filt_diff1_cpi = fedfunds_coef[1] + fedfunds_coef[2]*lag(ts(diff1_cpi_abb)) + fedfunds_coef[3]*lag(ts(diff1_cpi_abb),k=2) +
  fedfunds_coef[4]*lag(ts(diff1_cpi_abb),k=3) + fedfunds_coef[5]*lag(ts(diff1_cpi_abb),k=4) + fedfunds_coef[6]*lag(ts(diff1_cpi_abb),k=5) +
  fedfunds_coef[7]*lag(ts(diff1_cpi_abb),k=6) + fedfunds_coef[8]*lag(ts(diff1_cpi_abb),k=7) + fedfunds_coef[9]*lag(ts(diff1_cpi_abb),k=8) +
  fedfunds_coef[10]*lag(ts(diff1_cpi_abb),k=9) + fedfunds_coef[11]*lag(ts(diff1_cpi_abb),k=10) + fedfunds_coef[12]*lag(ts(diff1_cpi_abb),k=11) +
  fedfunds_coef[13]*lag(ts(diff1_cpi_abb),k=12) + fedfunds_coef[14]*lag(ts(diff1_cpi_abb),k=13)
time = 1:length(diff1_cpi_abb)
plot(time,diff1_cpi_abb,type='l',col='red')
lines(time[1:length(filt_diff1_cpi)],filt_diff1_cpi,col='blue')
## well that just looks odd, but those plots aren't really part of the analysis, I was just curious,
## step 5: look at cross-correlogram of "filtered" cpi and fedfunds residuals, the lags at which
## we see spikes are the lags we should use in our exogenous inclusion of the fedfunds residuals
## in our final regression, and those coefficients are what allow us to find the effect changes in
## interest rates have on inflation.
ccf(filt_diff1_cpi,fedfunds_res,main='CCF of Filtered CPI 1st Difference and Federal Funds Rate Model Residuals')
## thank the Lord Almighty there is a ccf function, I was this --> || close to doing it manually.
## appears that we have quite a few spikes, but the biggest ones are at lags 3, 8, 13, and 18.
## all that's left now is to run the AR(14) of the 1st difference of CPI with the 1st, 2nd,
## and 8th lags of fedfunds residuals as exogenous variables:
exo_var = cbind(lag(ts(fedfunds_res),k=3), lag(ts(fedfunds_res),k=8), lag(ts(fedfunds_res),k=13), lag(ts(fedfunds_res),k=18))
exo_var = exo_var[16:821,1:4]
model = arima(diff1_cpi_abb[1:806], order = c(14, 0, 0), xreg=exo_var)
model_coef = model[['coef']]
model_coef = model[['coef']][(length(model_coef)-2):length(model_coef)]
## print the coefficients of the fedfunds residuals + its lags
print(model_coef)
## get p-values
model_coef_var = model[['var.coef']]
model_coef_var = model_coef_var[(length(model_coef)-2):length(model_coef),(length(model_coef)-2):length(model_coef)]
model_coef_t_denom = sqrt(c(model_coef_var[1,1], model_coef_var[2,2], model_coef_var[3,3])/model[['nobs']])
model_t_stat = model_coef / model_coef_t_denom
model_p_vals = pt(abs(model_t_stat),df=(model[['nobs']]-length(model[['coef']])-1),lower.tail=FALSE)*2
## print p-values
print(model_p_vals)
## turns out to be a (mostly) negative and extremely significant effect of fedfunds innovations
## on inflation! Who would've thought!

## I wanted to see the "impulse response function," i.e. a plot of CPI as a function of the fedfunds residuals
print(fedfunds_coef)
print(1-fedfunds_coef)
cpi_impulse = fitted(model)
plot(fedfunds_res,cpi_impulse,type='l')




##################
## bad news: it appears the steps followed above were not quite what Ender recommends,
## below is what I /believe/ he recommends instead.

## Step 1: Estimate the optimal AR model for the intervention sequence (diff1_fedfunds),
## then store the residuals from this estimation.
## We in fact already did this step, so we can just reuse the content from above 
## (the fedfunds_model and fedfunds_res variables)
## Step 2: Identify plausible candidates for the optimal number of lags of the intervention sequence
## in an ARX model of cpi on fedfunds.
## This is done by a step we actually also already did above, where we "filter" the
## cpi sequence and look at the cross-correlogram between filt_diff1_cpi and fedfunds_res.
## We find the optimal lags are 3, 8, 13, and 18.
## Step 3: Identify plausible candidates for the optimal number of lags of the cpi sequence
## to use in the previously mentioned ARX model of cpi on fedfunds.
## Now /this/ step we have yet to do; regress diff1_cpi on the step 2 lags of diff1_fedfunds
## /without/ an intercept. We can then look at the ACF of the residuals from this estimation
## (as long as the series is white noise), and this will tell us which lags to use for cpi in the ARX model.
fedfunds_lags = cbind(lag(ts(diff1_fedfunds),k=3), lag(ts(diff1_fedfunds),k=8), lag(ts(diff1_fedfunds),k=13), lag(ts(diff1_fedfunds),k=18))[16:821,1:4]
AL_model = lm(diff1_cpi_abb[1:806] ~ fedfunds_lags - 1)
AL_res = AL_model[['residuals']]
print(Box.test(AL_res,type='Ljung'))
acf(AL_res,lag.max=60,main='ACF of Residuals from AR-X Model of 1st Difference CPI on\nOptimal Federal Funds Rate Lags')
pacf(AL_res,lag.max=60,main='PACF of Residuals from AR-X Model of 1st Difference CPI on\nOptimal Federal Funds Rate Lags')
## so... the residuals appear not to be white noise... a bit disheartening. Along with this,
## the ACF indicates that there is some moving average aspect to them as well (damped oscillations).
## Not exactly sure what to do from here! We can proceed with the full estimation if we want,
## we can come up with some form of lag polynomial (moving average component) for the residual
## process here and include that in our final model... But it's very much not white noise,
## and it doesn't seem like there are any obvious solutions to that from Ender. He seems
## to indicate that going down the moving average addition is the way to solve for this,
## and if we want to do that we'll say use up to the second lag of the residuals (MA coefficient is 2).
## We also need to use this ACF (and newly added PACF) plot to deduce the number of autoregressive components for cpi.
## It appears (definitely putting my eye-test skills to work here) that the ACF indicates
## an MA order of 2, and the PACF indicates an AR components at lags 1, 2, 10, 11, 24, 27, 
## It's entirely possible these are silly choices, but let's just try them out for now.
## Step 4: Simply combine the lag levels we've found in the final regression for ARIMAX of cpi on fedfunds;
final_exog_var = cbind(lag(ts(diff1_cpi_abb)), lag(ts(diff1_cpi_abb),k=2), lag(ts(diff1_cpi_abb),k=10), 
                       lag(ts(diff1_cpi_abb),k=11), lag(ts(diff1_cpi_abb),k=24), lag(ts(diff1_cpi_abb),k=27))[27:821,1:6]
final_exog_var = cbind(final_exog_var, fedfunds_lags[1:795,1:4])
print('--------------------------------')
final_model = arima(diff1_cpi_abb[3:797], order = c(0,0,2), xreg=final_exog_var,optim.control = list(maxit=200))
print(final_model[['coef']])
print(final_model[['aic']])
final_model = arima(diff1_cpi_abb[3:797], order = c(0,0,2), seasonal=list(order=c(0,0,12)),xreg=final_exog_var,optim.control = list(maxit=1000))
print(final_model[['coef']])
print(final_model[['aic']])
final_model = arima(diff1_cpi_abb[3:797], order = c(0,0,1), seasonal=list(order=c(0,0,12)),xreg=final_exog_var,optim.control = list(maxit=1000))
print(final_model[['coef']])
print(final_model[['aic']])
final_model = arima(diff1_cpi_abb[3:797], seasonal=list(order=c(0,0,1),period=12),xreg=final_exog_var,optim.control = list(maxit=200))
print(final_model[['coef']])
print(final_model[['aic']])
## tried out a few iterations of including a seasonal aspect, seems like based on AIC the third option is the best,
## which is nice because it by far makes the most sense out of the 4 of the models! The fedfunds coefficients
## are all /actually/ negative!

### SECOND METHOD ###

## process here is just to regress the final model starting at zero lags for each of CPI and fedfunds,
## then adding in lags up to some undetermined level and comparing by AIC. We'll test up to 15 lags
## for now for each of CPI and fedfunds, since it seems there are some issues when I tell it to go above that
## in a lot of cases.

print('-------------------------------------------------------------------')
aic_vec = c()
cpi_lag_vec = c()
ff_lag_vec = c()
max_lags = 0:8
for (cpi_lag in max_lags) {
  ff_exog = diff1_fedfunds
  for (ff_lag in max_lags) {
    if (ff_lag == 1) {
      new_col = diff1_fedfunds[(ff_lag+1):length(diff1_cpi_abb)]
      ff_exog = ff_exog[1:(length(diff1_cpi_abb)-ff_lag)]
      ff_exog = cbind(ff_exog, new_col)
      #print(head(ff_exog))
      #colnames(ff_exog) = ff_exog_colnames
    }
    if (ff_lag > 1) {
      new_col = diff1_fedfunds[(ff_lag+1):length(diff1_cpi_abb)]
      ff_exog = ff_exog[1:(length(diff1_cpi_abb)-ff_lag),]
      ff_exog = cbind(ff_exog, new_col)
      #print(head(ff_exog))
    }
    #print(head(ff_exog))
    test_model = arima(diff1_cpi_abb[1:(length(diff1_cpi_abb)-ff_lag)], order = c(cpi_lag,0,0),
                       seasonal=list(order=c(0,0,12)), xreg = ff_exog, method='ML', optim.control = list(maxit=1000))
    #cat('AIC FOR CPI_LAG = ', cpi_lag, ', FEDFUNDS_LAG = ', ff_lag, ':', test_model[['aic']])
    aic_vec = c(aic_vec, test_model[['aic']])
    cpi_lag_vec = c(cpi_lag_vec, cpi_lag)
    ff_lag_vec = c(ff_lag_vec, ff_lag)
    #cat('cpi_lag = ', cpi_lag,', ff_lag = ', ff_lag)
    #print('')
  }
}

min_aic = which(aic_vec==min(aic_vec))
min_cpi_lag = cpi_lag_vec[min_aic]
min_ff_lag = ff_lag_vec[min_aic]
cat('OPTIMAL ARDL MODEL HAS AIC ', min(aic_vec), ', AND THE NUMBER OF LAGS ON CPI AND THE FEDERAL FUNDS
    RATE ARE ', min_cpi_lag, 'AND ', min_ff_lag, ', RESPECTIVELY.')
## when max lags were set to 8, we got an 865.2145 AIC with cpi_lag = 14 and ff_lag = 15,
## which is much worse than the first method of construction, so I tried up to 20 (before possibly
## including some seasonal aspect to the model, since I imagine that'd make a huge difference).
## soooo turns out we were only including the 15th lag of fedfunds, not all lags up to 15, so that's unfortunate.
## I've adjusted it now to include all lags up to max_lags in the final loop, but turns out
## 15 is waaaay too computationally intensive, so I had to widdle it all the way down to just 8,
## in which case we get AIC = 906.082, cpi_lag = 7, and ff_lag = 8. Not so great! (with seasonality btw)

## well, max lags being set to 20 didn't work! too computationally intensive. Instead we just go
## straight to seasonal addition. We get an AIC of 809.5225 with cpi_lag = 10 and ff_lag = 9 (I
## had to decrease the max lags to 10 due to computational demands again! There must be a way around
## that problem...). I tried up to 12 lags for each and got an AIC of 793.1238 with cpi_lag = 10 (still)
## and ff_lag = 12. It seems like that 13th and 18th lag, and possibly the moving average aspects,
## are fairly important for the model. I won't add them though! I don't think it'll be enough of
## an improvement on this method to justify running another 15-20 minute double for loop (which
## is probably the root of my issues here, btw).


### FORECASTING ###

method1_model = arima(diff1_cpi_abb[3:797], order = c(0,0,1), seasonal=list(order=c(0,0,12)),xreg=final_exog_var,optim.control = list(maxit=1000))
coef1 = method1_model[['coef']]
for (rep in 1:10) {
  cat('FORECAST ', rep, '/10 FOR CHANGE IN CPI IN MARCH 2022 BY METHOD 1 MODEL: ', coef1[1]*mean(rnorm(200,mean=0,sd=1)) + coef1[2]*mean(rnorm(200,mean=0,sd=1)) + coef1[3]*mean(rnorm(200,mean=0,sd=1)) + 
        coef1[4]*mean(rnorm(200,mean=0,sd=1)) + coef1[5]*mean(rnorm(200,mean=0,sd=1)) + coef1[6]*mean(rnorm(200,mean=0,sd=1)) + coef1[7]*mean(rnorm(200,mean=0,sd=1)) + coef1[8]*mean(rnorm(200,mean=0,sd=1)) + 
        coef1[9]*mean(rnorm(200,mean=0,sd=1)) + coef1[10]*mean(rnorm(200,mean=0,sd=1)) + coef1[11]*mean(rnorm(200,mean=0,sd=1)) + coef1[12]*mean(rnorm(200,mean=0,sd=1)) + coef1[13]*mean(rnorm(200,mean=0,sd=1)) + 
        coef1[14] + coef1[15]*diff1_cpi_abb[length(diff1_cpi_abb)] + coef1[16]*diff1_cpi_abb[length(diff1_cpi_abb)-1] + 
        coef1[17]*diff1_cpi_abb[length(diff1_cpi_abb)-9] + coef1[18]*diff1_cpi_abb[length(diff1_cpi_abb)-10] + 
        coef1[19]*diff1_cpi_abb[length(diff1_cpi_abb)-23] + coef1[20]*diff1_cpi_abb[length(diff1_cpi_abb)-26] + 
        coef1[21]*0.25 + coef1[22]*diff1_fedfunds[length(diff1_fedfunds)-7] + coef1[23]*diff1_fedfunds[length(diff1_fedfunds)-12] + 
        coef1[24]*diff1_fedfunds[length(diff1_fedfunds)-17])
  print('')
}
## make correct ff_lag dataframe...
ff_exog = diff1_fedfunds
min_ff_lag = 8
for (ff_lag in 1:min_ff_lag) {
  if (ff_lag == 1) {
    new_col = diff1_fedfunds[(ff_lag+1):length(diff1_cpi_abb)]
    ff_exog = ff_exog[1:(length(diff1_cpi_abb)-ff_lag)]
    ff_exog = cbind(ff_exog, new_col)
  }
  if (ff_lag > 1) {
    new_col = diff1_fedfunds[(ff_lag+1):length(diff1_cpi_abb)]
    ff_exog = ff_exog[1:(length(diff1_cpi_abb)-ff_lag),]
    ff_exog = cbind(ff_exog, new_col)
  }
}
min_cpi_lag = 7
method2_model = arima(diff1_cpi_abb[1:(length(diff1_cpi_abb)-ff_lag)], order = c(min_cpi_lag,0,0),
                      seasonal=list(order=c(0,0,12)), xreg = ff_exog, method='ML', optim.control = list(maxit=1000))
coef2 = method2_model[['coef']]
for (rep in 1:10) {
  cat('FORECAST ', rep, '/10 FOR CHANGE CPI IN JANUARY 2022 BY METHOD 2 MODEL: ', coef2[1]*diff1_cpi_abb[length(diff1_cpi_abb)] +
        coef2[2]*diff1_cpi_abb[length(diff1_cpi_abb)-1] + coef2[3]*diff1_cpi_abb[length(diff1_cpi_abb)-2] +
        coef2[4]*diff1_cpi_abb[length(diff1_cpi_abb)-3] + coef2[5]*diff1_cpi_abb[length(diff1_cpi_abb)-4] +
        coef2[6]*diff1_cpi_abb[length(diff1_cpi_abb)-5] + coef2[7]*diff1_cpi_abb[length(diff1_cpi_abb)-6] + 
        coef2[8]*mean(rnorm(200,mean=0,sd=1)) + coef2[9]*mean(rnorm(200,mean=0,sd=1)) + coef2[10]*mean(rnorm(200,mean=0,sd=1)) + 
        coef2[11]*mean(rnorm(200,mean=0,sd=1)) + coef2[12]*mean(rnorm(200,mean=0,sd=1)) + coef2[13]*mean(rnorm(200,mean=0,sd=1)) + 
        coef2[14]*mean(rnorm(200,mean=0,sd=1)) + coef2[15]*mean(rnorm(200,mean=0,sd=1)) + coef2[16]*mean(rnorm(200,mean=0,sd=1)) + 
        coef2[17]*mean(rnorm(200,mean=0,sd=1)) + coef2[18]*mean(rnorm(200,mean=0,sd=1)) + coef2[19]*mean(rnorm(200,mean=0,sd=1)) + 
        coef2[20] + coef2[21]*0.25 + coef2[22]*diff1_fedfunds[length(diff1_fedfunds)-1] + coef2[23]*diff1_fedfunds[length(diff1_fedfunds)-1] +
        coef2[24]*diff1_fedfunds[length(diff1_fedfunds)-1] + coef2[25]*diff1_fedfunds[length(diff1_fedfunds)-1] + coef2[26]*diff1_fedfunds[length(diff1_fedfunds)-1] + 
        coef2[27]*diff1_fedfunds[length(diff1_fedfunds)-1] + coef2[28]*diff1_fedfunds[length(diff1_fedfunds)-1])
  print('')
}
print(coef2[(length(coef2)-8):length(coef2)])