# Start the R

In the Data folder and copy the script lines into that interface.


```
cd Data 
R
system( "ls" )
```
Run Script for the Stroke analysis directly in R:
```
system("Rscript /Users/ju5263ta/Github/Monocytes/Scripts_stroke/240605_Sublime_Github.R")
```

All output will be saved in a folder: TODAYS-DATE_Stroke_Results

Questions:
- Can this run like this for everyone? 
-


See if i still need this:

## Check imputation & Questions --------------------------------------------------

# Plots for imputation--------------------------------------------------
# Check convergence
#plot(finaldata)
# kernel density plots for imputed values 
#    densityplot(data_imputed)

### which one to look at, shows the spread of the imputated data
densityplot(finaldata, data=~dCT_Value | Subpopulation)
densityplot(data_imputed, data=~dCT_Value | Sex)

#choose imputation number 1
# better: choose the mean from 20 impoutations
check1<-lm(dCT_Mean~Gene+Age+Sex+Subpopulation+SampleID, data=finaldata, subset=.imp==1)
#plot(check1)
dataFINAL_imp <- filter(finaldata, .imp==1)

## different approach but oly works on the rE!
dataFINAL$rE_Value <- 2^((-1)*dataFINAL$dCT_Mean)
# with the original  data but without individuals!
glm(rE_Value~Gene+Age+Sex+Subpopulation, family = "gaussian", data = dataFINAL) %>% tab_model()

# Questions:
# how is the summary interpretated? estimate? SD error? p value?
# especially in this one! what does it mean for the one thats not shown of the levels in factors?
# lm(dCT_Mean~Gene+Age+Sex+Subpopulation, data = dataFINAL) %>% sjPlot::tab_model()
# Gene vs Gene1
# exclude failed genes / neg. controls?