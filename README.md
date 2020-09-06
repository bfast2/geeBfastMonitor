# Introduction
geeBFASTmonitor is a reimplementation of R bfastmonitor with only some slight variations.

# BFASTmonitor
BFASTmonitor is a variation of bfast for monitoring purposes. Since geeBFASTmonitor is mostly a reimplementation, the [original documentation](https://www.rdocumentation.org/packages/bfast/versions/1.5.7/topics/bfastmonitor) is still relevant.

# Running geeBFASTmonitor

It's as easy as including :

    var engine = require('users/andreim/geeMonitor:monitor.js');

in your gee script. Afterwards, you'll have access to the function ```engine.bfastMonitor(roi,historyStart,historyEnd,monitoringStart,monitoringEnd,h,period,alpha,magnitudeThreshold,harmonics)```

The parameters are similar to the R version of bfastMonitor:

 - roi: region of interest as a Feature
 - historyStart: The starting date of the stable period.
 - historyEnd: The end date of the stable period.
 - monitoringStart: The starting date of the monitoring period.
 - monitoringEnd: The end date of the monitoring period.
 - h: numeric scalar from interval (0,1) specifying the bandwidth relative to the sample size in MOSUM/ME monitoring processes
 - period: 
 - alpha: Significance level of the monitoring (and ROC, if selected) procedure, i.e., probability of type I error.
 - magnitudeThreshold: Threshold for magnitude level for which a change should be considered. Default is 0.
 - harmonics: Order of the harmonic term.

A typical way of using the function is to save its result:

```var result = engine.bfastMonitor(roi,historyStart,historyEnd,monitoringStart,monitoringEnd,h,period,alpha,magnitudeThreshold,harmonics)```

The most important results are:
 - the time of change found in bfastResults.timeCnk2
 - magnitude of change found in bfastResults.Cnk


# Applications

If you want to play with the algorithm interactively, please check:

https://andreim.users.earthengine.app/view/bfastmonitor 