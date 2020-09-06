# Introduction
geeBFASTmonitor is a reimplementation of R bfastmonitor with only some slight variations.

# BFASTmonitor
BFASTmonitor is a variation of bfast for monitoring purposes. Since geeBFASTmonitor is mostly a reimplementation, the [original documentation](https://www.rdocumentation.org/packages/bfast/versions/1.5.7/topics/bfastmonitor) is still relevant.

# Running geeBFASTmonitor

It's as easy as including :

    var engine = require('users/andreim/geeMonitor:monitor.js');

in your gee script. Afterwards, you'll have access to the function ```engine.bfastMonitor(roi,historyStart,historyEnd,monitoringStart,monitoringEnd,h,period,alpha,magnitudeThreshold,harmonics)```

A typical way of using the function is to save its result, which is an dictionary which includes:

var result = engine.bfastMonitor(roi,historyStart,historyEnd,monitoringStart,monitoringEnd,h,period,alpha,magnitudeThreshold,harmonics)

 - bfastResults: ImageCollection
 - residuals:
 - predictedValues:
 - mresiduals:

The most important of which is the bfastResults ImageCollection, with the following layers:

 - moRes
 - timeCnk2
 - Cnk
 - band_metrics_masked

# Applications

If you want to play with the algorithm interactively, please check:

https://andreim.users.earthengine.app/view/bfastmonitor 