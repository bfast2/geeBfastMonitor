// BFAST-MONITOR: BETA VERSION 2019

//BFASTMONITOR IS DEVELOPED BY: JAN VERBESSELT (e-mail: jan.verbesselt@wur.nl; WAGENINGEN UNIVERSITY) AND IS IMPLEMENTED IN BFAST PACKAGE, AN R PACKAGE.

//THIS GEE CODE WAS INITIALLY WRITTEN BY ELIAKIM HAMUNYELA IN 2014 WHEN HE WAS WITH THE WAGENINGEN UNIVERSITY (2013 - 2017),
//WITH A LOT OF HELP FROM ERIC ENGLE, DAVID THAU, NOEL GORELICK (GOOGLE). 
//SINCE 2018, ELIAKIM IS NOW WITH THE UNIVERSITY OF NAMIBIA, DEPARTMENT OF GEOGRAPHY.
//DUE TO ERRORS IN THE PREVIOUS IMPLEMENATION, ELIAKIM HAS NOW REVISED THE CODE (2018 -2019).

//THE PREVIOUS VERSION OF THE CODE WAS SHARED WIDELY WITH OTHER RESEARCHERS WHO REQUESTED IT.

//NB:THE USE OF THIS NEW VERSION OF THE CODE REQUIRES CONSENT FROM THE IMPLEMENTER OF IT (ELIAKIM: hamunyelae@unam.na)

//YEAR: 2019

///
//////////////////////
////SET PARAMETERS////
//////////////////////

/*
// Set the region of interest to a point.
//var newperu: lat2 =-69.58134, lat =-69.02241, lng =-12.8900, lng2 =-12.6526;
//var peru: lat2 =-69.58108, lat =-69.02646, lng =-12.88777, lng2 =-12.12106;
//var bolivia: lat2= -62.94023, lat = -61.80149, lng = -18.78287, lng2= -18.06762;
//var mozambique: lat2 =36.13087, lat = 37.56728, lng =-17.71388, lng2 =-16.53403 
//var newmozambique: lat2 =36.12234, lat = 36.5975, lng =-17.71089, lng2 =-16.52972 
var lng = -12.88777;//ymin
var lat = -69.02646; //xmax
var lng2= -12.12106; //ymax
var lat2= -69.58108;//xmin

var region = ee.Geometry.Polygon([[lat2, lng],[lat2, lng2], [ lat, lng2],[lat, lng] ]);
var midLat = lat2-((lat2 - lat)/2);
var midLon = lng2-((lng2-lng)/2);

//var roi = ee.Geometry.Point(midLat, midLon);
var roi = ee.Geometry.Point(-60.00, -14.33);
//var roio = ee.Geometry.Point(-69.211188, -12.220112);
var roio = roi
//var region = roi

//var region = ee.Geometry.Point([lng, lat]);

//var region = geometry2;
//var roi = region;
// set the history and monitoring period

var historyStart = '2013-01-01';
var monitoringEnd = '2018-12-31';
var historyEnd = "2016-12-31"
var monitoringStart = "2017-01-01";


//set Bfast-monitor parameters

var h = 0.25;
var hxh = h
var period = 10;
var alpha = 0.05;
var atAlpha = 1 - alpha;
var magnitudeThreshold = -0.000000000000000000000000000000121
var harmonics = 1;
*/

/////////////////////////////////////////////////////
////CHECK IF USER-DEFINED PARAMETERS ARE OK ////
/////////////////////////////////////////////////////
var bfastMonitor = function bfastMonitor(roi,historyStart,historyEnd,monitoringStart,monitoringEnd,h,period,alpha,magnitudeThreshold,harmonics){

roi = roi||ee.Geometry.Point(-60.00, -14.33)
var roio = roi
historyStart = historyStart||'2013-01-01'
monitoringEnd = monitoringEnd||'2018-12-31'
historyEnd = historyEnd||"2016-12-31"
monitoringStart = monitoringStart||"2017-01-01"
h = h||0.25
var hxh = h
period = period||10
alpha = alpha||0.05
var atAlpha = 1 - alpha;
magnitudeThreshold = magnitudeThreshold||-0.000000000000000000000000000000121
harmonics = harmonics||1

/// Check the h parameter
if (h == 0.25) {
  h = 0;
} else if (h == 0.5){
  h = 1;
}else if (h == 1){
  h = 2;
}else {
  alert("Error: h parameter can either be 0.25, 0.5, or 1");
}

///Check alpha. Alpha be equal or less than 0.05

if (alpha > 0.05) {
  alert("Error: alpha parameter set too large, try alpha = 0.05 or less");
}

var tTable = [0.950, 0.951, 0.952, 0.953, 0.954, 0.955, 0.956, 0.957, 0.958, 0.959, 0.960, 0.961, 0.962, 0.963, 0.964, 0.965, 0.966, 0.967, 0.968, 0.969,
0.970, 0.971 ,0.972, 0.973, 0.974, 0.975, 0.976, 0.977, 0.978, 0.979, 0.980, 0.981, 0.982, 0.983, 0.984, 0.985, 0.986, 0.987, 0.988, 0.989, 0.990 ,0.991 ,0.992,
0.993 ,0.994, 0.995, 0.996, 0.997, 0.998, 0.999];

var alphaIndex = tTable.indexOf(atAlpha)
if (alphaIndex < 0 || alphaIndex > 49) {
  alert("Error: for critical values, we only have for alpha parameter ranging from 0.001 to 0.05, try alpha from that range");
}

if (period == 2) {
   period  = 0;
} else if (period  == 4){
   period  = 1;
}else if (period  == 6){
   period  = 2;
}else if (period  == 8){
   period  = 3;
}else if (period  == 10){
   period  = 4;
}else {
  alert("Error: for period parameter, we only have 2, 4, 6, 8,10. Choose one of these values");
}


///////////////////////////////////////////////////////////////////////////
// DEFINE KEY FUNCTIONS FOR PRE=PROCESSING
///////////////////////////////////////////////////////////////////////////


// The dependent variable we are modeling.
var dependent = 'NDMI';

// Make a list of harmonic frequencies to model.
// These also serve as band name suffixes.
var harmonicFrequencies = ee.List.sequence(1, harmonics);

// Function to get a sequence of band names for harmonic terms.
var constructBandNames = function(base, list) {
  return ee.List(list).map(function(i) {
    return ee.String(base).cat(ee.Number(i).int());
  });
};

// Construct lists of names for the harmonic terms.
var cosNames = constructBandNames('cos_', harmonicFrequencies);
var sinNames = constructBandNames('sin_', harmonicFrequencies);

// Independent variables.
var independents = ee.List(['constant'])
  .cat(cosNames).cat(sinNames);
  


// Function to add an NDVI band, the dependent variable.
var addNDVI = function(image) {
  return image
    .addBands(image.normalizedDifference(['B4', 'B5'])
    .rename('NDMI'))
    .float();
};

var addNDVI8 = function(image) {
  return image
    .addBands(image.normalizedDifference(['B5', 'B6'])
    .rename('NDMI'))
    .float();
};

// Functions to add a time band.
var addDependents = function(image) {
  // Compute time in fractional years since the epoch.
  var years = image.date().difference('1970-01-01', 'year');
  var timeRadians = ee.Image(years.multiply(2 * Math.PI)).rename('t');
  var constant = ee.Image(1);
  return image.addBands(constant).addBands(timeRadians.float());
};


// Function that returns the year (with a fractional part to indicate progression within
// a year) that corresponds to the given unix time.
function unixToYear(unixTime) {
  return ee.Number(unixTime).divide(365.25 * 24 * 3600 * 1000).add(1970);
}

// Function to compute the specified number of harmonics
// and add them as bands.  Assumes the time band is present.
var addHarmonics = function(freqs) {
  return function(image) {
    // Make an image of frequencies.
    var frequencies = ee.Image.constant(freqs);
    // This band should represent time in radians.
    var time = ee.Image(image).select('t');
    // Get the cosine terms.
    var cosines = time.multiply(frequencies).cos().rename(cosNames);
    // Get the sin terms.
    var sines = time.multiply(frequencies).sin().rename(sinNames);
    return image.addBands(cosines).addBands(sines);
  };
};

////////////////////////////////////////////////
// LOAD THE IMAGES 
///////////////////////////////////////////////

// Load a collection of Landsat TOA reflectance images.
var landsat5Collection = ee.ImageCollection('LANDSAT/LT05/C01/T1_SR');
var landsat7Collection = ee.ImageCollection('LANDSAT/LE07/C01/T1_SR');
var landsat8Collection = ee.ImageCollection('LANDSAT/LC08/C01/T1_SR');

//mask the  clouds
var landsat5Collection1 = landsat5Collection
  .filterBounds(roi)
  .filterDate(historyStart, monitoringEnd)
  //.map(cloudMaskL457)
     .map(function(image) {return image.mask(image.select('pixel_qa')
                                               .remap([66,130,322,386],[1,1,1,1],0) // 66..386 -> 1 = not masked
                                          )
                                    //.clip(PolygonBuffer)
                        }
      )
   .map(addNDVI);
  
var landsat7Collection1 = landsat7Collection
 .filterBounds(roi)
 .filterDate(historyStart, monitoringEnd)
 // .map(cloudMaskL457)

   .map(function(image) {return image.mask(image.select('pixel_qa')
                                               .remap([66,130,322,386],[1,1,1,1],0) // 66..386 -> 1 = not masked
                                          )
                                    //.clip(PolygonBuffer)
                        }
      )
   .map(addNDVI);
  
var landsat8Collection1 = landsat8Collection
  .filterBounds(roi)
  .filterDate(historyStart, monitoringEnd)
  //.map(maskL8sr)
     .map(function(image) {return image.mask(image.select('pixel_qa')
                                               .remap([66,130,322,386],[1,1,1,1],0) // 66..386 -> 1 = not masked
                                          )
                                    //.clip(PolygonBuffer)
                        }
      )
  .map(addNDVI8);
  

var landsat5Collection01 = landsat5Collection1
.select('B1', 'B2','B3', 'B4','B5', 'B7','NDMI')
.map(function(image){
  return image.rename(['B1', 'B2','B3', 'B4','B5', 'B7','NDMI']);
})

var landsat7Collection01 = landsat7Collection1
.select('B1', 'B2','B3', 'B4','B5', 'B7','NDMI')
.map(function(image){
  return image.rename(['B1', 'B2','B3', 'B4','B5', 'B7','NDMI']);
})

var landsat8Collection01 = landsat8Collection1
.select('B2', 'B3','B4', 'B5','B6', 'B7','NDMI')
.map(function(image){
  return image.rename(['B1', 'B2','B3', 'B4','B5', 'B7', 'NDMI']);
})

// merge the collections and sort them by time
  
var collection_merge = ee.ImageCollection(landsat5Collection01.merge(landsat7Collection01));
var collection_merge2 = collection_merge.merge(landsat8Collection01);

//print (collection_merge2);

var collection_merge3 = collection_merge2
 .sort('system:time_start')
//print (collection_merge3.select('NDMI'));


// subset the collection into history and monitoring period
var histCollection = collection_merge3.filterDate(historyStart, historyEnd);
var moniCollection  = collection_merge3.filterDate(monitoringStart, monitoringEnd);

/////////////////////////
//// HISTORY PERIOD ////
////////////////////////

// Filter to the area of interest, mask clouds, add variables.
var harmonicLandsat = histCollection
  .map(addDependents)
  .map(addHarmonics(harmonicFrequencies));

// The output of the regression reduction is a 4x1 array image.
var harmonicTrend = harmonicLandsat
  .select(independents.add(dependent))
  .reduce(ee.Reducer.linearRegression(independents.length(), 1));

// Turn the array image into a multi-band image of coefficients.
var harmonicTrendCoefficients = harmonicTrend.select('coefficients')
  .arrayProject([0])
  .arrayFlatten([independents]);
  
// 3. CALCULATE COMMULATIVE SUM OF RESIDUALS

// Compute fitted values.
var fittedHarmonic = harmonicLandsat.map(function(image) {
  return image.addBands(
    image.select(independents)
      .multiply(harmonicTrendCoefficients)
      .reduce('sum')
      .rename('fitted'));
});

// Compute the residuals
var computResiduals = function(image){
  return image.addBands(
    image.select('NDMI')
    .subtract(image.select('fitted'))
    .reduce('sum')
    .rename('residual'));

}

var residuals = fittedHarmonic.map(computResiduals);



// Compute squared residuals
var computsquaredResiduals = function(image){
  return image.addBands(
    image.select('residual')
    .pow(2)
    .rename('squaredResidual'));
}

var squaredResiduals = residuals.map(computsquaredResiduals);
var rndmi = residuals.select('NDMI')
var nmdi = rndmi.toArray();
//calculate sigma i.e. sigma = sqrt(sum(e^2)/fm$df.residual)
var residu = residuals.select('residual')
var res = residu.toArray();
var resid = squaredResiduals.select('squaredResidual')
var residarray = resid.toArray();
var imageAxis = 0;
var ResidSum = residarray.arrayReduce("sum", [imageAxis]);
var coefArray = harmonicTrendCoefficients.toArray();
var coefLength = coefArray.arrayLength(imageAxis);
var histtimeseiesLength = residarray.arrayLength(imageAxis);
var dfresiduals = histtimeseiesLength.subtract(coefLength);
var ResidSumnorma = ResidSum.divide(dfresiduals);
var sigma = ResidSumnorma.sqrt();

 /// calculate K, which is N obs in historical sample multiplied by h parameter
var himage = ee.Image(hxh);
var ksize =  himage.multiply(histtimeseiesLength).floor();

// Calculate mosum of residuals for history period 

/// Concatenate two images: zer0 and residuals
var zer0 = ee.Image([0]).toArray().float();
var bandAxis = 1
var fi = zer0.arrayRepeat(bandAxis, 1);
var catImage = ee.ImageCollection([fi,res]);
var concat = catImage.toArrayPerBand(0);
var ImToArr = concat.toArray();

/// ...now calculate Cumsum
var histcumsum = ImToArr.arrayAccum(imageAxis, 'sum');
var histimZer0 = ee.Image(0).int();
var histdubImZer0 = histimZer0.arrayRepeat(imageAxis,ksize.int());
var histfx10 = histdubImZer0.arrayRepeat(bandAxis, 1);
var histverc= histfx10.toArray().float();
var histshiftedCumR = ee.ImageCollection([histverc, histcumsum]);
var histShfCumSumRight1 = histshiftedCumR.toArrayPerBand(imageAxis);
var histShfCumSumRight = histShfCumSumRight1.toArray()
var histshiftedCumL = ee.ImageCollection([histcumsum, histverc]);
var histShfCumSumLeft = histshiftedCumL.toArrayPerBand(imageAxis);
var histShiftDifference = histShfCumSumLeft.subtract(histShfCumSumRight);
var histmsk = histShiftDifference.arrayMask(histShfCumSumLeft);
var histter = ee.Image(0).float().toArray();
var histtesser = histter.arrayRepeat(bandAxis, 1);
var histtessermaske = ee.ImageCollection([histtesser, histmsk]);
var histtessermsk = histtessermaske.toArrayPerBand(imageAxis);

// ...now slash off unwanted first part of Cumsum of residuals to get MOSUM residuals
var histcotx2 = histtessermsk.arrayLength(imageAxis);
var histcumKsizeDiff = histcotx2.subtract(ksize).int();
var histknsiz1 = histimZer0.arrayRepeat(imageAxis, ksize.int());
var histknsiz2 = histknsiz1.arrayRepeat(bandAxis, 1);
var histknsiz =histknsiz2.toArray();

var histdcumn = ee.Image(1).multiply(histcumKsizeDiff);
var histimOnes = ee.Image(1).int();
var histknm1 = histimOnes.arrayRepeat(imageAxis, histdcumn);
var histknm2 = histknm1.arrayRepeat(bandAxis, 1);
var histknm = histknm2.toArray();

var histshmasker = ee.ImageCollection([histknsiz, histknm]);
var histShfmask = histshmasker.toArrayPerBand(imageAxis);
var histResmosum = histtessermsk.arrayMask(histShfmask);

/// ...standardise
var sigmaStanda = sigma.multiply(histtimeseiesLength.sqrt());

var duSigmaStanda = sigmaStanda.arrayRepeat(bandAxis,1);
var bleng = histResmosum.arrayLength(imageAxis);
var xSigmaStanda = duSigmaStanda.arrayRepeat(imageAxis,bleng).toArray();

var histResmosumStanda = histResmosum.divide(xSigmaStanda);


////////////////////////////
///CRITICAL VALUE TABLE ///
////////////////////////////

var criticalVTable = [[[1.227627,1.336231,1.341087,1.341657,1.341825],[1.687323,1.886331,1.899584,1.901299,1.902003],[2.224088,2.704437,2.737148,2.742879,2.745928]],
[[1.230670,1.338791,1.343916,1.344138,1.344391],[1.691601,1.890238,1.903723,1.905166,1.905759],[2.231672,2.713164,2.743205, 2.749723,2.753326]],
[[1.232765,1.341468,1.346242,1.346456,1.346603],[1.696334,1.895051,1.907863,1.909372,1.910032],[2.238556,2.722308,2.750227,2.757492,2.760331]],[[1.235641,1.344179,1.348399,1.348852,1.349151],
[1.701584,1.899687,1.912100,1.913716,1.914301],[2.246716,2.730136,2.757795,2.765451,2.767957]],[[1.238478,1.346586,1.351096,1.351554,1.351786],[1.705517,1.903961,1.916485,1.917941,1.918521],
[2.254955,2.737834,2.766094,2.772398,2.774493]],[[1.241981,1.349207,1.353765,1.353986,1.354179],[1.711073,1.908307, 1.920913,1.922321,1.923639],[ 2.263672,2.745290,2.772925,2.780656,2.783772]],
[[1.244816,1.352020,1.356149,1.356345,1.356684],[1.716266,1.912950,1.926254,1.927599,1.928130],[2.271981,2.753501,2.782581,2.788330,2.790409]],[[1.248790,1.354620,1.358915,1.359235,1.359487],
[1.720193,1.917516,1.931100,1.932476,1.933184],[2.280294,2.761655,2.789957,2.795878,2.797913]],[[1.252395,1.357197,1.361947,1.362265,1.362569],[1.725837,1.922129,1.936171,1.937710,1.938192],
[2.289980,2.769950,2.797966,2.804613,2.808125]],[[1.254989,1.360500,1.365236,1.365600,1.365772],[1.730801,1.927967,1.941870,1.943303,1.943724],[2.301392,2.777844,2.808374,2.813311,2.815859]],
[[1.258229,1.363725,1.368277,1.368505,1.368863],[1.736322,1.933243,1.947243,1.948883,1.949207],[2.310398,2.788013,2.816225,2.821322,2.824270]],[[1.262239,1.366679,1.371643,1.372149,1.372374],
[1.741964,1.939180,1.952111,1.953539,1.954109],[2.320316,2.796900,2.825448,2.831913,2.834508]],[[1.265966,1.370524,1.374604,1.374777,1.374852],[1.747394,1.945610,1.957768,1.959244,1.959426],
[2.329219,2.807824,2.835817,2.840917,2.843434]],[[1.269240,1.373742,1.377815,1.378134,1.378315],[1.753500,1.950991,1.963491,1.964871,1.965069],[2.339306,2.817101,2.846078,2.851046,2.853604]],
[[1.272641,1.376577,1.381153,1.381548,1.381751],[1.758805,1.957234,1.969616,1.970576,1.970974],[2.350632,2.828103,2.856533,2.860616,2.862433]],[[1.276040,1.380221,1.384750,1.385074,1.385378],
[1.765472,1.963210,1.974486,1.975609,1.975930],[2.362393,2.839029,2.866151,2.871058,2.872654]],[[1.279592,1.383943,1.388012,1.388368,1.388473],[1.772518,1.969691,1.980481,1.981478,1.981607],
[2.373682,2.849856,2.875078,2.878643,2.880942]],[[1.284859,1.387693,1.390913,1.391173,1.391456],[1.779402,1.975368,1.985724,1.987007,1.987355],[2.383993,2.859437,2.885461,2.889821,2.891402]],
[[1.289184,1.390484,1.394646,1.395081,1.395330],[1.788029,1.981488,1.992285,1.993116,1.993443],[2.396559,2.872062,2.896479,2.900194,2.901336]],[[1.293119,1.394082,1.398722,1.398860,1.399112],
[1.795241,1.987604,1.997985, 1.998916,1.999079],[2.409003,2.881957,2.906464,2.911006,2.912487]],[[1.297743,1.398436,1.402654,1.403074,1.403188],[1.802223,1.994707,2.004819,2.006503,2.006985],
[ 2.420513,2.894229,2.918304,2.920887,2.922340]],[[1.302930,1.402841,1.406762,1.407276,1.407490],[1.809158,2.001483,2.012325,2.013257,2.013485],[2.431878,2.905888,2.927930,2.931367,2.933102]],
[[1.307095,1.407343,1.411457,1.411639,1.411814],[1.816827,2.009906,2.019668,2.020551,2.020959],[2.442134,2.918303,2.939940,2.942373,2.943662]],[[1.311586,1.411896,1.415374,1.415582,1.415698],
[1.824857,2.018009,2.027714,2.028695, 2.029367],[2.455465,2.928892,2.951336,2.955135,2.955942]],[[1.317376,1.415971,1.419370,1.419562,1.419777],[1.833343,2.026390,2.035115,2.035950,2.036448],
[2.468240,2.941650,2.962908,2.965146,2.966898]],[[1.323352,1.420220,1.423625,1.423804,1.423819],[1.841864,2.034022,2.042662,2.044230,2.044388],[2.483054,2.955380,2.976538,2.979340, 2.980014]],
[[1.327864,1.425241,1.428682,1.428849,1.428957],[1.851771,2.042079,2.052127,2.053173,2.053381],[2.500143,2.968288,2.989569,2.992682,2.994808]],[[1.333507,1.430859, 1.433526,1.433628,1.433639],
[1.861211,2.052014,2.058428,2.059145,2.059377],[2.512775,2.983458,3.002817,3.007144,3.008677]],[[1.339192,1.435554,1.438252,1.438392,1.438405],[1.869505,2.058946,2.065555,2.066094,2.066162],
[2.527290,3.000886,3.017675,3.021515,3.022463]],[[1.344926,1.440531,1.442914,1.442974, 1.443076],[1.878587,2.066095,2.074048,2.074708,2.074738],[2.544121,3.014506,3.031640,3.033255,3.033942]],
[[1.350273,1.445440,1.448123,1.448228,1.448236],[1.886953,2.075533,2.082448,2.082848,2.082870],[2.559929,3.030217,3.045733,3.048442,3.049289]],[[1.356197,1.450892,1.453143,1.453262,1.453311],
[1.899271,2.085649,2.092158,2.092533,2.092569],[2.579905,3.045143,3.062709,3.065217,3.065598]],[[1.362590,1.456028,1.458898,1.458992, 1.459029],[1.909990,2.095269,2.101756,2.101928, 2.101952],
[2.599545,3.063189,3.081057,3.083965,3.085387]],[[1.369579,1.462782,1.465478,1.465547,1.465578],[1.920444,2.105089,2.111554,2.111633,2.111780],[2.622261,3.082044,3.099112,3.102763,3.103441]],
[[1.376458,1.469854,1.472349,1.472453,1.472531],[1.933259,2.114999,2.121438,2.121521,2.121702],[2.643282, 3.101187,3.118645,3.121287,3.121690]],[[1.384282,1.476966,1.480404,1.480559,1.480576],
[1.948464,2.126446,2.133841,2.133936,2.134194],[2.667571,3.121373,3.143574,3.145042,3.145468]],[[1.391945,1.485232,1.487427,1.487522,1.487593], [1.961357,2.138602,2.144182,2.144204,2.144253],
[2.689127, 3.146294,3.161756,3.163726,3.164096]],[[1.401008,1.492834,1.494943,1.495068,1.495171],[1.978307,2.151764,2.156617,2.156934,2.157615],[2.716153,3.169754,3.188961, 3.192023,3.193316]],
[[1.411670,1.500951,1.503252,1.503442,1.503842],[1.993349,2.165997,2.173272,2.173545,2.173771],[2.744999,3.198712,3.214096,3.216680,3.217122]],[[1.420877,1.510675,1.511970,1.511986,1.512084],
[2.010536,2.184522,2.190798,2.191140,2.191611],[2.769564,3.225578,3.237936,3.240050,3.240793]],[[1.433263,1.519837,1.521600,1.521629,1.521645],[2.031463, 2.201170,2.208535,2.208754,2.209073],
[2.799616,3.252830,3.274006,3.274860,3.276932]],[[1.443600,1.532697,1.533838,1.534082,1.534365],[2.052896,2.218109,2.224369,2.224911,2.225384],[2.842197,3.294485, 3.309482,3.310631,3.311442]],
[[1.455112,1.544403,1.545380,1.545547,1.545562],[2.073102,2.241080,2.246534,2.246710,2.247286],[2.882628,3.328245,3.339912,3.340922,3.341217]],[[1.471607,1.559106,1.560338,1.560338,1.560361],
[2.096964,2.264123,2.269144,2.269487,2.269782],[2.923850,3.368009,3.382317,3.383336,3.384313]],[[1.488860,1.575721,1.576582,1.576582,1.576732],[2.121520,2.290532,2.294708,2.295162,2.295703],
[2.971624,3.413466,3.423838,3.424968,3.425139]],[[1.507300,1.596956,1.597971,1.597971,1.597971],[2.151915,2.320520,2.325255,2.325522,2.325522],[3.029458,3.460251,3.473393,3.474227,3.474227]],
[[1.531756,1.618122,1.618397,1.618397,1.618397],[2.200397,2.355904,2.358182,2.359174,2.359174],[3.087042,3.515429,3.529342,3.529363,3.529363]],[[1.560381,1.649405,1.649405,1.649405,1.649405],
[2.247114,2.408982,2.411867,2.411867,2.411867],[3.172736,3.606243, 3.619386, 3.620481,3.620959]],[[1.604588,1.685943,1.685943,1.685943,1.685943],[2.308004,2.464537,2.465797,2.465797,2.465797],
[3.289425,3.718164,3.734348,3.736920,3.736979]],[[1.673977,1.745509, 1.745509,1.745509,1.745509],[2.434576,2.568862,2.570255,2.570255,2.570255],[3.454727,3.935357,3.941029,3.941029,3.941029]]]

// 1. GET  CRITICAL VALUE FROM THE CRITICAL VALUE TABLE. 
var tblindex = criticalVTable.slice(alphaIndex);
var tblind = tblindex[0];
var htablindex = tblind[h];
var criticalValue = htablindex[period];


///////////////////////////
//MONITORING  PERIOD
//////////////////////////

// Filter to the area of interest, mask clouds, add variables.
var MonitorLandsat = moniCollection
  //.map(addNDVI)
  .map(addDependents)
  .map(addHarmonics(harmonicFrequencies));
  
// Do the prediction.

var predictedValues = MonitorLandsat.map(function(image) {
  return image.addBands(
    image.select(independents)
      .multiply(harmonicTrendCoefficients)
      .reduce('sum')
      .rename('predicted'));
});

var mcomputResiduals = function(image){
  return image.addBands(
    image.select('NDMI')
    .subtract(image.select('predicted'))
    .reduce('sum')
    .rename('residual'));

}

// Calculate the MOSUM for monitoring period

///...calculate the Cumsum of residuals and do a shift trick to get MOSUM
var mresiduals = predictedValues.map(mcomputResiduals);
var mresid = mresiduals.select('residual')
var moResiduals = mresid.toArray();
var histResidual =  res.toArray();
var moMoresid = ee.ImageCollection([histResidual,moResiduals]);
var MonResid = moMoresid.toArrayPerBand(imageAxis);
var moRes = MonResid.toArray()

var cumSumRes = moRes.arrayAccum(imageAxis, 'sum'); 
var imZer0 = ee.Image(0).int();
var ksizeplu = ksize.add(1)
var dubImZer0 = imZer0.arrayRepeat(imageAxis,ksize.int());
var fx10 = dubImZer0.arrayRepeat(bandAxis, 1);
var verc= fx10.toArray().float();
var shiftedCumR = ee.ImageCollection([verc, cumSumRes]);
var ShfCumSumRight1 = shiftedCumR.toArrayPerBand(imageAxis);
var ShfCumSumRight = ShfCumSumRight1.toArray()
var shiftedCumL = ee.ImageCollection([cumSumRes, verc]);
var ShfCumSumLeft = shiftedCumL.toArrayPerBand(imageAxis);
var ShiftDifference = ShfCumSumLeft.subtract(ShfCumSumRight);
var msk = ShiftDifference.arrayMask(ShfCumSumLeft);
var ter = ee.Image(0).float().toArray();
var tesser = ter.arrayRepeat(bandAxis, 1);
var tessermaske = ee.ImageCollection([tesser, msk]);
var tessermsk = tessermaske.toArrayPerBand(imageAxis);

///... now slash off unwanted first part of Cumsum of residuals to get MOSUM residuals
var cotx2 = tessermsk.arrayLength(imageAxis);
var cumKsizeDiff = cotx2.subtract(ksize).int();
var knsiz1 = imZer0.arrayRepeat(imageAxis, ksize.int());
var knsiz2 = knsiz1.arrayRepeat(bandAxis, 1);
var knsiz =knsiz2.toArray();


var dcumn = ee.Image(1).multiply(cumKsizeDiff);
var imOnes = ee.Image(1).int();
var knm1 = imOnes.arrayRepeat(imageAxis, dcumn);
var knm2 = knm1.arrayRepeat(bandAxis, 1);
var knm = knm2.toArray();

var shmasker = ee.ImageCollection([knsiz, knm]);
var Shfmask = shmasker.toArrayPerBand(imageAxis);
var moResmosum = tessermsk.arrayMask(Shfmask);

var dumoSigmaStanda = sigmaStanda.arrayRepeat(bandAxis,1);
var mobleng = moResmosum.arrayLength(imageAxis);
var moSigmaStanda = dumoSigmaStanda.arrayRepeat(imageAxis,mobleng).toArray();

///... now standardise the mosum
var moResmosumStanda1 = moResmosum.divide(moSigmaStanda);
var moResmosumStanda = moResmosumStanda1.abs();
///...Time stamp...needed for time of change later

// Number of milliseconds in a day.
var MSEC_PER_DAY = 365.25 * 24 * 3600 * 1000;

function makeDate(image) {
  // Create a band containing the time of each image,
  // converted to days since the epoch (1/1/1970).
  var md = image.metadata("system:time_start").divide(MSEC_PER_DAY).add(1970);
  // Copy the mask from band 0.
  var mask = image.select(0).mask();
  return md.mask(mask);
}

var vsr = MonitorLandsat.map(makeDate);


var time = vsr.toArray();

///...magnitude ..needed for magnitude change later

var magnTudmas = moResiduals;



/// 3.Calculate the critical borders

///.. generating a sequence of values  from N obs in historical sample to N obs in monitoring sample with 
 //interval of 1 interval and this is not a straight forward procedure
 
 
var monisize1 = moResiduals.arrayLength(imageAxis);
var signb = ee.Image(1).toArray()
var kplushistzie = signb.add(histtimeseiesLength);

//logPlus
var fo1 =  kplushistzie.divide(histtimeseiesLength);
var foLogplus1 = fo1.log().max(1);

//border
var foLop1 = foLogplus1.multiply(2);
var foLopSq1 = foLop1.sqrt();
var foLopSqx = foLopSq1.arrayRepeat(imageAxis,monisize1.int());
var criticalBorder1 = foLopSqx.multiply(ee.Image(criticalValue));
var criticalBorder = criticalBorder1.arrayRepeat(bandAxis, 1);
 

///...now create masks to mask out parts that belong to the history period

var yeOnes = ee.Image(1).int();
//var timeMa = time.arrayLength(imageAxis)
var moRepeat1 = yeOnes.arrayRepeat(imageAxis, monisize1.int()).toArray();
var moRepeat  = moRepeat1.arrayRepeat(bandAxis, 1);
var rMa = moResmosumStanda.arrayLength(imageAxis);
var rt = rMa.subtract(monisize1);
var rtx = ee.Image(0).int();
var rhRepeat1 = rtx.arrayRepeat(imageAxis, rt.int()).toArray();
var rhRepeat  = rhRepeat1.arrayRepeat(bandAxis, 1);

var rZOnesx = ee.ImageCollection([rhRepeat,moRepeat]);
var rMZOnesxo = rZOnesx.toArrayPerBand(0);
var rMZexOnesx = rMZOnesxo.arrayRepeat(bandAxis, 1);

/// ...now start masking 
var xmoResmosumStanda = moResmosumStanda.arrayMask(rMZexOnesx);
var  chiefMasker = xmoResmosumStanda.divide(criticalBorder).int();

/// 4. Check if a break is detected

///... Now mask out the positions where critical boundary is not closed
// you do this for change, time and magnitude of change image
var change = xmoResmosumStanda.arrayMask(chiefMasker);
var timechange = time.arrayMask(chiefMasker);
var magniofchange = magnTudmas.arrayMask(chiefMasker);

///...mask out the positions where the magnitude of change is positive

//var imangx = magniofchange.divide(magniofchange.abs());
//var msMang = imangx.max(0);
//var msChangex = msMang.subtract(1);
//var fichange = change.arrayMask(msChangex);
//var timxchange = timechange.arrayMask(msChangex);
//var tMagofchange = magniofchange.arrayMask(msChangex);

///...Use the magnitude threshold to mask out positions where the change magnitude is
//smaller then the magnitude threshold
//var tmagofange = tMagofchange.divide(magnitudeThreshold).int();
//var vichange =fichange.arrayMask(tmagofange);
//var ctimevichange =timxchange.arrayMask(tmagofange);
//var xtmagofchange =tMagofchange.arrayMask(tmagofange);


///... now slice out the first position the change has occured
var firstChange = change.arraySlice(0,0,1).toArray();
var timeOfchange = timechange.arraySlice(0,0,1).toArray();
var magnitude = magniofchange.arraySlice(0,0,1).toArray();
var criticalBorder1 = criticalBorder.arraySlice(0,0,1).toArray();

///...mask out locations where no change is detected;
var time0x = ee.Image(0).int().toArray();
var timZx = time0x.arrayRepeat(bandAxis,1);
var mas = ee.ImageCollection([magnitude,timZx]);
var timOChx = ee.ImageCollection([timeOfchange,timZx]);
var firstChx = ee.ImageCollection([firstChange,timZx]);

var timeOfCvx = timOChx.toArrayPerBand(0);
var tMa = mas.toArrayPerBand(0);
var firsttCha = firstChx.toArrayPerBand(0);

var cox = tMa.arrayReduce('sum', [imageAxis]);
var Timecox = timeOfCvx.arrayReduce("max", [imageAxis]);
var firsttChab = firsttCha.arrayReduce('max', [imageAxis]);


var Cnk1 = cox.mask(cox.arrayGet([0,0]).neq(0)).toArray();
var timeCnk1 = Timecox.mask(Timecox.arrayGet([0,0]).neq(0)).toArray();
var firtsChnk = firsttChab.mask(firsttChab.arrayGet([0,0]).neq(0)).toArray();
var firstCha = firtsChnk.subtract(criticalBorder1);

///... flatten the array to allow for use of palletes
var timeCnk2 = timeCnk1.arrayFlatten([['x'],['y']]);
var Cnk = Cnk1.arrayFlatten([['x'],['y']]);
var firstChaDif = firstCha.arrayFlatten([['x'],['y']]);

timeCnk2 = timeCnk2.select(['.*'],['breakTime']).set({'bfast:label' : 'Time of change', 'bfast:result' : 'breakTime'})
Cnk = Cnk.select(['.*'],['breakMagnitude']).set({'bfast:label' : 'Magnitude of change', 'bfast:result' : 'breakMagnitude'})
// Compute temporal metrics for prediction
var b1_median_mon = MonitorLandsat.select('B1')
    .median();
var b1_median_his = histCollection.select('B1')
    .median();
    
var b2_median_mon = MonitorLandsat.select('B2')
    .median();
var b2_median_his = histCollection.select('B2')
    .median();   
    
var b3_median_mon = MonitorLandsat.select('B3')
    .median();
var b3_median_his = histCollection.select('B3')
    .median(); 
    
var b4_median_mon = MonitorLandsat.select('B4')
    .median();
var b4_median_his = histCollection.select('B4')
    .median(); 

var b5_median_mon = MonitorLandsat.select('B5')
    .median();
var b5_median_his = histCollection.select('B5')
    .median(); 
    
var b7_median_mon = MonitorLandsat.select('B7')
    .median();
var b7_median_his = histCollection.select('B7')
    .median();  
 var ndmi_median_mon = MonitorLandsat.select('NDMI')
    .median();
var ndmi_median_his = histCollection.select('NDMI')
    .median();  
var monrest = mresiduals.select('residual')
    .median()
    .float(); 
var hisrest = residuals.select('residual')
    .median()
    .float(); 
// Compute temporal metrics for prediction
var b1_ratio_median = b1_median_his.divide(b1_median_mon);
var b2_ratio_median = b2_median_his.divide(b2_median_mon);
var b3_ratio_median = b3_median_his.divide(b3_median_mon);
var b4_ratio_median = b4_median_his.divide(b4_median_mon);
var b5_ratio_median = b5_median_his.divide(b5_median_mon);
var b7_ratio_median = b7_median_his.divide(b7_median_mon);
var ndmires_ratio_median = hisrest.divide(monrest);
var ndmi_ratio_median = ndmi_median_his.divide(ndmi_median_mon);



//var ndmi_perc1 = ndmi_perc.mask(Timecox.arrayGet([0,0]));

var band_metrics = ee.Image(b1_ratio_median)
.addBands(b2_ratio_median)
.addBands(b3_ratio_median)
.addBands(b4_ratio_median)
.addBands(b5_ratio_median)
.addBands(b7_ratio_median)
.addBands(ndmi_ratio_median)
.addBands(ndmires_ratio_median)
.addBands(firstChaDif.toFloat())
.addBands(Cnk.toFloat());

var band_metrics_masked = band_metrics.mask(timeCnk2);


// historyStart,monitoringEnd,historyEnd,monitoringStart,h,period,alpha,magnitudeThreshold,harmonics

var bfastResults = ee.ImageCollection([moRes,timeCnk2,Cnk,criticalBorder1,band_metrics_masked]).set({
  'bfast:historyStart' : historyStart,
  'bfast:historyEnd' : historyEnd,
  'bfast:monitoringStart' : monitoringStart,
  'bfast:monitoringEnd' : monitoringEnd,
  'bfast:h' : h,
  'bfast:period' : period,
  'bfast:alpha' : alpha,
  'bfast:magnitudeThreshold' : magnitudeThreshold,
  'bfast:harmonics' : harmonics
})

var results = {
  'bfastResults' : bfastResults,
  'residuals' : residuals,
  'predictedValues' :predictedValues,
  'mresiduals' : mresiduals
}
return results
}
exports.bfastMonitor = bfastMonitor;

// testing
// var test = bfastMonitor()
// print(test.bfastResults)

//Map.addLayer(test.bfastResults)
// // Export pixel time series
/// Export.table(LakeICE, 'LakeICE_peruGee', {fileFormat: 'CSV'});
/// Export.table(LakeDAM, 'LakeDAM_peruGee', {fileFormat: 'CSV'});


//Map.centerObject(roi, 11);
//Map.addLayer(roi, {}, 'ROI');
//Map.addLayer(region, {}, 'Image2');
//Map.addLayer(moRes, {}, 'Image');
//Map.addLayer(timeCnk2, {min: 2017, max: 2018,'palette' : '00BFFF,CC2EFA,A901DB,6A0888,5858FA,0101DF,2E2EFE,0B0B61'},"Time of change");
//Map.addLayer(Cnk, {min: -0.4, max:0,'palette' : 'F4FA58,FFFF00,F7D358,F7D358,FFBF00,FF4000,B43104,8A0808'}, 'Magnitude of change');
//Map.addLayer(criticalBorder1 , {}, 'Time of change');

