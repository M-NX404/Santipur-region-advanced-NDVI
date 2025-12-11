/*******************************************************
Fixed NDVI script — robust Sentinel-2 masking (no QA60 crash)
- Handles missing QA60 / SCL bands safely
- Builds yearly median NDVI, computes slope, plots chart
*******************************************************/

// ---------- USER PARAMETERS ----------
var sensor = 'S2';             // 'S2' or 'L8'
var yearStart = 2018;
var yearEnd = 2023;
var region = ee.Geometry.Rectangle([88.35, 23.05, 88.55, 23.35]);
var studyScale = (sensor === 'S2') ? 10 : 30;
Map.centerObject(region, 12);

// ---------- Robust S2 mask ----------
function maskS2SR(image) {
  // Band names available on this image
  var bnames = image.bandNames();

  // If SCL exists, use it (preferred)
  var sclExists = bnames.contains('SCL');
  var sclMask = ee.Image(ee.Algorithms.If(
    sclExists,
    // keep main land classes: Vegetation(4), Not vegetated(5), Bare soils(6), Water(7)
    image.select('SCL').eq(4).or(image.select('SCL').eq(5))
                        .or(image.select('SCL').eq(6)).or(image.select('SCL').eq(7)),
    // else return a constant image of 1 (no mask from SCL)
    ee.Image.constant(1)
  ));

  // If QA60 exists, use it (0 = no clouds)
  var qaExists = bnames.contains('QA60');
  var qaMask = ee.Image(ee.Algorithms.If(
    qaExists,
    image.select('QA60').eq(0),
    // else, if MSK_CLDPRB exists, use MSK_CLDPRB < threshold (e.g. 50)
    ee.Image(ee.Algorithms.If(
      bnames.contains('MSK_CLDPRB'),
      image.select('MSK_CLDPRB').lt(50),
      // fallback: keep all pixels
      ee.Image.constant(1)
    ))
  ));

  // Combine masks (SCL mask OR qaMask) to be permissive (if either says ok, keep)
  var combined = sclMask.or(qaMask);

  // Update mask and return; keep original bands
  return image.updateMask(combined);
}

// ---------- Robust L8 mask (unchanged) ----------
function maskL8SR(image) {
  var qa = image.select('QA_PIXEL');
  var cloud = 1 << 3;
  var cloudShadow = 1 << 4;
  var mask = qa.bitwiseAnd(cloud).eq(0).and(qa.bitwiseAnd(cloudShadow).eq(0));
  return image.updateMask(mask);
}

// ---------- NDVI function ----------
function addNDVIBands(img) {
  if (sensor === 'S2') {
    return img.addBands(img.normalizedDifference(['B8', 'B4']).rename('NDVI'));
  } else {
    return img.addBands(img.normalizedDifference(['SR_B5', 'SR_B4']).rename('NDVI'));
  }
}

// ---------- Helper to load collection for a year ----------
function collectionForYear(y) {
  var start = ee.Date.fromYMD(y,1,1);
  var end = start.advance(1,'year');
  if (sensor === 'S2') {
    return ee.ImageCollection('COPERNICUS/S2_SR')
             .filterDate(start, end)
             .filterBounds(region)
             .map(maskS2SR)
             .map(addNDVIBands);
  } else {
    return ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
             .filterDate(start, end)
             .filterBounds(region)
             .map(maskL8SR)
             .map(addNDVIBands);
  }
}

// ---------- Build yearly images safely using server-side map ----------
var years = ee.List.sequence(yearStart, yearEnd);

// Print image counts per year (debug)
var counts = years.map(function(y){
  y = ee.Number(y);
  var c = collectionForYear(y).size();
  return c;
});
print('Image counts per year:', counts);

// Create yearly NDVI images (with nodata -9999 if empty)
var yearlyList = years.map(function(y){
  y = ee.Number(y);
  var col = collectionForYear(y);
  var med = ee.Image(ee.Algorithms.If(
    col.size().gt(0),
    col.select('NDVI').median().clip(region),
    ee.Image.constant(-9999).clip(region).rename('NDVI')
  ));
  med = med.addBands(ee.Image.constant(y).rename('t')).set('year', y);
  return med;
});

var yearCollection = ee.ImageCollection.fromImages(yearlyList);
print('Yearly NDVI collection size:', yearCollection.size());

// Visual check of last year NDVI
var ndviVis = {min:-0.3, max:0.9, palette:['ffffff','ce7e45','df923d','f1b555','fcd163','99b718','74a901','66a000','529400','3e8601']};
var lastYear = yearCollection.filter(ee.Filter.eq('year', yearEnd)).first();
Map.addLayer(lastYear.select('NDVI'), ndviVis, 'NDVI ' + yearEnd);

// ---------- Trend: linear slope per year ----------
var paired = yearCollection.select(['t','NDVI']);
var fit = paired.reduce(ee.Reducer.linearFit());
var slope = fit.select('scale').rename('NDVI_slope_per_year');

// Mask nodata-only pixels
var anyValid = yearCollection.select('NDVI').reduce(ee.Reducer.min()).neq(-9999);
slope = slope.updateMask(anyValid);

Map.addLayer(slope, {min:-0.05, max:0.05, palette:['red','white','green']}, 'NDVI slope (per year)');

// ---------- Regional mean chart ----------
var meanPerYearFC = yearCollection.map(function(img){
  var yr = ee.Number(img.get('year'));
  var mean = img.select('NDVI').reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: region,
    scale: studyScale,
    bestEffort: true,
    maxPixels: 1e13
  });
  return ee.Feature(null, {'year': yr, 'meanNDVI': mean.get('NDVI')});
});
print('Mean NDVI per year (FeatureCollection):', meanPerYearFC);

var chart = ui.Chart.feature.byFeature(meanPerYearFC, 'year', ['meanNDVI'])
  .setOptions({title: 'Mean NDVI per year', vAxis:{title:'NDVI'}, hAxis:{title:'Year'}, lineWidth:2, pointSize:4});
print(chart);

// ---------- Area stats ----------
var pixelArea = ee.Image.pixelArea();
print('Area pos trend (m^2):', slope.gt(0).multiply(pixelArea).reduceRegion({
  reducer: ee.Reducer.sum(), geometry: region, scale: studyScale, maxPixels:1e13, bestEffort:true
}));
print('Area neg trend (m^2):', slope.lt(0).multiply(pixelArea).reduceRegion({
  reducer: ee.Reducer.sum(), geometry: region, scale: studyScale, maxPixels:1e13, bestEffort:true
}));

print('Finished. If another error appears, paste the first line(s) of the Console error here.');
