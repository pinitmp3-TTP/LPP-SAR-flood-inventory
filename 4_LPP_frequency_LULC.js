/**********************************************************************
 * LPP SAR — FREQUENCY x LULC cross-tabulation (บน Watershed_v2)
 * เวอร์ชันเบา: reduce 3 รอบ + tileScale 16 (กัน capacity) + print อ่านง่าย
 *   A) ท่วม >=1 ครั้ง แยกตาม LULC   -> Table 6 ซ้าย
 *   B) recurrent (3 ครั้ง) แยก LULC -> Table 6 ขวา
 *   C) พื้นที่แยก 1/2/3 ครั้ง        -> Table 5
 * ทุกพื้นที่บน Watershed_v2 ; Otsu บน aoi
 **********************************************************************/

var aoi = ee.Geometry.Rectangle([101.28, 14.05, 102.42, 14.95]);
var WS  = ee.FeatureCollection('projects/lamphra-pleong-project/assets/Watershed_v2');
var ws  = WS.geometry();
var PASS='DESCENDING', POLAR='VV';
var LULC = ee.ImageCollection('ESA/WorldCover/v200').first().rename('lc'); // 2021

function otsu(h){h=ee.Dictionary(h);var c=ee.Array(h.get('histogram')),m=ee.Array(h.get('bucketMeans'));
var s=m.length().get([0]),t=c.reduce(ee.Reducer.sum(),[0]).get([0]),su=m.multiply(c).reduce(ee.Reducer.sum(),[0]).get([0]),mn=su.divide(t);
var idx=ee.List.sequence(1,s);var bss=idx.map(function(i){var ac=c.slice(0,0,i);var aC=ac.reduce(ee.Reducer.sum(),[0]).get([0]);
var aM=m.slice(0,0,i).multiply(ac).reduce(ee.Reducer.sum(),[0]).get([0]).divide(aC);var bC=t.subtract(aC);var bM=su.subtract(aC.multiply(aM)).divide(bC);
return aC.multiply(aM.subtract(mn).pow(2)).add(bC.multiply(bM.subtract(mn).pow(2)));});return m.sort(ee.Array(bss)).get([-1]);}

var hand=ee.Image('MERIT/Hydro/v1_0_1').select('hnd');
var slope=ee.Algorithms.Terrain(ee.Image('WWF/HydroSHEDS/03VFDEM')).select('slope');
var swater=ee.Image('JRC/GSW1_4/GlobalSurfaceWater').select('seasonality');

function detectFlood(name,b0,b1,a0,a1){
  var coll=ee.ImageCollection('COPERNICUS/S1_GRD').filter(ee.Filter.eq('instrumentMode','IW'))
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation',POLAR))
    .filter(ee.Filter.eq('orbitProperties_pass',PASS)).filter(ee.Filter.eq('resolution_meters',10))
    .select(POLAR).filterBounds(aoi);
  var after_f=coll.filterDate(a0,a1).mosaic().clip(aoi).focal_mean(50,'circle','meters');
  var low=after_f.updateMask(hand.lt(10)).updateMask(slope.lt(3));
  var thr=ee.Number(otsu(low.reduceRegion({reducer:ee.Reducer.histogram(255),
    geometry:aoi,scale:30,maxPixels:1e13,bestEffort:true}).get(POLAR)));
  print('Otsu '+name+' (dB) =',thr);
  var fl=after_f.lt(thr).selfMask().where(swater.gte(10),0).selfMask();
  fl=fl.updateMask(fl.connectedPixelCount(25).gte(8)).updateMask(slope.lt(5)).updateMask(hand.lt(15));
  return fl.unmask(0).rename('flood');
}

var freq = detectFlood('Nangka2020','2020-09-20','2020-10-06','2020-10-13','2020-10-19')
  .add(detectFlood('Dianmu2021','2021-09-06','2021-09-20','2021-09-24','2021-10-02'))
  .add(detectFlood('Noru2022','2022-09-01','2022-09-24','2022-09-28','2022-10-16'));

Map.centerObject(ws,10);
Map.addLayer(freq.selfMask(),{min:1,max:3,palette:['#fff7bc','#fec44f','#d95f0e']},'Frequency');

// ---- print กลุ่มเป็นข้อความแบน ๆ (ha, ทศนิยม 1 ตำแหน่ง) ----
function show(label, groups, field){
  var lines=ee.List(groups).map(function(o){o=ee.Dictionary(o);
    return ee.String(field+' ').cat(ee.Number(o.get(field)).format('%d'))
      .cat(' = ').cat(ee.Number(o.get('sum')).divide(10000).format('%.1f')).cat(' ha');});
  print(label, lines);
}

// A) ท่วม >=1 : พื้นที่แยกตาม LULC  (lc = รหัส WorldCover)
var gA=ee.Image.pixelArea().addBands(LULC).updateMask(freq.gte(1)).reduceRegion({
  reducer:ee.Reducer.sum().group({groupField:1,groupName:'lc'}),
  geometry:ws,scale:10,maxPixels:1e13,tileScale:16});
show('A) FLOODED >=1  [Table 6 ซ้าย]', gA.get('groups'), 'lc');

// B) recurrent (3 ครั้ง) : พื้นที่แยกตาม LULC
var gB=ee.Image.pixelArea().addBands(LULC).updateMask(freq.eq(3)).reduceRegion({
  reducer:ee.Reducer.sum().group({groupField:1,groupName:'lc'}),
  geometry:ws,scale:10,maxPixels:1e13,tileScale:16});
show('B) RECURRENT 3x [Table 6 ขวา]', gB.get('groups'), 'lc');

// C) Table 5 : พื้นที่แยก 1/2/3 ครั้ง
var gC=ee.Image.pixelArea().addBands(freq.rename('f')).updateMask(freq.gte(1)).reduceRegion({
  reducer:ee.Reducer.sum().group({groupField:1,groupName:'f'}),
  geometry:ws,scale:10,maxPixels:1e13,tileScale:16});
show('C) FREQUENCY   [Table 5]', gC.get('groups'), 'f');

/**********************************************************************
 * รหัส LULC: 10=Tree 20=Shrub 30=Grassland 40=Cropland 50=Built-up
 *            60=Bare 80=Water 90=Wetland 95=Mangrove 100=Moss
 * ผลจะ print เป็น "lc 40 = 1807.0 ha" ฯลฯ — copy ทั้งชุดส่งผมได้เลย
 **********************************************************************/

// ============ Export Figure 4 : แผนที่ frequency (ไปทำรูปใน QGIS/ArcGIS) ============
Export.image.toDrive({image: freq.unmask(0).toByte(), description:'LPP_frequency_3events',
  folder:'GEE_LPP_flood', region: aoi, scale:10, maxPixels:1e13});
