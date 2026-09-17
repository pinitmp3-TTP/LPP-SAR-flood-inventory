/**********************************************************************
 * LPP SAR — A1 : InSAR coherence (บน Watershed_v2)  -> Table 9
 * ------------------------------------------------------------------
 *  ต่อ 1 event: cohLoss = ref - co ; disturbance ในเขตเมือง ; cross-check กับ amplitude
 *  แก้รีวิวข้อ 6: รายงาน %confirmed "เฉพาะ built-up" (ไม่มีเงื่อนไข backscatter = independent จริง)
 *  พื้นที่นับบน Watershed_v2 ; Otsu (amplitude) บน aoi
 *
 *  coherence assets (Path 164, HyP3 GAMMA) ที่มีอยู่แล้วใน project:
 *    Nangka ref = nangka_corr_ref2608 (26Sep-08Oct,12d) · co = nangka_corr_co0814 (08-14Oct,6d)
 *    Dianmu ref = corr_ref_0921       (09-21Sep,12d)    · co = corr_ref_2127       (21-27Sep,6d)
 *    Noru   ref = noru_corr_ref1628    (16-28Sep,12d)    · co = noru_corr_co2810    (28Sep-10Oct,12d)
 **********************************************************************/

// ============ ตั้งค่า ============
var aoi = ee.Geometry.Rectangle([101.28,14.05,102.42,14.95]);
var WS  = ee.FeatureCollection('projects/lamphra-pleong-project/assets/Watershed_v2');
var ws  = WS.geometry();
var PASS='DESCENDING', POLAR='VV';
var PROJ='projects/lamphra-pleong-project/assets/';     // path นำหน้า asset
var LULC = ee.ImageCollection('ESA/WorldCover/v200').first();
var builtup = LULC.eq(50);                              // class 50 = เขตเมือง/สิ่งปลูกสร้าง
var LOSS = 0.3;                                         // เกณฑ์ coherence loss หลัก (Table 9)
var SCALE = 20;                                         // A1 นับพื้นที่ที่ 20 m (coherence หยาบกว่า 10 m)

// event config: coherence assets + วันภาพ amplitude (ช่วงพายุ/peak เหมือน Table 4)
var EV = {
  Nangka2020: {ref:'nangka_corr_ref2608', co:'nangka_corr_co0814',
               a0:'2020-10-13', a1:'2020-10-19'},
  Dianmu2021: {ref:'corr_ref_0921',        co:'corr_ref_2127',
               a0:'2021-09-24', a1:'2021-10-02'},
  Noru2022:   {ref:'noru_corr_ref1628',    co:'noru_corr_co2810',
               a0:'2022-09-28', a1:'2022-10-16'}
};

// ============ ฟังก์ชัน ============
function otsu(h){h=ee.Dictionary(h);var c=ee.Array(h.get('histogram')),m=ee.Array(h.get('bucketMeans'));
var s=m.length().get([0]),t=c.reduce(ee.Reducer.sum(),[0]).get([0]),su=m.multiply(c).reduce(ee.Reducer.sum(),[0]).get([0]),mn=su.divide(t);
var idx=ee.List.sequence(1,s);var bss=idx.map(function(i){var ac=c.slice(0,0,i);var aC=ac.reduce(ee.Reducer.sum(),[0]).get([0]);
var aM=m.slice(0,0,i).multiply(ac).reduce(ee.Reducer.sum(),[0]).get([0]).divide(aC);var bC=t.subtract(aC);var bM=su.subtract(aC.multiply(aM)).divide(bC);
return aC.multiply(aM.subtract(mn).pow(2)).add(bC.multiply(bM.subtract(mn).pow(2)));});return m.sort(ee.Array(bss)).get([-1]);}

var hand=ee.Image('MERIT/Hydro/v1_0_1').select('hnd');
var slope=ee.Algorithms.Terrain(ee.Image('WWF/HydroSHEDS/03VFDEM')).select('slope');
var swater=ee.Image('JRC/GSW1_4/GlobalSurfaceWater').select('seasonality');

// amplitude flood (C2) + after image + Otsu threshold ของ event
function amplitude(a0,a1){
  var coll=ee.ImageCollection('COPERNICUS/S1_GRD').filter(ee.Filter.eq('instrumentMode','IW'))
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation',POLAR))
    .filter(ee.Filter.eq('orbitProperties_pass',PASS)).filter(ee.Filter.eq('resolution_meters',10))
    .select(POLAR).filterBounds(aoi);
  var after=coll.filterDate(a0,a1).mosaic().clip(aoi).focal_mean(50,'circle','meters');
  var low=after.updateMask(hand.lt(10)).updateMask(slope.lt(3));
  var thr=ee.Number(otsu(low.reduceRegion({reducer:ee.Reducer.histogram(255),
    geometry:aoi,scale:30,maxPixels:1e13,bestEffort:true}).get(POLAR)));   // Otsu บน aoi
  var fl=after.lt(thr).selfMask().where(swater.gte(10),0).selfMask();
  fl=fl.updateMask(fl.connectedPixelCount(25).gte(8)).updateMask(slope.lt(5)).updateMask(hand.lt(15));
  return {flood:fl.unmask(0), after:after, thr:thr};
}

// พื้นที่ (ha) ของ mask บน ws
function areaHa(mask){
  return ee.Number(ee.Image.pixelArea().updateMask(mask).reduceRegion({
    reducer:ee.Reducer.sum(), geometry:ws, scale:SCALE, maxPixels:1e13, tileScale:8
  }).get('area')).divide(10000);
}

// ============ ประมวลผลราย event ============
function runA1(name){
  var e=EV[name];
  var refC=ee.Image(PROJ+e.ref).select(0);           // coherence 0-1 (band แรก)
  var coC =ee.Image(PROJ+e.co).select(0);
  var cohLoss=refC.subtract(coC);                    // ค่าบวก = coherence ตก

  var A=amplitude(e.a0,e.a1);
  var ampFlood=A.flood.eq(1);                        // amplitude flood (C2)
  var lowBack =A.after.lt(A.thr);                    // เกณฑ์ backscatter ต่ำ (คล้ายน้ำ)

  // coherence-flagged flood:
  //   ในเมือง = cohLoss>0.3 (ไม่มีเงื่อนไข backscatter)  <- โซน independent
  //   นอกเมือง = cohLoss>0.3 AND backscatter ต่ำ (กัน decorrelation ในนา)
  var cohBU   = builtup.and(cohLoss.gt(LOSS));
  var cohNonBU= builtup.not().and(cohLoss.gt(LOSS)).and(lowBack);
  var cohFlood= cohBU.or(cohNonBU);

  // Table 9
  var ampUrban = areaHa(ampFlood.and(builtup));                  // amplitude urban flood (ha)
  var urbanDist= areaHa(cohBU);                                  // urban coherence disturbance @0.3 (ha)
  var confAll  = areaHa(cohFlood.and(ampFlood)).divide(areaHa(cohFlood)).multiply(100);      // %conf ทั้งหมด
  var confBU   = areaHa(cohBU.and(ampFlood)).divide(areaHa(cohBU)).multiply(100);            // %conf เฉพาะเมือง (independent)

  print('==== '+name+' : A1 (Watershed_v2) ====');
  print('  Amplitude urban flood (ha)          =', ampUrban);
  print('  Urban coherence disturbance @0.3 (ha)=', urbanDist);
  print('  Coherence confirmed by amplitude — ALL (%)        =', confAll);
  print('  Coherence confirmed by amplitude — BUILT-UP only (%)=', confBU);  // <- ใช้ในบทความ (independent)

  // (ทางเลือก Figure 8a) disturbance ที่เกณฑ์ 0.2 / 0.4 — uncomment ถ้าจะ regen Figure 8
  // print('  Urban disturbance @0.2 / @0.4 (ha) =', areaHa(builtup.and(cohLoss.gt(0.2))), areaHa(builtup.and(cohLoss.gt(0.4))));

  Map.addLayer(cohLoss.updateMask(cohLoss.gt(LOSS)).clip(ws),
    {min:0.3,max:0.8,palette:['#fee08b','#f46d43','#a50026']}, name+' cohLoss>0.3', false);
  return null;
}

Map.centerObject(ws,10);
Object.keys(EV).forEach(runA1);

/**********************************************************************
 * ถ้า error "band ... not found" ที่ .select(0):
 *   coherence asset อาจมีชื่อ band เฉพาะ — พิมพ์ ee.Image(PROJ+'nangka_corr_ref2608').bandNames()
 *   ดูชื่อจริง แล้วเปลี่ยน .select(0) เป็น .select('<ชื่อ band>')
 * ถ้า "capacity exceeded": เพิ่ม SCALE เป็น 30 (บรรทัดบนสุด)
 **********************************************************************/
