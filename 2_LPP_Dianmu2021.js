/**********************************************************************
 * LPP SAR FLOOD — Dianmu2021  : ทุกอย่างในไฟล์เดียว
 * ------------------------------------------------------------------
 *  ส่วนที่ 3 (DETECTION)  -> Otsu + พื้นที่น้ำท่วม   = Table 4
 *  ส่วนที่ 4 (VALIDATION) -> C0/C0.5/C1/C2 + CI     = Table 7/8 + Figure 6
 *
 *  หลักการ:  Otsu คิดบนฉากกว้าง (aoi)  → threshold ดี
 *            พื้นที่ + การสุ่ม validation คิดบนลุ่มน้ำ (Watershed_v2)
 **********************************************************************/

// ============ ส่วนที่ 1 : ตั้งค่า (เฉพาะพายุนี้) ============
var EVENT = 'Dianmu2021';
var aoi = ee.Geometry.Rectangle([101.28, 14.05, 102.42, 14.95]);  // ฉากกว้าง (โหลดภาพ + Otsu)
var WS  = ee.FeatureCollection('projects/lamphra-pleong-project/assets/Watershed_v2');
var ws  = WS.geometry();                                          // ลุ่มน้ำ (นับพื้นที่ + สุ่ม)

// วันที่ (before = ก่อนพายุ , afterDet = ช่วงพายุ/peak)
var before_start='2021-09-06', before_end='2021-09-20';
var aDet_start ='2021-09-24', aDet_end ='2021-10-02';
// หา S2 ฟ้าใสสำหรับ validation (+ ครึ่งหน้าต่าง SAR รอบวัน S2)
var s2_start='2021-09-24', s2_end='2021-10-25';
var HALF_WINDOW=6;

// ค่าคงที่ ablation
var FIXED_DB=-14.6;   // C0.5 = ค่า threshold คงที่ค่าเดียวทั้ง 3 พายุ (= mean ของ Otsu)
var RATIO_C0=1.25;    // C0 = change detection ratio แบบเดิม
var SEED=42, NPC=1000;
var PASS='DESCENDING', POLAR='VV';

// ============ ส่วนที่ 2 : ฟังก์ชัน + ชั้นข้อมูลร่วม ============
function otsu(h){h=ee.Dictionary(h);var c=ee.Array(h.get('histogram')),m=ee.Array(h.get('bucketMeans'));
var s=m.length().get([0]),t=c.reduce(ee.Reducer.sum(),[0]).get([0]),su=m.multiply(c).reduce(ee.Reducer.sum(),[0]).get([0]),mn=su.divide(t);
var idx=ee.List.sequence(1,s);var bss=idx.map(function(i){var ac=c.slice(0,0,i);var aC=ac.reduce(ee.Reducer.sum(),[0]).get([0]);
var aM=m.slice(0,0,i).multiply(ac).reduce(ee.Reducer.sum(),[0]).get([0]).divide(aC);var bC=t.subtract(aC);var bM=su.subtract(aC.multiply(aM)).divide(bC);
return aC.multiply(aM.subtract(mn).pow(2)).add(bC.multiply(bM.subtract(mn).pow(2)));});return m.sort(ee.Array(bss)).get([-1]);}

var hand   = ee.Image('MERIT/Hydro/v1_0_1').select('hnd');
var slope  = ee.Algorithms.Terrain(ee.Image('WWF/HydroSHEDS/03VFDEM')).select('slope');
var swater = ee.Image('JRC/GSW1_4/GlobalSurfaceWater').select('seasonality');
var permW  = swater.gte(10);
var permWmask = swater.gte(10).unmask(0);  // ★ 1=น้ำถาวร / 0=อื่น ๆ (valid ทุก pixel — mask ไม่ลาม)

// โหลด S1 (dB) + speckle 50 m
function s1(d0,d1){
  return ee.ImageCollection('COPERNICUS/S1_GRD')
    .filter(ee.Filter.eq('instrumentMode','IW'))
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation',POLAR))
    .filter(ee.Filter.eq('orbitProperties_pass',PASS))
    .filter(ee.Filter.eq('resolution_meters',10))
    .select(POLAR).filterBounds(aoi)
    .filterDate(d0,d1).mosaic().clip(aoi).focal_mean(50,'circle','meters');
}
// Otsu threshold (สุ่ม histogram บนที่ราบลุ่ม HAND<10 & slope<3 ; คิดบน "aoi")
function otsuThr(img){
  var low=img.updateMask(hand.lt(10)).updateMask(slope.lt(3));
  return ee.Number(otsu(low.reduceRegion({reducer:ee.Reducer.histogram(255),
    geometry:aoi, scale:30, maxPixels:1e13, bestEffort:true}).get(POLAR)));
}

// ============ ส่วนที่ 3 : DETECTION → flood extent + Otsu + พื้นที่ (Table 4) ============
var afterDet = s1(aDet_start, aDet_end);          // ภาพช่วงพายุ (peak)
var thrDet   = otsuThr(afterDet);                 // Otsu (จดลง Table 4)
var flood = afterDet.lt(thrDet).selfMask();       // total-water
flood = flood.where(permW,0).selfMask();          // หักน้ำถาวร
flood = flood.updateMask(flood.connectedPixelCount(25).gte(8))  // ตัด noise จุดเล็ก
             .updateMask(slope.lt(5)).updateMask(hand.lt(15));  // ตัดที่ลาดชัน/สูง
var floodBin = flood.unmask(0);

// พื้นที่ (นับบนลุ่มน้ำ)
var ha = ee.Number(flood.multiply(ee.Image.pixelArea()).reduceRegion({
  reducer:ee.Reducer.sum(), geometry:ws, scale:10, maxPixels:1e13, bestEffort:true}).get(POLAR)).divide(10000);
print('======== '+EVENT+' : DETECTION (Table 4) ========');
print('Otsu threshold (dB) =', thrDet);
print('พื้นที่น้ำท่วม (ha)  =', ha);
print('พื้นที่น้ำท่วม (rai) =', ha.multiply(6.25));

Map.centerObject(ws,10);
Map.addLayer(afterDet,{min:-25,max:0},'SAR หลัง (VV dB)',false);
Map.addLayer(flood,{palette:['#0000ff']},EVENT+' flood',true);

Export.image.toDrive({image:floodBin.toByte(), description:'LPP_flood_'+EVENT,
  folder:'GEE_LPP_flood', region:aoi, scale:10, maxPixels:1e13});

// ============ ส่วนที่ 4 : VALIDATION + ABLATION (Table 7/8 + Figure 6) ============
// 4.1 หา S2 ฟ้าใสที่สุด → กำหนดหน้าต่าง SAR สำหรับ validation
function maskS2(img){var scl=img.select('SCL');
  return img.updateMask(scl.neq(3).and(scl.neq(8)).and(scl.neq(9)).and(scl.neq(10)).and(scl.neq(11)));}
var s2col=ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED').filterBounds(aoi)
  .filterDate(s2_start,s2_end).filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE',60))
  .sort('CLOUDY_PIXEL_PERCENTAGE');
var s2best=ee.Image(s2col.first());
var s2date=ee.Date(s2best.get('system:time_start'));
print('======== '+EVENT+' : VALIDATION ========');
print('S2 อ้างอิง วันที่ =', s2date.format('YYYY-MM-dd'),' เมฆ% =', s2best.get('CLOUDY_PIXEL_PERCENTAGE'));

// 4.2 ภาพ SAR รอบวัน S2 (ใช้กับ validation) + ภาพก่อนพายุ (ใช้ทำ C0)
var afterVal = s1(s2date.advance(-HALF_WINDOW,'day'), s2date.advance(HALF_WINDOW,'day'));
var beforeF  = s1(before_start, before_end);
var thrVal   = otsuThr(afterVal);

// 4.3 สร้าง 4 configs จาก afterVal/beforeF/thrVal เดียวกัน
var dbGap = ee.Number(10).multiply(ee.Number(RATIO_C0).log10());   // ratio 1.25 -> ~0.969 dB
var c0  = beforeF.subtract(afterVal).gt(dbGap).and(permWmask.not()).unmask(0).rename('c0');  // เดิม
var c05 = afterVal.lt(FIXED_DB).and(permWmask.not()).unmask(0).rename('c05');                // total-water คงที่
var c1  = afterVal.lt(thrVal).and(permWmask.not()).unmask(0).rename('c1');                   // Otsu ต่อพายุ
var c2m = afterVal.lt(thrVal).and(permWmask.not()).selfMask();                                        // C1 + masking = product
c2m = c2m.updateMask(c2m.connectedPixelCount(25).gte(8)).updateMask(slope.lt(5)).updateMask(hand.lt(15));
var c2  = c2m.unmask(0).rename('c2');

// 4.4 reference จาก S2 (MNDWI>0, หักน้ำถาวร, จำกัดที่ราบลุ่ม)
var s2img=maskS2(s2best).clip(aoi);
var ref=s2img.normalizedDifference(['B3','B11']).gt(0).and(permWmask.not()).updateMask(hand.lt(20)).rename('ref');
print('DEBUG: ref นับ pixel แต่ละ class =', ref.reduceRegion({reducer:ee.Reducer.frequencyHistogram(), geometry:aoi, scale:100, maxPixels:1e13, bestEffort:true}));

// 4.5 สุ่มจุดชุดเดียว (seed คงที่) บนลุ่มน้ำ → วัดทุก config
var stack=ref.addBands([c0,c05,c1,c2]);
var samples=stack.stratifiedSample({numPoints:NPC, classBand:'ref',
  classValues:[0,1], classPoints:[NPC,NPC], region:aoi, scale:10, seed:SEED,
  dropNulls:true, tileScale:4});
print('จำนวนจุดตัวอย่าง =', samples.size());
print('Otsu(validation image, dB) =', thrVal);

function wilson(p){var z=1.96,n=2*NPC,zz=z*z;p=ee.Number(p);
  var d=ee.Number(1).add(zz/n);var cc=p.add(zz/(2*n)).divide(d);
  var hh=p.multiply(ee.Number(1).subtract(p)).divide(n).add(zz/(4*n*n)).sqrt().multiply(z).divide(d);
  return ee.Dictionary({lo:cc.subtract(hh),hi:cc.add(hh)});}

['c0','c05','c1','c2'].forEach(function(k){
  var cm=samples.errorMatrix('ref',k);
  var oa=cm.accuracy(), kp=cm.kappa();
  var pa=ee.Number(cm.producersAccuracy().get([1,0]));   // flood = class 1
  var ua=ee.Number(cm.consumersAccuracy().get([0,1]));
  var f1=pa.multiply(ua).multiply(2).divide(pa.add(ua));
  var ci=wilson(oa);
  print(EVENT+' | '+k, ee.Dictionary({OA:oa,OA_CI_lo:ci.get('lo'),OA_CI_hi:ci.get('hi'),
    Kappa:kp,PA_flood:pa,UA_flood:ua,F1_flood:f1}));
});
// C2 confusion สำหรับ regen Figure 6 -> [[TN,FP],[FN,TP]]
print(EVENT+' | C2 confusion [[TN,FP],[FN,TP]] =', samples.errorMatrix('ref','c2').array());

// 4.6 agreement map (เขียว=ตรง แดง=จับเกิน ส้ม=จับขาด) + export
var refB=ref.unmask(0);
var agree=ee.Image(0).where(c2.eq(1).and(refB.eq(1)),1).where(c2.eq(1).and(refB.eq(0)),2)
  .where(c2.eq(0).and(refB.eq(1)),3).updateMask(ref.mask()).selfMask().rename('agreement');
Map.addLayer(agree,{min:1,max:3,palette:['#2ca25f','#de2d26','#fec44f']},'Agreement',false);
Export.image.toDrive({image:agree.unmask(0).toByte(), description:'LPP_agreement_'+EVENT,
  folder:'GEE_LPP_flood', region:aoi, scale:10, maxPixels:1e13});
