import ee, math

# Calculate True azimuth direction for  the near range image edge
def getDESCCorners(f: ee.image.Image):
    # Get the coords as a transposed array
    coords = ee.Array(f.geometry().coordinates().get(0)).transpose()

    crdLons = ee.List(coords.toList().get(0))
    crdLats = ee.List(coords.toList().get(1))

    minLon = crdLons.sort().get(0)
    maxLon = crdLons.sort().get(-1)
    minLat = crdLats.sort().get(0)
    maxLat = crdLats.sort().get(-1)

    azimuth = ee.Number(crdLons.get(crdLats.indexOf(minLat))).subtract(minLon) \
    .atan2(ee.Number(crdLats.get(crdLons.indexOf(minLon))).subtract(minLat)) \
    .multiply(180.0/math.pi) \
    .add(180.0)

    geom = ee.Geometry.LineString([
        crdLons.get(crdLats.indexOf(maxLat)), 
        maxLat, 
        minLon, 
        crdLats.get(crdLons.indexOf(minLon))
    ])

    return ee.Feature(geom, { 'azimuth': azimuth}).copyProperties(f)

def toDB(img):
    return ee.Image(ee.Image(img).log10().multiply(10.0).copyProperties(img,['system:time_start','sliceNumber']))

