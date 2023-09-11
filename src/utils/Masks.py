import ee, math



class Masks:

    def __init__(self) -> None:
        pass


    def getForestMask(asset) -> ee.image.Image:
        mask =  ee.Image(asset)
        return mask
    
    def getDeforestationMask(asset) -> ee.image.Image:
        mask = ee.Image(asset)
        return mask

    def getMaskBorders(self, img):
        az = ee.Number(self._getDESCCorners(img).get('azimuth')).multiply(math.pi/180)
        total_d = ee.Number(1000)
        
        mask = img.select(0).unmask(0).gt(0).selfMask()
        
        displacement = ee.Image(total_d.multiply(az.sin())).addBands(total_d.multiply(az.cos()))
        
        maskD1 = mask.displace(displacement, 'nearest_neighbor')
        maskD2 = mask.displace(displacement.multiply(-1), 'nearest_neighbor')
        
        
        return ee.Image(img.updateMask(maskD1.And(maskD2)))


    
    def getAngleMask(self, image) -> ee.image.Image:
        angle = image.select(['angle'])
        return image.updateMask(angle.lt(45.23993).And(angle.gt(30.63993)))


    # Calculate True azimuth direction for  the near range image edge
    def _getDESCCorners(self, f: ee.image.Image):
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