import ee, math

class CorrectionLIA:

    MDT_IMAGE = ee.Image('USGS/SRTMGL1_003')

    def __init__(self, image: ee.image, azimuth_edge: ee.feature) -> None:
        self.image = image
        self.azimuth_edge = azimuth_edge

    # Compute local incidence angle by Felix Greifeneder, Guido Lemoine 
    def getLIA(self) -> ee.image:

        angleS1 = self.image.select('angle')
        
        azimuthS1 = ee.Terrain.aspect(angleS1)\
            .reduceRegion(ee.Reducer.mean(), angleS1.get('system:footprint'), 100) \
            .get('aspect')
        
        # This should be some degree off the South direction (180), due to Earth rotation
        trueAzimuth = self.azimuth_edge.get('azimuth')

        # rotationFromSouth = ee.Number(trueAzimuth).subtract(180.0)  
        # azimuthS1Number = ee.Number(azimuthS1).add(rotationFromSouth)

        slopeMDT = ee.Terrain.slope(self.MDT_IMAGE).select('slope')
        aspectMDT = ee.Terrain.aspect(slopeMDT).select('aspect')

        slopeProjected = slopeMDT.multiply(
            ee.Image.constant(trueAzimuth)\
            .subtract(90.0)\
            .subtract(aspectMDT)\
            .multiply(math.pi/180)\
            .cos()
        )

        lia = angleS1\
            .subtract(ee.Image.constant(90)\
            .subtract(ee.Image.constant(90)\
            .subtract(slopeProjected)))\
            .abs()
        
        return lia