import ee

# This class will return a filtered collection
# the filter is the proposed by Quegan&Yu (2001),
# and uses a custom filter to find <J>
# Warning: will only work on natural numbers (no DB).
class QueganYuFilter:

    def __init__(self, collection: ee.imagecollection.ImageCollection, kernel_size) -> None:
        self.collection = collection
        self.kernel_size = kernel_size


    def _medianFilter(self, image: ee.image):
        kernelWindow = ee.Kernel.square((self.kernel_size-1)/2, 'pixels', False)
        return ee.Image(image.reduceNeighborhood(ee.Reducer.median(), kernelWindow)\
            .copyProperties(image,['system:time_start','sliceNumber']))            


    def apply(self):

        collectionMedian = self.collection.map(self._medianFilter)
        
        correctionFactorCol = self.collection.map(
            lambda img: img.divide(self._medianFilter(img))
        )

        correctionFactor = correctionFactorCol.sum().divide(self.collection.count())


        return collectionMedian.map(
            lambda img: ee.Image(img.multiply(correctionFactor).copyProperties(img,['system:time_start'])).select(
                ["VHg0_median", "VVg0_median", "LIA_median"], ["VHg0", "VVg0", "LIA"]
            )
        )