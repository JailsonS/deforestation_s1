import ee

class LeeFilterRefined:


    def apply(image):

        # by Guido Lemoine
        # img must be in natural units, i.e. not in dB!
        img = ee.Image(image)

        # Set up 3x3 kernels 
        weights3 = ee.List.repeat(ee.List.repeat(1,3), 3)
        kernel3 = ee.Kernel.fixed(3,3, weights3, 1, 1, False)

        mean3 = img.reduceNeighborhood(ee.Reducer.mean(), kernel3)
        variance3 = img.reduceNeighborhood(ee.Reducer.variance(), kernel3)

        # Use a sample of the 3x3 windows inside a 7x7 windows to determine gradients and directions
        sampleWeights = ee.List([
            [0,0,0,0,0,0,0], 
            [0,1,0,1,0,1,0],
            [0,0,0,0,0,0,0], 
            [0,1,0,1,0,1,0], 
            [0,0,0,0,0,0,0], 
            [0,1,0,1,0,1,0],
            [0,0,0,0,0,0,0]
        ])

        sampleKernel = ee.Kernel.fixed(7, 7, sampleWeights, 3, 3, False)

        # Calculate mean and variance for the sampled windows and store as 9 bands
        sampleMean = mean3.neighborhoodToBands(sampleKernel)
        sampleVar = variance3.neighborhoodToBands(sampleKernel)


        # Determine the 4 gradients for the sampled windows
        gradients = sampleMean.select(1).subtract(sampleMean.select(7)).abs()
        gradients = gradients.addBands(sampleMean.select(6).subtract(sampleMean.select(2)).abs())
        gradients = gradients.addBands(sampleMean.select(3).subtract(sampleMean.select(5)).abs())
        gradients = gradients.addBands(sampleMean.select(0).subtract(sampleMean.select(8)).abs())

        # And find the maximum gradient amongst gradient bands
        maxGradient = gradients.reduce(ee.Reducer.max())
        
        # Create a mask for band pixels that are the maximum gradient
        gradmask = gradients.eq(maxGradient)
        
        # duplicate gradmask bands: each gradient represents 2 directions
        gradmask = gradmask.addBands(gradmask)
        
    
        # Determine the 8 directions
        directions = sampleMean.select(1).subtract(sampleMean.select(4)).gt(sampleMean.select(4).subtract(sampleMean.select(7))).multiply(1)
        directions = directions.addBands(sampleMean.select(6).subtract(sampleMean.select(4)).gt(sampleMean.select(4).subtract(sampleMean.select(2))).multiply(2))
        directions = directions.addBands(sampleMean.select(3).subtract(sampleMean.select(4)).gt(sampleMean.select(4).subtract(sampleMean.select(5))).multiply(3))
        directions = directions.addBands(sampleMean.select(0).subtract(sampleMean.select(4)).gt(sampleMean.select(4).subtract(sampleMean.select(8))).multiply(4))
        
        # The next 4 are the not() of the previous 4
        directions = directions.addBands(directions.select(0).Not().multiply(5))
        directions = directions.addBands(directions.select(1).Not().multiply(6))
        directions = directions.addBands(directions.select(2).Not().multiply(7))
        directions = directions.addBands(directions.select(3).Not().multiply(8))
        
        # Mask all values that are not 1-8
        directions = directions.updateMask(gradmask)
        
        # "collapse" the stack into a singe band image (due to masking, each pixel has just one value (1-8) in it's directional band, and is otherwise masked)
        directions = directions.reduce(ee.Reducer.sum())

                
        sampleStats = sampleVar.divide(sampleMean.multiply(sampleMean))
        
        # Calculate localNoiseVariance
        sigmaV = sampleStats.toArray().arraySort().arraySlice(0,0,5).arrayReduce(ee.Reducer.mean(), [0])
        
        # Set up the 7*7 kernels for directional statistics
        rectWeights = ee.List.repeat(ee.List.repeat(0,7),3).cat(ee.List.repeat(ee.List.repeat(1,7),4))
        
        diagWeights = ee.List([
            [1,0,0,0,0,0,0], 
            [1,1,0,0,0,0,0], 
            [1,1,1,0,0,0,0], 
            [1,1,1,1,0,0,0], 
            [1,1,1,1,1,0,0], 
            [1,1,1,1,1,1,0], 
            [1,1,1,1,1,1,1]
        ])

        rectKernel = ee.Kernel.fixed(7,7, rectWeights, 3, 3, False)
        diagKernel = ee.Kernel.fixed(7,7, diagWeights, 3, 3, False)

        # Create stacks for mean and variance using the original kernels. Mask with relevant direction.
        dirMean = img.reduceNeighborhood(ee.Reducer.mean(), rectKernel).updateMask(directions.eq(1))
        dirVar = img.reduceNeighborhood(ee.Reducer.variance(), rectKernel).updateMask(directions.eq(1))
        
        dirMean = dirMean.addBands(img.reduceNeighborhood(ee.Reducer.mean(), diagKernel).updateMask(directions.eq(2)))
        dirVar = dirVar.addBands(img.reduceNeighborhood(ee.Reducer.variance(), diagKernel).updateMask(directions.eq(2)))

        # and add the bands for rotated kernels
        for i in [1,2,3]:
            dirMean = dirMean.addBands(img.reduceNeighborhood(ee.Reducer.mean(), rectKernel.rotate(i)).updateMask(directions.eq(2*i+1)))
            dirVar = dirVar.addBands(img.reduceNeighborhood(ee.Reducer.variance(), rectKernel.rotate(i)).updateMask(directions.eq(2*i+1)))
            dirMean = dirMean.addBands(img.reduceNeighborhood(ee.Reducer.mean(), diagKernel.rotate(i)).updateMask(directions.eq(2*i+2)))
            dirVar = dirVar.addBands(img.reduceNeighborhood(ee.Reducer.variance(), diagKernel.rotate(i)).updateMask(directions.eq(2*i+2)))
        
        # "collapse" the stack into a single band image (due to masking, each pixel has just one value in it's directional band, and is otherwise masked)
        dirMean = dirMean.reduce(ee.Reducer.sum())
        dirVar = dirVar.reduce(ee.Reducer.sum())
        
        # A finally generate the filtered value
        varX = dirVar.subtract(dirMean.multiply(dirMean).multiply(sigmaV)).divide(sigmaV.add(1.0))
        
        b = varX.divide(dirVar)
        
        result = dirMean.add(b.multiply(img.subtract(dirMean)))
        return result.arrayFlatten([['sum']]).rename(img.bandNames()).copyProperties(img,['system:time_start'])
