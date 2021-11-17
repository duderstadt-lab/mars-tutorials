#@ ImagePlus image
#@output ImagePlus twoChannelImage

import ij.IJ
import ij.ImagePlus

twoChannelImage = IJ.createHyperStack(image.getTitle(), image.getWidth(), image.getHeight(), 2, 1, (int)(image.getNFrames()/2), 16)

newFrameIndex = 1
for (int frame=1; frame < image.getNFrames(); frame+=2) {
	image.setPosition(frame)
	twoChannelImage.setPosition(1, 1, newFrameIndex)
	for (int x=0; x < image.getWidth(); x++)
		for (int y=0; y < image.getHeight(); y++)
			 twoChannelImage.getProcessor().set(x, y, image.getProcessor().get(x, y))
	
	image.setPosition(frame+1)
	twoChannelImage.setPosition(2, 1, newFrameIndex)
	for (int x=0; x < image.getWidth(); x++)
		for (int y=0; y < image.getHeight(); y++)
			 twoChannelImage.getProcessor().set(x, y, image.getProcessor().get(x, y))
			 
	newFrameIndex++
}