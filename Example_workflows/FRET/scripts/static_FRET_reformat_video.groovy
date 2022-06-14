#@ DatasetService datasetService
#@ Dataset dataset
#@OUTPUT ImgPlus imgPlus

import net.imglib2.view.Views
import net.imagej.axis.AxisType
import net.imagej.axis.Axes
import net.imagej.ImgPlus

img = dataset.getImgPlus().getImg()

viewOfCh1 = Views.subsample(img, 1, 1, 2)
viewOfCh2 = Views.subsample(Views.interval(img, new long[]{0, 0, 1}, new long[]{511, 511, 999}), 1, 1, 2)
rai = Views.stack(viewOfCh1, viewOfCh2)

axes = new AxisType[]{ Axes.X, Axes.Y, Axes.TIME, Axes.CHANNEL}
imgPlus = new ImgPlus( datasetService.create( rai  ), dataset.getName() + "_withChannels", axes )
