/*******************************************************************************
 * Copyright (C) 2022, Duderstadt Lab
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 ******************************************************************************/
//written by: Karl E. Duderstadt

//This script accompanies the 'FRET dataset analysis using Mars' example workflows as described on the mars docs:
//https://duderstadt-lab.github.io/mars-docs/examples/Static_FRET

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
