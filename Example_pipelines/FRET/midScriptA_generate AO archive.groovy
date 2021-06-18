/*******************************************************************************
 * Copyright (C) 2021, Duderstadt Lab
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
//written by: Nadia M. Huisjes, MSc. (Duderstadt lab)

//This script accompanies the 'FRET dataset analysis using Mars' example pipeline as described on the mars docs.
//https://duderstadt-lab.github.io/mars-docs/examples/FRET

#@ Dataset dataset
#@ ImageJ ij
#@ RoiManager roiManager

import de.mpg.biochem.mars.roi.commands.*
import de.mpg.biochem.mars.molecule.*
import de.mpg.biochem.mars.table.*
import de.mpg.biochem.mars.util.*
import org.scijava.table.*
import de.mpg.biochem.mars.image.commands.*
import de.mpg.biochem.mars.image.*
import java.util.HashMap
import java.util.Map
import java.io.File
import java.lang.*

//Set the name of the path and file that will be generated
def file_name =  "0_AOarchive.yama"
def path_name = "/Users/nadiahuisjes/Desktop/Mars Paper/Analysis scripts/210611 - 1-mid Analysis/Archives/"


//Run the peak finder command for the AO archive
PeakFinderCommand peakFinder = new PeakFinderCommand()

peakFinder.setContext(ij.getContext())
peakFinder.setDataset(dataset)
peakFinder.setRoiManager(roiManager)  
peakFinder.setUseRoi(true)
peakFinder.setRoi(x = 6, y = 9, width = 299, height = 490)
peakFinder.setChannel(0)
peakFinder.setUseDogFiler(true)
peakFinder.setDogFilterRadius(5.5)
peakFinder.setThreshold(30)
peakFinder.setMinimumDistance(10)
peakFinder.setT(2)
peakFinder.setFindNegativePeaks(false)
peakFinder.setGeneratePeakCountTable(false)
peakFinder.setGeneratePeakTable(false)
peakFinder.setAddToRoiManager(true)
peakFinder.setProcessAllFrames(false)
peakFinder.setFitPeaks(false)
peakFinder.setFitRadius(4)
peakFinder.setMinimumRsquared(0)
peakFinder.setIntegrate(false)
peakFinder.setIntegrationInnerRadius(1)
peakFinder.setIntegrationOuterRadius(3)
peakFinder.setVerboseOutput(false)
peakFinder.setRoiType("circle")

peakFinder.run()

//Do the ROI Transformation of the AO peaks
TransformROIsCommand transformROIs = new TransformROIsCommand();

transformROIs.setContext(ij.getContext())
transformROIs.setDataset(dataset)
transformROIs.setRoiManager(roiManager)  
transformROIs.setM00(0.986)
transformROIs.setM01(0)
transformROIs.setM02(268.64)
transformROIs.setM10(0)
transformROIs.setM11(0.993)
transformROIs.setM12(1.170)
transformROIs.setTransformationDirection("Long Wavelength to Short Wavelength")
transformROIs.setColocalize(true)
transformROIs.setChannel(1)
transformROIs.setT(2)
transformROIs.setUseDogFilter(true)
transformROIs.setDogFilterRadius(5.5)
transformROIs.setThreshold(20)
transformROIs.setFilterOriginalRois(true)
transformROIs.setFilterColocalizingRois(true)
transformROIs.setColocalizeSearchRadius(2)

transformROIs.run()

//Run the molecule integrator to make the Archive 
#@OUTPUT MoleculeArchive(label="AOarchive.yama") archive

final MoleculeIntegratorDualCommand moleculeIntegrator = new MoleculeIntegratorDualCommand()

moleculeIntegrator.setDataset(dataset)
moleculeIntegrator.setContext(ij.getContext())
moleculeIntegrator.initialize()
moleculeIntegrator.setRoiManager(roiManager)
moleculeIntegrator.setInnerRadius(3)
moleculeIntegrator.setOuterRadius(5)
moleculeIntegrator.setLONGx0(0)
moleculeIntegrator.setLONGy0(0)
moleculeIntegrator.setLONGWidth(252)
moleculeIntegrator.setLONGHeight(512)
moleculeIntegrator.setSHORTx0(257)
moleculeIntegrator.setSHORTy0(0)
moleculeIntegrator.setSHORTWidth(255)
moleculeIntegrator.setSHORTHeight(512)
moleculeIntegrator.setMicroscope("NatureMethodsTraining")
moleculeIntegrator.setFretShortName("Green")
moleculeIntegrator.setFretLongName("Red")
moleculeIntegrator.setIntegrationChannel(0, "Long")
moleculeIntegrator.setIntegrationChannel(1, "Both")
moleculeIntegrator.setThreads(8)

moleculeIntegrator.run()
archive = moleculeIntegrator.getArchive()

//Tag the archive as AO, save the archive, close the ROI manager
archive.metadata().forEach{metadata ->
	metadata.addTag("AO")
}

archive.saveAs(new File(path_name+file_name))
roiManager.close()
