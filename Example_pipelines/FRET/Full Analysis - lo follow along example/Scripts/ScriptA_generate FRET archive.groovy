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
def file_name =  "0_FRETarchive.yama"
def path_name = "/Users/nadiahuisjes/Desktop/Mars Paper/Analysis scripts/210608 - 1-lo Analysis/Archives/"


//Run the Peak Finder command for the FRET archive
PeakFinderCommand peakFinder2 = new PeakFinderCommand()

peakFinder2.setContext(ij.getContext())
peakFinder2.setDataset(dataset)
peakFinder2.setRoiManager(roiManager)  
peakFinder2.setUseRoi(true)
peakFinder2.setRoi(x = 6, y = 9, width = 299, height = 490)
peakFinder2.setChannel(0)
peakFinder2.setUseDogFiler(true)
peakFinder2.setDogFilterRadius(5.5)
peakFinder2.setThreshold(500)
peakFinder2.setMinimumDistance(10)
peakFinder2.setT(2)
peakFinder2.setFindNegativePeaks(false)
peakFinder2.setGeneratePeakCountTable(false)
peakFinder2.setGeneratePeakTable(false)
peakFinder2.setAddToRoiManager(true)
peakFinder2.setProcessAllFrames(false)
peakFinder2.setFitPeaks(false)
peakFinder2.setFitRadius(4)
peakFinder2.setMinimumRsquared(0)
peakFinder2.setIntegrate(false)
peakFinder2.setIntegrationInnerRadius(1)
peakFinder2.setIntegrationOuterRadius(3)
peakFinder2.setVerboseOutput(false)
peakFinder2.setRoiType("circle")

peakFinder2.run()


//Do the ROI Transformation of the FRET peaks
TransformROIsCommand transformROIs2 = new TransformROIsCommand();

transformROIs2.setContext(ij.getContext())
transformROIs2.setDataset(dataset)
transformROIs2.setRoiManager(roiManager)  
transformROIs2.setM00(0.986)
transformROIs2.setM01(0)
transformROIs2.setM02(268.64)
transformROIs2.setM10(0)
transformROIs2.setM11(0.993)
transformROIs2.setM12(1.170)
transformROIs2.setTransformationDirection("Long Wavelength to Short Wavelength")
transformROIs2.setColocalize(true)
transformROIs2.setChannel(1)
transformROIs2.setT(2)
transformROIs2.setUseDogFilter(true)
transformROIs2.setDogFilterRadius(5.5)
transformROIs2.setThreshold(500)
transformROIs2.setFilterOriginalRois(true)
transformROIs2.setFilterColocalizingRois(false)
transformROIs2.setColocalizeSearchRadius(2)

transformROIs2.run()

//Run the molecule integrator to make the Archive 
#@OUTPUT MoleculeArchive(label="FRETarchive.yama") archive

final MoleculeIntegratorDualCommand moleculeIntegrator2 = new MoleculeIntegratorDualCommand()

moleculeIntegrator2.setDataset(dataset)
moleculeIntegrator2.setContext(ij.getContext())
moleculeIntegrator2.initialize()
moleculeIntegrator2.setRoiManager(roiManager)
moleculeIntegrator2.setInnerRadius(3)
moleculeIntegrator2.setOuterRadius(5)
moleculeIntegrator2.setLONGx0(0)
moleculeIntegrator2.setLONGy0(0)
moleculeIntegrator2.setLONGWidth(252)
moleculeIntegrator2.setLONGHeight(512)
moleculeIntegrator2.setSHORTx0(257)
moleculeIntegrator2.setSHORTy0(0)
moleculeIntegrator2.setSHORTWidth(255)
moleculeIntegrator2.setSHORTHeight(512)
moleculeIntegrator2.setMicroscope("NatureMethodsTraining")
moleculeIntegrator2.setFretShortName("Green")
moleculeIntegrator2.setFretLongName("Red")
moleculeIntegrator2.setIntegrationChannel(0, "Long")
moleculeIntegrator2.setIntegrationChannel(1, "Both")
moleculeIntegrator2.setThreads(8)

moleculeIntegrator2.run()
archive = moleculeIntegrator2.getArchive()

//Tag the archive as FRET, save the archive, close the ROI manager
archive.metadata().forEach{metadata ->
	metadata.addTag("FRET")
}


archive.saveAs(new File(path_name+file_name))
roiManager.close()


