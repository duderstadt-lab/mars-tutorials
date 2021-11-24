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
//written by: prof. dr. Karl E. Duderstadt (Duderstadt lab)

//This script accompanies the 'FRET dataset analysis using Mars' example pipeline as described on the mars docs.
//https://duderstadt-lab.github.io/mars-docs/examples/FRET


#@ Dataset rr_iMap
#@ Dataset og_iMap
#@ MoleculeArchive archive
#@ OpService opService

import net.imglib2.view.Views

int maxRR = opService.stats().max(rr_iMap).getInteger()
int maxOG = opService.stats().max(og_iMap).getInteger()

def rr_iMap_ra = Views.extendMirrorSingle(rr_iMap.getImgPlus()).randomAccess()
def og_iMap_ra = Views.extendMirrorSingle(og_iMap.getImgPlus()).randomAccess()

archive.molecules().forEach{ molecule ->
	molecule.getTable().rows().forEach{ row ->
		double rr = row.getValue("0")
		int rr_x = (int) row.getValue("0_X")
		int rr_y = (int) row.getValue("0_Y")
		int og_x = (int) row.getValue("1_Green_X")
		int og_y = (int) row.getValue("1_Green_Y")

		double rr_corr = rr*(og_iMap_ra.setPositionAndGet(og_x, og_y).getRealDouble()/maxOG)/(rr_iMap_ra.setPositionAndGet(rr_x, rr_y).getRealDouble()/maxRR)
		row.setValue("0_Profile_Corrected", rr_corr)
	}
}
